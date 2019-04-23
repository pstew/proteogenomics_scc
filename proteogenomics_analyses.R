library(stringr)
library(gplots)
library(LSD)
library(devtools)
library(ConsensusClusterPlus)
library(VennDiagram)
library(enrichR)
library(multiOmicsViz)
library(biomaRt)
library(ggbiplot)
library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(data.table)
library(survminer)
library(survival)
library(Rtsne)
setwd("~/projects/proteogenomics")
source("~/scripts/data_processing_functions.R") #primarily proteomics data processing functions
source("scripts/main/heatmap.4.R") #slight modification of https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R to make heatmap font size adjustable

#read input ----

protein_expression = read.delim("data/main/original_data/v2.1/tmt29_processed_data_uncorrected.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) 
protein_annotation = read.delim("data/main/original_data/v2.1/tmt29_protein_group_annotation.txt", header = TRUE, stringsAsFactors = FALSE) 
protein_expression = merge(protein_annotation, protein_expression, by = "ProteinGroup")
protein_expression = subset(protein_expression, All_Reverse_Flag == 0)
tmt_meta = read.delim("data/main/original_data/v2.1/tmt29_sample_annotation.txt", header = TRUE, stringsAsFactors = FALSE)

#need to take out the following samples for treatment prior to surgery; very last sample ...7303 has problem with DNA so skip
#these 06s identifiers are internal Moffitt identifers. Samples have been renamed in supplemental data for the manuscript, e.g. SCC001, SCC002, etc.
skip_samples = c("06S15101140", "06S15101158", "06S15101159", "06S15101176", "06S15101180", "06S15101185", "06S15101193", "06S10347307")
tmt_meta = subset(tmt_meta, !(TissueID %in% skip_samples))

#Get rid of pools; going to call the sample names directly from the metadata column
tumor_tmt_colnames = subset(tmt_meta, TissueID != "Pool")$TissueID 
protein_expression[protein_expression == 0] = NA #no imputation except in select cases (e.g. PCA)

wilkerson_genes_centroids = read.delim("data/main/wilkerson.scc/predictor.centroids.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE) #used for getting gene symbols for looking at specific Wilkerson SCC signature genes and their correlations
gene_expression = read.delim("data/main/proteogenomics_cel_file_expression_sampleidmap_08-15-16.txt", header = TRUE, stringsAsFactors = FALSE) #microarray gene expression
gene_expression_best = read.delim("data/main/gene_expression_best_03-24-17.txt", header = TRUE, stringsAsFactors = FALSE) #already log2'd; processed with get_best_probe() from data_processing_functions.R
gene_expression_06s_mapping = read.delim("data/main/proteogenomics_cel_file_meta_03-15-17.txt", header = TRUE, stringsAsFactors = FALSE)
mutation_status_prelim = read.delim("data/main/preliminary_mutations_12-20-16/LUSC_PG_genes_all.list_minGQ15.status_mapped_sample.txt", header = TRUE, stringsAsFactors = FALSE) #this is called "preliminary" but hasn't changed and was used in the final analysis for the manuscript
mutation_gene_list = names(mutation_status_prelim)[names(mutation_status_prelim) != "tumor_06s" & names(mutation_status_prelim) != "dna_sample"] 
tumor_mutation_burden = read.delim("data/main/tmb_07-21-17.txt", header = TRUE, stringsAsFactors = FALSE) # used to look if there is any association with proteomic subtype and mutation burden
mutation_status_prelim_tmb = merge(mutation_status_prelim, tumor_mutation_burden, by.x = "dna_sample", by.y = "sample")
sqlc_classification_combined = read.delim("data/main/squamous_subtype_classifications_microarray_rnaseq_2018-05-10.txt", header = TRUE, stringsAsFactors = FALSE)
clinical = read.delim("data/main/sqlc116_clinical_05-23-17.txt", header = TRUE, stringsAsFactors = FALSE)
rnaseq_meta = read.delim("data/main/v2.2_lusq_proteogenomics_rnaseq_metadata_with_batches_2018-05-09.txt", header = TRUE, stringsAsFactors = FALSE) #updated metadata for remapped batch 3
rnaseq_meta = subset(rnaseq_meta, Redo != "Original" & Sample.Family.ID != "---") #don't include original samples that needed to be rerun
smoking = read.delim("data/main/smoking_trimmed_merged_06-16-17.txt", header = TRUE, stringsAsFactors = FALSE)
lsqcc_tcc_signatures = read.delim("data/main/original_data/LSQCC-signatures_TCC_from_steven_10-04-17.txt", header = TRUE, stringsAsFactors = FALSE)

#scaling factors used for multi-plex normalization
scaling_factors = read.delim("data/main/original_data/v2.1/tmt29_replicate_scaling_factors_uncorrected..txt", header = TRUE, stringsAsFactors = FALSE)

scaling_factors = scaling_factors[str_detect(scaling_factors$Plex, "01$"), ]
scaling_factors$Plex = str_replace(scaling_factors$Plex, "\\_01$", "") #use first technical replicate 

scaling_factors$TMT.Sample.Name = paste(scaling_factors$Plex, scaling_factors$SampleName, sep = "_")
scaling_factors = subset(scaling_factors, SampleName != "TMT-126")


#combine into master metadata file to line up data by sample and for easy reference
tmt_meta = merge(tmt_meta, gene_expression_06s_mapping, by = "TissueID", all = TRUE)
tmt_meta = merge(tmt_meta, mutation_status_prelim_tmb, by.x = "TissueID", by.y = "tumor_06s", all = TRUE)
tmt_meta = merge(tmt_meta, sqlc_classification_combined, by.x = "TissueID", by.y = "tumor_06s", all = TRUE)
tmt_meta = merge(tmt_meta, clinical, by.x = "TissueID", by.y = "Aliquot.Sample.ID")
tmt_meta = merge(tmt_meta, rnaseq_meta, by = "TissueID")
tmt_meta = merge(tmt_meta, smoking, by.x = "TissueID", by.y = "Aliquot.ID")
tmt_meta = merge(tmt_meta, lsqcc_tcc_signatures, by.x = "TissueID", by.y = "row.names", all = TRUE)
tmt_meta = merge(tmt_meta, scaling_factors, by = "TMT.Sample.Name")

number_of_samples = nrow(tmt_meta) #use this for calculations needing a percentage of the total number of samples


rnaseq_expression = read.delim("data/main/original_data/v2.2_lusq116_rnaseq_iron_untilt_filt3_debatched.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

#get TMT columns, gene columns, and patients; order is the same since these are coming from the meta file
tumor_06s = subset(tmt_meta, TissueID != "Pool")$TissueID


#rnaseq consensus cluster ----
rnaseq_colnames = tmt_meta$GenomicsSample
rnaseq_expression$sd = apply(rnaseq_expression[, rnaseq_colnames], 1, function(x) {
	sd(na.omit(x))
})
rnaseq_expression$mean = apply(rnaseq_expression[, rnaseq_colnames], 1, function(x) {
	mean(na.omit(x))
})
rnaseq_expression$mads = apply(rnaseq_expression[, rnaseq_colnames], 1, function(x) {
	mad(na.omit(x))
})
rnaseq_expression$na_count = rowSums(is.na(rnaseq_expression[, rnaseq_colnames]))


rnaseq_expression_heatmap = subset(rnaseq_expression, na_count < number_of_samples * 0.1 ) #less than 10% missingness

#get top 1000 variable genes
rnaseq_expression_heatmap = rnaseq_expression_heatmap[rev(order(rnaseq_expression_heatmap$mads))[1:1000], rnaseq_colnames]

title_rnaseq = tempdir()
#title_rnaseq = "output/main/cc_results/rnaseq_cc_complete_madtop1k_10pct_missing_1krep_2018-05-10"
consensus_rnaseq_matrix = as.matrix(rnaseq_expression_heatmap[, rnaseq_colnames])
consensus_rnaseq_matrix = sweep(consensus_rnaseq_matrix, 1, apply(consensus_rnaseq_matrix, 1, median, na.rm = T))
rnaseq_results = ConsensusClusterPlus(consensus_rnaseq_matrix, maxK = 13, reps = 1000, pItem = 0.8, pFeature = 1, title = title_rnaseq, clusterAlg = "hc", distance = "pearson", seed = 1234, plot = "png", innerLinkage = "complete", finalLinkage = "complete", corUse = "pairwise.complete.obs")

#extract consensus clustering 
tmt_meta = merge(tmt_meta, as.data.frame(rnaseq_results[[5]]["consensusClass"]), by.x = "GenomicsSample", by.y = "row.names")
names(tmt_meta)[names(tmt_meta) == "consensusClass"] = "consensusClassK5_rnaseq"

#coloring to mirror the cc plots
setcolors = brewer.pal(5, "Set2")

k5_rnaseq_cc_group1_color = setcolors[1]
k5_rnaseq_cc_group2_color = setcolors[2]
k5_rnaseq_cc_group3_color = setcolors[3]
k5_rnaseq_cc_group4_color = setcolors[4]
k5_rnaseq_cc_group5_color = setcolors[5]

#order RNAseq meta by cc clustering
tmt_meta = tmt_meta[order(tmt_meta$consensusClassK5_rnaseq), ] #order according to consensus class for displaying mutations and classifications
gene_mutations_for_heatmap = c("TP53", "MLL2", "NFE2L2", "CDKN2A", "NOTCH1", "APC", "KEAP1", "PIK3CA", "PTEN", "RB1")
k5_rnaseq_coloring = tmt_meta[, rev(c(gene_mutations_for_heatmap, "Classification_rnaseq", "consensusClassK5_rnaseq"))]


k5_rnaseq_coloring$consensusClassK5_rnaseq = str_replace(k5_rnaseq_coloring$consensusClassK5_rnaseq, "^1$", fixed(k5_rnaseq_cc_group1_color))
k5_rnaseq_coloring$consensusClassK5_rnaseq = str_replace(k5_rnaseq_coloring$consensusClassK5_rnaseq, "^2$", fixed(k5_rnaseq_cc_group2_color))
k5_rnaseq_coloring$consensusClassK5_rnaseq = str_replace(k5_rnaseq_coloring$consensusClassK5_rnaseq, "^3$", fixed(k5_rnaseq_cc_group3_color))
k5_rnaseq_coloring$consensusClassK5_rnaseq = str_replace(k5_rnaseq_coloring$consensusClassK5_rnaseq, "^4$", fixed(k5_rnaseq_cc_group4_color))
k5_rnaseq_coloring$consensusClassK5_rnaseq = str_replace(k5_rnaseq_coloring$consensusClassK5_rnaseq, "^5$", fixed(k5_rnaseq_cc_group5_color))
k5_rnaseq_coloring$Classification_rnaseq = str_replace(k5_rnaseq_coloring$Classification_rnaseq, "classical", "red")
k5_rnaseq_coloring$Classification_rnaseq = str_replace(k5_rnaseq_coloring$Classification_rnaseq, "basal", "blue")
k5_rnaseq_coloring$Classification_rnaseq = str_replace(k5_rnaseq_coloring$Classification_rnaseq, "secretory", "green")
k5_rnaseq_coloring$Classification_rnaseq = str_replace(k5_rnaseq_coloring$Classification_rnaseq, "primitive", "black")
k5_rnaseq_coloring[, gene_mutations_for_heatmap] = apply(k5_rnaseq_coloring[, gene_mutations_for_heatmap], 2, function(x) { str_replace(x, "mut", "aquamarine4")})
k5_rnaseq_coloring[, gene_mutations_for_heatmap] = apply(k5_rnaseq_coloring[, gene_mutations_for_heatmap], 2, function(x) { str_replace(x, "wt", "lightgrey")})
k5_rnaseq_coloring[, gene_mutations_for_heatmap] = apply(k5_rnaseq_coloring[, gene_mutations_for_heatmap], 2, function(x) { str_replace(x, "trunc", "aquamarine")})

colnames(k5_rnaseq_coloring)[colnames(k5_rnaseq_coloring) == "consensusClassK5_rnaseq"] = "RNAseq"
colnames(k5_rnaseq_coloring)[colnames(k5_rnaseq_coloring) == "Classification_rnaseq"] = "Wilkerson"

heatmap_obj = heatmap.4(rnaseq_expression_heatmap[, tmt_meta$GenomicsSample], Colv = FALSE, col = bluered, trace = "none", density.info = "none", main = "", scale = "row", dendrogram = "none", labRow = FALSE, labCol = FALSE, ColSideColors = as.matrix(k5_rnaseq_coloring), margins = c(2,2), ColSideColorsSize = 16, xlab = "Patients", ylab = "Genes", na.color = "darkgrey", breaks = c(-2, seq(-2, 0, length = 50), seq(0, 2, length = 50)), distfun = function(c){as.dist(1 - cor(t(c), use = "pairwise.complete.obs"))})
legend("bottomleft", legend = c("WT", "Mut", "Trunc", "Classical", "Basal", "Secretory", "Primitive"), fill = c("lightgrey", "aquamarine4", "aquamarine", "red", "blue", "green", "black"), cex = 1)

#protein consensus cluster ----
protein_expression$sd = apply(protein_expression[, tumor_tmt_colnames], 1, function(x) {
	sd(na.omit(x))
})
protein_expression$mean = apply(protein_expression[, tumor_tmt_colnames], 1, function(x) {
	mean(na.omit(x))
})
protein_expression$mads = apply(protein_expression[, tumor_tmt_colnames], 1, function(x) {
	mad(na.omit(x))
})
protein_expression$na_count = rowSums(is.na(protein_expression[, tumor_tmt_colnames]))

protein_expression_10pct = subset(protein_expression, na_count < number_of_samples * 0.1 ) #10% missingness
protein_expression_heatmap = protein_expression_10pct[rev(order(protein_expression_10pct$mads))[1:1000], tumor_tmt_colnames]

row.names(protein_expression_heatmap) = protein_expression_10pct[rev(order(protein_expression_10pct$mads))[1:1000], "Symbol"] 
protein_expression_heatmap = protein_expression_10pct[rev(order(protein_expression_10pct$mads))[1:1000], c("Accession", "Symbol", tumor_tmt_colnames)] 


title_protein = tempdir()

consensus_protein_matrix = as.matrix(protein_expression_heatmap[, tumor_tmt_colnames])
consensus_protein_matrix = sweep(consensus_protein_matrix, 1, apply(consensus_protein_matrix, 1, median, na.rm= T))
protein_results = ConsensusClusterPlus(consensus_protein_matrix, maxK = 13, reps= 1000, pItem= 0.8, pFeature= 1, title = title_protein, clusterAlg= "hc", distance= "pearson", seed= 1234, plot= "png", innerLinkage = "complete", finalLinkage = "complete", corUse = "pairwise.complete.obs")

#extract consensus clustering into tmt_meta file
tmt_meta = merge(tmt_meta, as.data.frame(protein_results[[4]]["consensusClass"]), by.x = "TissueID", by.y = "row.names")
names(tmt_meta)[names(tmt_meta) == "consensusClass"] = "consensusClassK4_protein"
tmt_meta = merge(tmt_meta, as.data.frame(protein_results[[5]]["consensusClass"]), by.x = "TissueID", by.y = "row.names")
names(tmt_meta)[names(tmt_meta) == "consensusClass"] = "consensusClassK5_protein"

#protein subtype 
tmt_meta$subtype = NA
for (i in 1:nrow(tmt_meta)) {
	if (tmt_meta$consensusClassK5_protein[i] == 1 | tmt_meta$consensusClassK5_protein[i] == 3) {
		tmt_meta$subtype[i] = "Inflamed" #immune lump group
	} 
	if (tmt_meta$consensusClassK5_protein[i] == 2 | tmt_meta$consensusClassK5_protein[i] == 4) {
		tmt_meta$subtype[i] = "Redox" #metabolism lump group
	} 
	if (tmt_meta$consensusClassK5_protein[i] == 5) {
		tmt_meta$subtype[i] = "Mixed" #mixed lump group
	} 
}

#protein subtype numeric
tmt_meta$subtype_numeric = NA
for (i in 1:nrow(tmt_meta)) {
	if (tmt_meta$consensusClassK5_protein[i] == 1 | tmt_meta$consensusClassK5_protein[i] == 3) {
		tmt_meta$subtype_numeric[i] = 1 #immune lump group
	} 
	if (tmt_meta$consensusClassK5_protein[i] == 2 | tmt_meta$consensusClassK5_protein[i] == 4) {
		tmt_meta$subtype_numeric[i] = 2 #metabolism lump group
	} 
	if (tmt_meta$consensusClassK5_protein[i] == 5) {
		tmt_meta$subtype_numeric[i] = 3 #mixed lump group
	} 
}


#protein heatmaps ----


k4_protein_cc_group1_color = "#33A22B"
k4_protein_cc_group2_color = "#B1DF8A"
k4_protein_cc_group3_color = "#A5CEE2"
k4_protein_cc_group4_color = "#F89998"

darkcolors = brewer.pal(5, "Dark2")

k5_protein_cc_group1_color = darkcolors[1]
k5_protein_cc_group3_color = darkcolors[1]
k5_protein_cc_group2_color = darkcolors[2]
k5_protein_cc_group4_color = darkcolors[2]
k5_protein_cc_group5_color = darkcolors[3]


#re-cluster groups defined by consesnsus clustering for visualizing expression in the heatmap
protein_expression_heatmap_group1 = scale(protein_expression_heatmap[, subset(tmt_meta, consensusClassK5_protein == 1)$TissueID], center = TRUE)
protein_expression_heatmap_group1_dist = as.dist(1 - cor(protein_expression_heatmap_group1, use = "pairwise.complete.obs"))
protein_expression_heatmap_group1_clust = hclust(protein_expression_heatmap_group1_dist, method = "complete")
protein_expression_heatmap_group1_order = colnames(protein_expression_heatmap_group1)[protein_expression_heatmap_group1_clust$order]
protein_expression_heatmap_group2 = scale(protein_expression_heatmap[, subset(tmt_meta, consensusClassK5_protein == 2)$TissueID], center = TRUE)
protein_expression_heatmap_group2_dist = as.dist(1 - cor(protein_expression_heatmap_group2, use = "pairwise.complete.obs"))
protein_expression_heatmap_group2_clust = hclust(protein_expression_heatmap_group2_dist, method = "complete")
protein_expression_heatmap_group2_order = colnames(protein_expression_heatmap_group2)[protein_expression_heatmap_group2_clust$order]
protein_expression_heatmap_group3 = scale(protein_expression_heatmap[, subset(tmt_meta, consensusClassK5_protein == 3)$TissueID], center = TRUE)
protein_expression_heatmap_group3_dist = as.dist(1 - cor(protein_expression_heatmap_group3, use = "pairwise.complete.obs"))
protein_expression_heatmap_group3_clust = hclust(protein_expression_heatmap_group3_dist, method = "complete")
protein_expression_heatmap_group3_order = colnames(protein_expression_heatmap_group3)[protein_expression_heatmap_group3_clust$order]
protein_expression_heatmap_group4 = scale(protein_expression_heatmap[, subset(tmt_meta, consensusClassK5_protein == 4)$TissueID], center = TRUE)
protein_expression_heatmap_group4_dist = as.dist(1 - cor(protein_expression_heatmap_group4, use = "pairwise.complete.obs"))
protein_expression_heatmap_group4_clust = hclust(protein_expression_heatmap_group4_dist, method = "complete")
protein_expression_heatmap_group4_order = colnames(protein_expression_heatmap_group4)[protein_expression_heatmap_group4_clust$order]
protein_expression_heatmap_group5 = scale(protein_expression_heatmap[, subset(tmt_meta, consensusClassK5_protein == 5)$TissueID], center = TRUE)
protein_expression_heatmap_group5_dist = as.dist(1 - cor(protein_expression_heatmap_group5, use = "pairwise.complete.obs"))
protein_expression_heatmap_group5_clust = hclust(protein_expression_heatmap_group5_dist, method = "complete")
protein_expression_heatmap_group5_order = colnames(protein_expression_heatmap_group5)[protein_expression_heatmap_group5_clust$order]
final_heatmap_ordering = c(protein_expression_heatmap_group1_order, protein_expression_heatmap_group3_order, protein_expression_heatmap_group2_order, protein_expression_heatmap_group4_order, protein_expression_heatmap_group5_order)

tmt_meta = tmt_meta[match(final_heatmap_ordering, tmt_meta$TissueID), ] #match the final ordering to the metadata


#matrix 
k5_protein_coloring = tmt_meta[, rev(c(gene_mutations_for_heatmap, "Classification_rnaseq", "consensusClassK5_rnaseq", "consensusClassK5_protein"))]


k5_protein_coloring$consensusClassK5_rnaseq = str_replace(k5_protein_coloring$consensusClassK5_rnaseq, "^1$", fixed(k5_rnaseq_cc_group1_color))
k5_protein_coloring$consensusClassK5_rnaseq = str_replace(k5_protein_coloring$consensusClassK5_rnaseq, "^2$", fixed(k5_rnaseq_cc_group2_color))
k5_protein_coloring$consensusClassK5_rnaseq = str_replace(k5_protein_coloring$consensusClassK5_rnaseq, "^3$", fixed(k5_rnaseq_cc_group3_color))
k5_protein_coloring$consensusClassK5_rnaseq = str_replace(k5_protein_coloring$consensusClassK5_rnaseq, "^4$", fixed(k5_rnaseq_cc_group4_color))
k5_protein_coloring$consensusClassK5_rnaseq = str_replace(k5_protein_coloring$consensusClassK5_rnaseq, "^5$", fixed(k5_rnaseq_cc_group5_color))

k5_protein_coloring$Classification_rnaseq = str_replace(k5_protein_coloring$Classification_rnaseq, "classical", "red")
k5_protein_coloring$Classification_rnaseq = str_replace(k5_protein_coloring$Classification_rnaseq, "basal", "blue")
k5_protein_coloring$Classification_rnaseq = str_replace(k5_protein_coloring$Classification_rnaseq, "secretory", "green")
k5_protein_coloring$Classification_rnaseq = str_replace(k5_protein_coloring$Classification_rnaseq, "primitive", "black")

k5_protein_coloring[, gene_mutations_for_heatmap] = apply(k5_protein_coloring[, gene_mutations_for_heatmap], 2, function(x) { str_replace(x, "mut", "#4DAF4A")})
k5_protein_coloring[, gene_mutations_for_heatmap] = apply(k5_protein_coloring[, gene_mutations_for_heatmap], 2, function(x) { str_replace(x, "wt", "lightgrey")})
k5_protein_coloring[, gene_mutations_for_heatmap] = apply(k5_protein_coloring[, gene_mutations_for_heatmap], 2, function(x) { str_replace(x, "trunc", "#984EA3")})

k5_protein_coloring$consensusClassK5_protein = str_replace(k5_protein_coloring$consensusClassK5_protein, "^1$", fixed(k5_protein_cc_group1_color))
k5_protein_coloring$consensusClassK5_protein = str_replace(k5_protein_coloring$consensusClassK5_protein, "^2$", fixed(k5_protein_cc_group2_color))
k5_protein_coloring$consensusClassK5_protein = str_replace(k5_protein_coloring$consensusClassK5_protein, "^3$", fixed(k5_protein_cc_group3_color))
k5_protein_coloring$consensusClassK5_protein = str_replace(k5_protein_coloring$consensusClassK5_protein, "^4$", fixed(k5_protein_cc_group4_color))
k5_protein_coloring$consensusClassK5_protein = str_replace(k5_protein_coloring$consensusClassK5_protein, "^5$", fixed(k5_protein_cc_group5_color))

k5_protein_coloring$CD3 = str_replace(k5_protein_coloring$CD3, "high", "darkred")
k5_protein_coloring$CD3 = str_replace(k5_protein_coloring$CD3, "low", "steelblue")
k5_protein_coloring$CD4 = str_replace(k5_protein_coloring$CD4, "high", "darkred")
k5_protein_coloring$CD4 = str_replace(k5_protein_coloring$CD4, "low", "steelblue")
k5_protein_coloring$CD8 = str_replace(k5_protein_coloring$CD8, "high", "darkred")
k5_protein_coloring$CD8 = str_replace(k5_protein_coloring$CD8, "low", "steelblue")
k5_protein_coloring$PDL1 = str_replace(k5_protein_coloring$PDL1, "high", "darkred")
k5_protein_coloring$PDL1 = str_replace(k5_protein_coloring$PDL1, "low", "steelblue")
k5_protein_coloring$PD1 = str_replace(k5_protein_coloring$PD1, "high", "darkred")
k5_protein_coloring$PD1 = str_replace(k5_protein_coloring$PD1, "low", "steelblue")
k5_protein_coloring$TMB_status = str_replace(k5_protein_coloring$TMB_status, "high", "darkred")
k5_protein_coloring$TMB_status = str_replace(k5_protein_coloring$TMB_status, "low", "steelblue")

colnames(k5_protein_coloring)[colnames(k5_protein_coloring) == "consensusClassK5_rnaseq"] = "RNAseq"
colnames(k5_protein_coloring)[colnames(k5_protein_coloring) == "consensusClassK5_protein"] = "Proteomics Group"
colnames(k5_protein_coloring)[colnames(k5_protein_coloring) == "Classification_rnaseq"] = "Wilkerson et al."


heatmap_obj = heatmap.4(protein_expression_heatmap[, tmt_meta$TissueID], Colv = FALSE, col = bluered, trace = "none", density.info = "none", main = "", scale = "row", dendrogram = "none", labRow = FALSE, labCol = FALSE, ColSideColors = as.matrix(k5_protein_coloring[, c(19:ncol(k5_protein_coloring))]), ColSideColorsSize = 20, cexCol = 1.5, margins = c(2,2), na.color = "darkgrey", breaks = c(-2, seq(-2, 0, length = 50), seq(0, 2, length = 50)), xlab = "Patients", ylab = "Protein Expression",  distfun = function(c){as.dist(1 - cor(t(c), use = "pairwise.complete.obs"))}, colsep = c(23, 43, 55, 94), sepcolor = "black")
legend("bottomleft", legend = c("WT", "Mut", "Trunc", "Classical", "Basal", "Secretory", "Primitive"), fill = c("lightgrey", "#4DAF4A", "#984EA3", "red", "blue", "green", "black"), cex = 1)



#heatmap object for dendrogram needs labeled rows or it defaults to the numeric row numbering from the original protein expression matrix
row.names(protein_expression_heatmap) = protein_expression_heatmap$Accession
heatmap_obj = heatmap.4(protein_expression_heatmap[, tmt_meta$TissueID], Colv = FALSE, col = bluered, trace = "none", density.info = "none", main = "", scale = "row", dendrogram = "none", labRow = protein_expression_heatmap$Accession, labCol = FALSE, ColSideColors = as.matrix(k5_protein_coloring[, c(19:ncol(k5_protein_coloring))]), ColSideColorsSize = 20, cexCol = 1.5, margins = c(2,2), na.color = "darkgrey", breaks = c(-2, seq(-2, 0, length = 50), seq(0, 2, length = 50)), xlab = "Patients", ylab = "Protein Expression",  distfun = function(c){as.dist(1 - cor(t(c), use = "pairwise.complete.obs"))}, colsep = c(23, 43, 55, 94), sepcolor = "black")


#extract dendrogram
heat_hclust = heatmap_obj$rowDendrogram
#plot(heat_hclust, labels = FALSE, axes = FALSE)
#labels_colors(heat_hclust) = c("red", "black", "blue", "green", "darkred")
plot(heat_hclust)
#plot(color_branches(heat_hclust, k = 7))
#heat_hclust_cut = cutree(heat_hclust, k = 7)
#heat_hclust_cut = cutree(heat_hclust, h = 1.4)
cut_at = 1.5
plot(color_branches(heat_hclust, h = cut_at, labels = F))
cutree(heat_hclust, h = cut_at)[heatmap_obj$rowInd]
heat_cut_table = as.data.frame(cutree(heat_hclust, h = cut_at)[heatmap_obj$rowInd])
names(heat_cut_table) = "cluster"
paste(row.names(subset(heat_cut_table, heat_cut_table$cluster == 1)), collapse = ",")
#reversed order of heatmap 
unique(heat_cut_table$cluster)
#1.4
#[1] 5 6 3 2 7 4 1

#1.5
#[1] 4 3 2 5 1
heat_cut_table = merge(heat_cut_table, protein_expression_10pct[, c("Symbol", "Accession")], by.x = "row.names", by.y = "Accession", sort = FALSE)
names(heat_cut_table) = c("Acession", "cluster", "Symbol")
heat_cut_table = heat_cut_table[, c(3, 1, 2)]

subset(heat_cut_table, cluster == 1)$Symbol
#sorted by padj, get top hit
#4 bicarbonate transport (GO:0015701) padj 0.0003474
#3 glutathione metabolic process (GO:0006749) padj 4.053e-10
#2 platelet degranulation (GO:0002576) padj 1.453e-19
#5 extracellular matrix organization (GO:0030198) padj 9.660e-22
#1 neutrophil degranulation (GO:0043312) padj 4.944e-21


#correlation: protein/rnaseq ----

protein_rnaseq_expression_for_correlations = merge(protein_expression_10pct , rnaseq_expression, by.x = "Symbol", by.y = "Symbol")

correlation_table_rnaseq_by_gene = data.frame(Accession = protein_rnaseq_expression_for_correlations$Symbol, spearman_corr = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(protein_rnaseq_expression_for_correlations)) {
	corr_result = tryCatch({
		proteins_for_corr = as.numeric(protein_rnaseq_expression_for_correlations[i, tmt_meta$TissueID])
		genes_for_corr = as.numeric(protein_rnaseq_expression_for_correlations[i, tmt_meta$GenomicsSample])
		if (sum(!is.na(proteins_for_corr)) >= number_of_samples * 0.1 & sum(!is.na(genes_for_corr)) >= number_of_samples * 0.1) {
			cor.test(proteins_for_corr, genes_for_corr, method = "spearman", use = "pairwise.complete.obs")
		} else {
			NA
		}
	}, error = function(e) {
		NA
	})
	if (is.na(corr_result)) {
		correlation_table_rnaseq_by_gene$spearman_corr[i] = NA
		correlation_table_rnaseq_by_gene$pval[i] = NA
	} else {
		correlation_table_rnaseq_by_gene$spearman_corr[i] = corr_result$estimate
		correlation_table_rnaseq_by_gene$pval[i] = corr_result$p.value
	}
}
correlation_table_rnaseq_by_gene$padj = p.adjust(correlation_table_rnaseq_by_gene$pval, method = "BH")


correlation_table_rnaseq_by_gene$plot_color = NA
correlation_table_rnaseq_by_gene$plot_color[correlation_table_rnaseq_by_gene$spearman_corr > 0] = "red"
correlation_table_rnaseq_by_gene$plot_color[correlation_table_rnaseq_by_gene$spearman_corr <= 0] = "blue"

hist(correlation_table_rnaseq_by_gene$spearman_corr, breaks = 50, col = c(rep("blue", 12), rep("red", 37)), main = "", xlab = "", xlim = c(-1, 1), yaxt = "n", ylab = "", ylim = c(0, 400), cex.axis = 3, cex.lab = 3)
axis(side = 2, at = seq(0, 400, 100), cex.axis = 3)
meanspear = mean(na.omit(correlation_table_rnaseq_by_gene$spearman_corr))
abline(v = meanspear, lwd = 2, lty = "dashed")
text(x = meanspear - 0.45, y = 400, paste("Mean = ", round(meanspear, 2)), pos = 1, cex = 3)



correlation_table_rnaseq_by_gene = correlation_table_rnaseq_by_gene[order(correlation_table_rnaseq_by_gene$spearman_corr), ]

png("output/main/pathway_correlation_waterfall_plot_2018-12-06.png", width = 1600, height = 1600)
plot(0, type = "n", xlim = c(1, nrow(correlation_table_rnaseq_by_gene)), ylim = c(-12, 1), xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
axis(side = 2, at = c(-1, 0, 1), cex.axis = 2)
segments(x0 = 1:nrow(correlation_table_rnaseq_by_gene), y0 = 0, x1 = 1:nrow(correlation_table_rnaseq_by_gene), y1 = correlation_table_rnaseq_by_gene$spearman_corr, col = correlation_table_rnaseq_by_gene$plot_color)
#mtext("SpearmanCorrelation", side = 2, adj = 1, line = 2)

text_start_y = -1.5
y0_start = text_start_y - 0.25
y1_start = text_start_y - 0.75


text(x = 0, y = text_start_y, "neutrophil degranulation (GO:0043312)", offset = 0, pos = 4, cex = 3)
neutrophil_genes = pathway_list$`neutrophil degranulation (GO:0043312)`
neutrophil_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% neutrophil_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% neutrophil_df_for_corr_plot$Accession), y0 = y0_start, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% neutrophil_df_for_corr_plot$Accession), y1 = y1_start, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% neutrophil_df_for_corr_plot$Accession)])

text(x = 0, y = text_start_y - 1, "extracellular matrix organization (GO:0030198)", offset = 0, pos = 4, cex = 3)
ecm_genes = pathway_list$`extracellular matrix organization (GO:0030198)`
ecm_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% ecm_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% ecm_df_for_corr_plot$Accession), y0 = y0_start - 1, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% ecm_df_for_corr_plot$Accession), y1 = y1_start - 1, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% ecm_df_for_corr_plot$Accession)])

text(x = 0, y = text_start_y - 2, "negative regulation of apoptotic process (GO:0043066)", offset = 0, pos = 4, cex = 3)
apoptosis_genes = pathway_list$`negative regulation of apoptotic process (GO:0043066)`
apoptosis_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% apoptosis_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% apoptosis_df_for_corr_plot$Accession), y0 = y0_start - 2, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% apoptosis_df_for_corr_plot$Accession), y1 = y1_start - 2, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% apoptosis_df_for_corr_plot$Accession)])

text(x = 0, y = text_start_y - 3, "cellular response to oxidative stress (GO:0034599)", offset = 0, pos = 4, cex = 3)
oxstress_genes = pathway_list$`cellular response to oxidative stress (GO:0034599)`
oxstress_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% oxstress_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% oxstress_df_for_corr_plot$Accession), y0 = y0_start - 3, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% oxstress_df_for_corr_plot$Accession), y1 = y1_start - 3, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% oxstress_df_for_corr_plot$Accession)])

text(x = 0, y = text_start_y - 4, "canonical glycolysis (GO:0061621)", offset = 0, pos = 4, cex = 3)
glycolysis_genes = pathway_list$`canonical glycolysis (GO:0061621)`
glycolysis_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% glycolysis_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% glycolysis_df_for_corr_plot$Accession), y0 = y0_start - 4, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% glycolysis_df_for_corr_plot$Accession), y1 = y1_start - 4, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% glycolysis_df_for_corr_plot$Accession)])


#well correlated
# neutrophil degranulation (GO:0043312) 
# extracellular matrix organization (GO:0030198)
# negative regulation of apoptotic process (GO:0043066)
# cellular response to oxidative stress (GO:0034599)
# canonical glycolysis (GO:0061621)

abline(h = text_start_y - 4.8, lty = 2)
text(x = 0, y = text_start_y - 5, "SRP targeting to membrane (GO:0006614)", offset = 0, pos = 4, cex = 3)
srp_genes = pathway_list$`SRP-dependent cotranslational protein targeting to membrane (GO:0006614)`
srp_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% srp_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% srp_df_for_corr_plot$Accession), y0 = y0_start - 5, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% srp_df_for_corr_plot$Accession), y1 = y1_start - 5, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% srp_df_for_corr_plot$Accession)])

text(x = 0, y = text_start_y - 6, "nonsense-mediated decay (GO:0000184)", offset = 0, pos = 4, cex = 3)
nmd_genes = pathway_list$`nuclear-transcribed mRNA catabolic process, nonsense-mediated decay (GO:0000184)`
nmd_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% nmd_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% nmd_df_for_corr_plot$Accession), y0 = y0_start - 6, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% nmd_df_for_corr_plot$Accession), y1 = y1_start - 6, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% nmd_df_for_corr_plot$Accession)])

text(x = 0, y = text_start_y - 7, "rRNA processing (GO:0006364)", offset = 0, pos = 4, cex = 3)
rrna_genes = pathway_list$`rRNA processing (GO:0006364)`
rrna_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% rrna_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% rrna_df_for_corr_plot$Accession), y0 = y0_start - 7, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% rrna_df_for_corr_plot$Accession), y1 = y1_start - 7, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% rrna_df_for_corr_plot$Accession)])

text(x = 0, y = text_start_y - 8, "translational initiation (GO:0006413)", offset = 0, pos = 4, cex = 3)
transinit_genes = pathway_list$`translational initiation (GO:0006413)`
transinit_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% transinit_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% transinit_df_for_corr_plot$Accession), y0 = y0_start - 8, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% transinit_df_for_corr_plot$Accession), y1 = y1_start - 8, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% transinit_df_for_corr_plot$Accession)])

text(x = 0, y = text_start_y - 9, "translation (GO:0006412)", offset = 0, pos = 4, cex = 3)
transl_genes = pathway_list$`translation (GO:0006412)`
transl_df_for_corr_plot = subset(correlation_table_rnaseq_by_gene, Accession %in% transl_genes)
segments(x0 = which(correlation_table_rnaseq_by_gene$Accession %in% transl_df_for_corr_plot$Accession), y0 = y0_start - 9, x1 = which(correlation_table_rnaseq_by_gene$Accession %in% transl_df_for_corr_plot$Accession), y1 = y1_start - 9, col = correlation_table_rnaseq_by_gene$plot_color[which(correlation_table_rnaseq_by_gene$Accession %in% transl_df_for_corr_plot$Accession)])
dev.off()

#poorly correlated
#SRP-dependent cotranslational protein targeting to membrane (GO:0006614)
#nuclear-transcribed mRNA catabolic process, nonsense-mediated decay (GO:0000184)
#rRNA processing (GO:0006364)
#translational initiation (GO:0006413)
#MAPK cascade (GO:0000165)


#pairwise comparison of sample groups ----


differential_expression = function(expression_in, condition1_samples_in, condition2_samples_in, accession_column = "Symbol", plot_title = "", plot_location = "output/main/") {
    condition1_v_condition2_results_table = data.frame(Accession = expression_in[, accession_column], stringsAsFactors = FALSE)
    condition1_v_condition2_results_table$condition1_log2_mean = apply(expression_in, 1, function(x) {
        mean(na.omit(as.numeric(x[condition1_samples_in])))
    })
    condition1_v_condition2_results_table$condition2_log2_mean = apply(expression_in, 1, function(x){
        mean(na.omit(as.numeric(x[condition2_samples_in])))
    })
    condition1_v_condition2_results_table$condition2_condition1_log2fc = condition1_v_condition2_results_table$condition2_log2_mean - condition1_v_condition2_results_table$condition1_log2_mean
    condition1_v_condition2_results_table$wilcoxon_p = apply(expression_in, 1, function(x) {
        tryCatch({wilcox.test(as.numeric(na.omit(x[condition1_samples_in])), as.numeric(na.omit(x[condition2_samples_in])))$p.value
        },error = function(e) {
            NA
        })
    })
    condition1_v_condition2_results_table$wilcoxon_padj = p.adjust(condition1_v_condition2_results_table$wilcoxon_p, method = "BH")
    condition1_v_condition2_results_table$volcano_colors = NA
    condition1_v_condition2_results_table$volcano_text = NA
    return(condition1_v_condition2_results_table)
}

#redox vs all with fold change
metabolism_proteincc_samples = subset(tmt_meta, consensusClassK5_protein == 2 | consensusClassK5_protein == 4)$TissueID
notmetabolism_proteincc_samples = subset(tmt_meta, consensusClassK5_protein != 2 & consensusClassK5_protein != 4)$TissueID
metabolism_notmetabolism_de = differential_expression(expression_in = protein_expression_10pct, condition1_samples_in = notmetabolism_proteincc_samples, condition2_samples_in = metabolism_proteincc_samples, plot_title = "metabolism vs notmetabolism", plot_location = "output/main/", accession_column = "Symbol")
paste(subset(metabolism_notmetabolism_de, abs(metabolism_notmetabolism_de$condition2_condition1_log2fc) >= log2(1.5) & metabolism_notmetabolism_de$wilcoxon_padj <= 0.05)$Accession, collapse = ",")
metabolism_de_genes = subset(metabolism_notmetabolism_de, abs(metabolism_notmetabolism_de$condition2_condition1_log2fc) >= log2(1.5) & metabolism_notmetabolism_de$wilcoxon_padj <= 0.05)$Accession
metabolismhigher_de_genes = subset(metabolism_notmetabolism_de, metabolism_notmetabolism_de$condition2_condition1_log2fc >= log2(1.5) & metabolism_notmetabolism_de$wilcoxon_padj <= 0.05)$Accession
metabolism_notmetabolism_de_enrichr = enrichr(subset(metabolism_notmetabolism_de, abs(metabolism_notmetabolism_de$condition2_condition1_log2fc) >= log2(1.5) &  metabolism_notmetabolism_de$wilcoxon_padj <= 0.05)$Accession, "GO_Biological_Process_2017")
png(filename = "output/main/pathway_enrichment/metabolism_vs_all_protein_fc.png", width = 800, height = 400)
par(mar =  c(5.1, 38, 2, 2.1))
plot_enrichr_results(enrichr_in = metabolism_notmetabolism_de_enrichr$GO_Biological_Process_2017, number_of_pathways_to_plot = 5, xlim_in = c(0, 40), ylim_in = c(0, 40), plot_title = "redox vs not redox")
dev.off()


metabolismhigher_notmetabolism_de_enrichr = enrichr(subset(metabolism_notmetabolism_de, metabolism_notmetabolism_de$condition2_condition1_log2fc >= log2(1.5) & metabolism_notmetabolism_de$wilcoxon_padj <= 0.05)$Accession, "GO_Biological_Process_2017")
paste(subset(metabolism_notmetabolism_de, metabolism_notmetabolism_de$condition2_condition1_log2fc >= log2(1.5) & metabolism_notmetabolism_de$wilcoxon_padj <= 0.05)$Accession, collapse = ",")
metabolismhigher_notmetabolism_de_proteins = subset(metabolism_notmetabolism_de, metabolism_notmetabolism_de$condition2_condition1_log2fc >= log2(1.5) & metabolism_notmetabolism_de$wilcoxon_padj <= 0.05)$Accession
metabolism_notmetabolismhigher_de_enrichr = enrichr(subset(metabolism_notmetabolism_de, metabolism_notmetabolism_de$condition2_condition1_log2fc <= log2(0.6667) & metabolism_notmetabolism_de$wilcoxon_padj <= 0.05)$Accession, "GO_Biological_Process_2017")


png(filename = "output/main/pathway_enrichment/metabolismhigher_vs_notmetabolism_protein_fc.png", width = 900, height = 200)
par(mar =  c(5.1, 38, 2, 2.1))
plot_enrichr_results(enrichr_in = metabolismhigher_notmetabolism_de_enrichr$GO_Biological_Process_2017, number_of_pathways_to_plot = 5, xlim_in = c(0, 300), ylim_in = c(0, 10), plot_title = "")
dev.off()


#inflamed vs all with fold change
immune_proteincc_samples = subset(tmt_meta, consensusClassK5_protein == 1 | consensusClassK5_protein == 3)$TissueID
notimmune_proteincc_samples = subset(tmt_meta, consensusClassK5_protein != 1 & consensusClassK5_protein != 3)$TissueID
immune_notimmune_de = differential_expression(expression_in = protein_expression_10pct, condition1_samples_in = notimmune_proteincc_samples, condition2_samples_in = immune_proteincc_samples, plot_title = "immune vs notimmune", plot_location = "output/main/", accession_column = "Symbol")
paste(subset(immune_notimmune_de, abs(immune_notimmune_de$condition2_condition1_log2fc) >= log2(1.5) & immune_notimmune_de$wilcoxon_padj <= 0.05)$Accession, collapse = ",")
immune_de_genes = subset(immune_notimmune_de, abs(immune_notimmune_de$condition2_condition1_log2fc) >= log2(1.5) & immune_notimmune_de$wilcoxon_padj <= 0.05)$Accession
immunehigher_de_genes = subset(immune_notimmune_de, immune_notimmune_de$condition2_condition1_log2fc >= log2(1.5) & immune_notimmune_de$wilcoxon_padj <= 0.05)$Accession
immune_notimmune_de_enrichr = enrichr(subset(immune_notimmune_de, abs(immune_notimmune_de$condition2_condition1_log2fc) >= log2(1.5) &  immune_notimmune_de$wilcoxon_padj <= 0.05)$Accession, "GO_Biological_Process_2017")
png(filename = "output/main/pathway_enrichment/immune_vs_all_protein_fc.png", width = 800, height = 400)
par(mar =  c(5.1, 38, 2, 2.1))
plot_enrichr_results(enrichr_in = immune_notimmune_de_enrichr$GO_Biological_Process_2017, number_of_pathways_to_plot = 5, xlim_in = c(0, 40), ylim_in = c(0, 10), plot_title = "inflamed vs not inflamed")
dev.off()

immunehigher_notimmune_de_enrichr = enrichr(subset(immune_notimmune_de, immune_notimmune_de$condition2_condition1_log2fc >= log2(1.5) & immune_notimmune_de$wilcoxon_padj <= 0.05)$Accession, "GO_Biological_Process_2017")
immunehigher_notimmune_de_proteins = subset(immune_notimmune_de, immune_notimmune_de$condition2_condition1_log2fc >= log2(1.5) & immune_notimmune_de$wilcoxon_padj <= 0.05)$Accession
paste(subset(immune_notimmune_de, immune_notimmune_de$condition2_condition1_log2fc >= log2(1.5) & immune_notimmune_de$wilcoxon_padj <= 0.05)$Accession, collapse = ",")
immune_notimmunehigher_de_enrichr = enrichr(subset(immune_notimmune_de, immune_notimmune_de$condition2_condition1_log2fc <= log2(0.6667) & immune_notimmune_de$wilcoxon_padj <= 0.05)$Accession, "GO_Biological_Process_2017")

png(filename = "output/main/pathway_enrichment/immunehigher_vs_notimmune_protein_fc.png", width = 900, height = 200)
par(mar =  c(5.1, 38, 2, 2.1))
plot_enrichr_results(enrichr_in = immunehigher_notimmune_de_enrichr$GO_Biological_Process_2017, number_of_pathways_to_plot = 5, xlim_in = c(0, 300), ylim_in = c(0, 10), plot_title = "")
dev.off()


#Mixed vs all with fold change
Mixed_proteincc_samples = subset(tmt_meta, consensusClassK5_protein == 5)$TissueID
notMixed_proteincc_samples = subset(tmt_meta, consensusClassK5_protein != 5)$TissueID
Mixed_notMixed_de = differential_expression(expression_in = protein_expression_10pct, condition1_samples_in = notMixed_proteincc_samples, condition2_samples_in = Mixed_proteincc_samples, plot_title = "Mixed vs notMixed", plot_location = "output/main/", accession_column = "Symbol")
paste(subset(Mixed_notMixed_de, abs(Mixed_notMixed_de$condition2_condition1_log2fc) >= log2(1.5) & Mixed_notMixed_de$wilcoxon_padj <= 0.05)$Accession, collapse = ",")
Mixed_de_genes = subset(Mixed_notMixed_de, abs(Mixed_notMixed_de$condition2_condition1_log2fc) >= log2(1.5) & Mixed_notMixed_de$wilcoxon_padj <= 0.05)$Accession
Mixedhigher_de_genes = subset(Mixed_notMixed_de, Mixed_notMixed_de$condition2_condition1_log2fc >= log2(1.5) & Mixed_notMixed_de$wilcoxon_padj <= 0.05)$Accession
Mixed_notMixed_de_enrichr = enrichr(subset(Mixed_notMixed_de, abs(Mixed_notMixed_de$condition2_condition1_log2fc) >= log2(1.5) &  Mixed_notMixed_de$wilcoxon_padj <= 0.05)$Accession, "GO_Biological_Process_2017")
png(filename = "output/main/pathway_enrichment/Mixed_vs_all_protein_fc.png", width = 800, height = 200)
par(mar =  c(5.1, 38, 2, 2.1))
plot_enrichr_results(enrichr_in = Mixed_notMixed_de_enrichr$GO_Biological_Process_2017, number_of_pathways_to_plot = 5, xlim_in = c(0, 40), ylim_in = c(0, 40), plot_title = "Mixed vs all")
dev.off()

Mixedhigher_notMixed_de_enrichr = enrichr(subset(Mixed_notMixed_de, Mixed_notMixed_de$condition2_condition1_log2fc >= log2(1.5) & Mixed_notMixed_de$wilcoxon_padj <= 0.05)$Accession, "GO_Biological_Process_2017")
Mixedhigher_notMixed_de_proteins = subset(Mixed_notMixed_de, Mixed_notMixed_de$condition2_condition1_log2fc >= log2(1.5) & Mixed_notMixed_de$wilcoxon_padj <= 0.05)$Accession
Mixed_notMixedhigher_de_enrichr = enrichr(subset(Mixed_notMixed_de, Mixed_notMixed_de$condition2_condition1_log2fc <= log2(0.6667) & Mixed_notMixed_de$wilcoxon_padj <= 0.05)$Accession, "GO_Biological_Process_2017")

png(filename = "output/main/pathway_enrichment/Mixedhigher_vs_notMixed_protein_fc.png", width = 900, height = 200)
par(mar =  c(5.1, 38, 2, 2.1))
plot_enrichr_results(enrichr_in = Mixedhigher_notMixed_de_enrichr$GO_Biological_Process_2017, number_of_pathways_to_plot = 5, xlim_in = c(0, 300), ylim_in = c(0, 10), plot_title = "")
dev.off()

#ESTIMATE scores ----
#for headers with characters that R doesn't like, need to edit the score names by hand and replace some "." with "-"; checknames argument is not being run in the above ESTIMATE commands so the RNA sample names in the header might not match the ESTIMATE output

estimate_scores = read.delim("output/main/estimate_rnaseq_output_scores_2018-05-10.gct", header = TRUE, row.names = 1, stringsAsFactors = FALSE, skip = 2, check.names = FALSE)
estimate_scores$Description = NULL

Inflamed_stromal = estimate_scores["StromalScore", tmt_meta$GenomicsSample[tmt_meta$consensusClassK5_protein == 1 | tmt_meta$consensusClassK5_protein == 3]]
Inflamed_immune = estimate_scores["ImmuneScore", tmt_meta$GenomicsSample[tmt_meta$consensusClassK5_protein == 1 | tmt_meta$consensusClassK5_protein == 3]]
Inflamed_estimate = estimate_scores["ESTIMATEScore", tmt_meta$GenomicsSample[tmt_meta$consensusClassK5_protein == 1 | tmt_meta$consensusClassK5_protein == 3]]
Redox_stromal = estimate_scores["StromalScore", tmt_meta$GenomicsSample[tmt_meta$consensusClassK5_protein == 4  | tmt_meta$consensusClassK5_protein == 4]]
Redox_immune = estimate_scores["ImmuneScore", tmt_meta$GenomicsSample[tmt_meta$consensusClassK5_protein == 4  | tmt_meta$consensusClassK5_protein == 4]]
Redox_estimate = estimate_scores["ESTIMATEScore", tmt_meta$GenomicsSample[tmt_meta$consensusClassK5_protein == 4  | tmt_meta$consensusClassK5_protein == 4]]
Mixed_stromal = estimate_scores["StromalScore", tmt_meta$GenomicsSample[tmt_meta$consensusClassK5_protein == 5]]
Mixed_immune = estimate_scores["ImmuneScore", tmt_meta$GenomicsSample[tmt_meta$consensusClassK5_protein == 5]]
Mixed_estimate = estimate_scores["ESTIMATEScore", tmt_meta$GenomicsSample[tmt_meta$consensusClassK5_protein == 5]]


png(filename = "output/main/ESTIMATE_2018-06-18.png", width = 900, height = 500)
par(mfrow = c(1,3))
par(mar = c(8.1, 5.1, 4.1, 2.1))

boxplot(list("Inflamed" = as.numeric(Inflamed_immune), "Redox" = as.numeric(Redox_immune), "Mixed" = as.numeric(Mixed_immune)), ylab = "Score", main = "", cex.lab = 2, cex.axis = 2, col = c(darkcolors[1], darkcolors[2], darkcolors[3]))

boxplot(list("Inflamed" = as.numeric(Inflamed_stromal), "Redox" = as.numeric(Redox_stromal), "Mixed" = as.numeric(Mixed_stromal)), ylab = "Score", main = "", cex.lab = 2, cex.axis = 2, col = c(darkcolors[1], darkcolors[2], darkcolors[3]))

boxplot(list("Inflamed" = as.numeric(Inflamed_estimate), "Redox" = as.numeric(Redox_estimate), "Mixed" = as.numeric(Mixed_estimate)), ylab = "Score", main = "", cex.lab = 2, cex.axis = 2, col = c(darkcolors[1], darkcolors[2], darkcolors[3]))
dev.off()


#CNV data and prep RNA and protein for integration with CNV ----
cnv_data = read.delim("data/main/Gene_CNV_400kbp.50counts_noXYgene_2018-02-12.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
cnv_meta = read.delim("data/main/original_data/cytoscan_run_batch_meta_09-25-17.txt", header = TRUE, stringsAsFactors = FALSE)
cnv_meta$duplicated = duplicated(cnv_meta$TissueID)
cnv_meta_duplicated = subset(cnv_meta, duplicated == TRUE)
cnv_meta_duplicated = subset(cnv_meta_duplicated, !(TissueID == "06S15101193" & Attempt == 2)) #duplicated twice in attempt 2 and attempt 3; remove
cnv_meta = subset(cnv_meta, !(TissueID %in% cnv_meta_duplicated$TissueID) & Attempt == 1 & DNA.Aliquot != "Positive_control") #get non-duplicated and non-control entries
cnv_meta = rbind(cnv_meta, cnv_meta_duplicated)



names(cnv_data) = str_replace(names(cnv_data), "^X0", "0")
tmt_meta_for_cnv = subset(tmt_meta, TissueID %in% names(cnv_data))
tmt_meta_for_cnv = merge(tmt_meta_for_cnv, cnv_meta, by = "TissueID")
cnv_data = cnv_data[, tmt_meta_for_cnv$TissueID]
rnaseq_expression_for_cnv = rnaseq_expression
row.names(rnaseq_expression_for_cnv) = rnaseq_expression_for_cnv$Symbol
rnaseq_expression_for_cnv$Symbol = NULL
rnaseq_expression_for_cnv = rnaseq_expression_for_cnv[, tmt_meta_for_cnv$GenomicsSample]
names(cnv_data) = tmt_meta_for_cnv$TissueID
names(rnaseq_expression_for_cnv) = tmt_meta_for_cnv$TissueID
protein_expression_for_cnv = protein_expression_10pct[, c("Symbol", tmt_meta_for_cnv$TissueID)]
protein_expression_for_cnv = protein_expression_for_cnv[!duplicated(protein_expression_for_cnv$Symbol), ]
row.names(protein_expression_for_cnv) = protein_expression_for_cnv$Symbol

cnv_protein_merged = merge(cnv_data, protein_expression_for_cnv, by.x = "row.names", by.y = "Symbol")

#CNV correlation with protein ----
cnv_protein_merged[cnv_protein_merged == 0] = NA
correlation_table_cnv_protein = data.frame(Accession = cnv_protein_merged$Row.names, spearman_corr = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(cnv_protein_merged)) {
	corr_result = tryCatch({
		cnv_for_corr = as.numeric(cnv_protein_merged[i, paste(tmt_meta_for_cnv$TissueID, "x", sep = ".")])
		proteins_for_corr = as.numeric(cnv_protein_merged[i, paste(tmt_meta_for_cnv$TissueID, "y", sep = ".")])
		if (sum(!is.na(proteins_for_corr)) >= number_of_samples * 0.1 & sum(!is.na(cnv_for_corr)) >= number_of_samples * 0.1) {
			cor.test(proteins_for_corr, cnv_for_corr, method = "spearman", use = "pairwise.complete.obs")
		} else {
			NA
		}
	}, error = function(e) {
		NA
	})
	if (is.na(corr_result)) {
		correlation_table_cnv_protein$spearman_corr[i] = NA
		correlation_table_cnv_protein$pval[i] = NA
	} else {
		correlation_table_cnv_protein$spearman_corr[i] = corr_result$estimate
		correlation_table_cnv_protein$pval[i] = corr_result$p.value
	}
}
correlation_table_cnv_protein$padj = p.adjust(correlation_table_cnv_protein$pval, method = "BH")

correlation_table_cnv_protein$hist_color = NA
correlation_table_cnv_protein$hist_color[correlation_table_cnv_protein$spearman_corr > 0] = "red"
correlation_table_cnv_protein$hist_color[correlation_table_cnv_protein$spearman_corr <= 0] = "blue"

hist(correlation_table_cnv_protein$spearman_corr, breaks = 50, col = c(rep("blue", 20), rep("red", 37)), main = "", xlab = "Spearman's Correlation", xlim = c(-1, 1), yaxt = "n", ylab = "Count", ylim = c(0, 300), cex.axis = 1.5, cex.lab = 1.5)
axis(side = 2, at = seq(0, 500, 100), cex.axis = 1.5)
meanspear = mean(na.omit(correlation_table_cnv_protein$spearman_corr))
abline(v = meanspear, lwd = 2, lty = "dashed")
text(x = meanspear + .4, y = 315, paste("Mean = ", round(meanspear, 2)), pos = 1, cex = 1.5)


#poster

png("output/main/rnaseq_protein_cnv_protein_correlation_histograms_2018-06-21.png", width = 1200, height = 800)
par(mfrow = c(1, 2))

hist(correlation_table_cnv_protein$spearman_corr, breaks = 50, col = c(rep("blue", 28), rep("red", 37)), main = "", xlab = "", xlim = c(-1, 1), yaxt = "n", ylab = "", ylim = c(0, 400), cex.axis = 2, cex.lab = 2)
axis(side = 2, at = seq(0, 400, 100), cex.axis = 2)
meanspear = mean(na.omit(correlation_table_cnv_protein$spearman_corr))

hist(correlation_table_rnaseq_by_gene$spearman_corr, breaks = 50, col = c(rep("blue", 21), rep("red", 37)), main = "", xlab = "", xlim = c(-1, 1), yaxt = "n", ylab = "", ylim = c(0, 400), cex.axis = 2, cex.lab = 2)
axis(side = 2, at = seq(0, 400, 100), cex.axis = 2)

dev.off()


#combine rna-protein and cn-protein tables for paper
correlation_table_rna_protein_cn_protein = merge(correlation_table_rnaseq_by_gene, correlation_table_cnv_protein, by = "Accession", all = TRUE)


correlation_table_cnv_protein_positive = subset(correlation_table_cnv_protein, spearman_corr > 0.5 & padj < 0.25) 
correlation_table_rnaseq_by_gene_positive = subset(correlation_table_rnaseq_by_gene, spearman_corr > 0.5 & padj < 0.25) 
correlation_table_cnv_protein_positive_matched = subset(correlation_table_cnv_protein_positive, Accession %in% correlation_table_rnaseq_by_gene_positive$Accession)
dim(correlation_table_cnv_protein_positive_matched)


metabolism_correlations = correlation_table_cnv_protein_positive_matched[correlation_table_cnv_protein_positive_matched$hgnc_symbol %in% metabolismhigher_de_genes, ]

metabolism_correlations[order(metabolism_correlations$chromosome_name), ]

#TCGA mut status comparisons ----

tcga_meta = read.delim("data/main/tcga_classification_with_mut_2018-03-06.txt", header = TRUE, stringsAsFactors = FALSE)
tcga_meta[is.na(tcga_meta)] = "wt"
tcga_meta[tcga_meta == ""] = "wt"

table(subset(tcga_meta, tcga_meta$Classification == "secretory")$TP53 == "wt")
table(subset(tcga_meta, tcga_meta$Classification != "secretory")$TP53 == "wt")

table(subset(tmt_meta, subtype == "Inflamed")$TP53 == "wt")
table(subset(tmt_meta, subtype != "Inflamed")$TP53 == "wt")

table(subset(tcga_meta, Classification == "classical")$NFE2L2 != "wt")

table(subset(tmt_meta, subtype == "Redox")$NFE2L2 != "wt")
table(subset(tmt_meta, subtype == "Redox")$KEAP1 != "wt")

nrow(subset(tmt_meta, NFE2L2 != "wt" | KEAP1 != "wt"))
table(subset(tmt_meta, NFE2L2 != "wt" | KEAP1 != "wt")$subtype)

#boxplots protein ----

inflamed_proteincc_samples = subset(tmt_meta, consensusClassK5_protein == 1 | consensusClassK5_protein == 3)$TissueID
redox_proteincc_samples = subset(tmt_meta, consensusClassK5_protein == 2 | consensusClassK5_protein == 4)$TissueID
mixed_proteincc_samples = subset(tmt_meta, consensusClassK5_protein == 5)$TissueID

gene_symbol_boxplot = "CD40"
boxplot(list(Inflamed = as.numeric(protein_expression[protein_expression$Symbol == gene_symbol_boxplot, inflamed_proteincc_samples]), Redox = as.numeric(protein_expression[protein_expression$Symbol == gene_symbol_boxplot, redox_proteincc_samples]), Mixed = as.numeric(protein_expression[protein_expression$Symbol == gene_symbol_boxplot, mixed_proteincc_samples])), main = paste(gene_symbol_boxplot, "", sep = " "), ylab = "Log2 Expression")


#CNV regions for heatmaps  ----
cnv_long = read.delim("data/main/Standard_400kbp.50counts_updated_processed_20171220.txt", header = TRUE, stringsAsFactors = FALSE)
cnv_long_3q26.1 = cnv_long[str_detect(cnv_long$chromosome_band, "^3\\_q26\\.1"), ]
tmt_meta$cnv_3q261 = tmt_meta$TissueID %in% cnv_long_3q26.1$TissueID
cnv_long_3q26 = cnv_long[str_detect(cnv_long$chromosome_band, "^3\\_q26"), ]
tmt_meta$cnv_3q26 = tmt_meta$TissueID %in% cnv_long_3q26$TissueID
cnv_long_3q2x = cnv_long[str_detect(cnv_long$chromosome_band, "^3\\_q2[5-7]"), ]
tmt_meta$cnv_3q2x = tmt_meta$TissueID %in% cnv_long_3q2x$TissueID
cnv_long_1q31 = cnv_long[str_detect(cnv_long$chromosome_band, "^1\\_q31"), ]
tmt_meta$cnv_1q31 = tmt_meta$TissueID %in% cnv_long_1q31$TissueID
cnv_long_5p1x = cnv_long[str_detect(cnv_long$chromosome_band, "^5\\_p1[3-4]"), ]
tmt_meta$cnv_5p1x = tmt_meta$TissueID %in% cnv_long_5p1x$TissueID


#cnv count plot ----
gt_x_samples = 0
chromosome_coords = read.delim("data/main/chromosome_start_end_from_hg19_2018-04-05.txt", header = TRUE, stringsAsFactors = FALSE)
cnv_long = read.delim("data/main/Standard_400kbp.50counts_updated_processed_20171220.txt", header = TRUE, stringsAsFactors = FALSE)

#subset by subtype needs to be commented out to get summary gains/losses for manuscript

#pick one
#cnv_long = subset(cnv_long, TissueID %in% subset(tmt_meta, subtype == "Redox")$TissueID)
#cnv_long = subset(cnv_long, TissueID %in% subset(tmt_meta, subtype == "Inflamed")$TissueID)
#cnv_long = subset(cnv_long, TissueID %in% subset(tmt_meta, subtype == "Mixed")$TissueID)


cnv_long = subset(cnv_long, TissueID %in% tmt_meta$TissueID)
cnv_long$chromosome_band = paste("0", cnv_long$chromosome_band, sep = "")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "010", "10")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "011", "11")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "012", "12")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "013", "13")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "014", "14")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "015", "15")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "016", "16")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "017", "17")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "018", "18")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "019", "19")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "020", "20")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "021", "21")
cnv_long$chromosome_band = str_replace(cnv_long$chromosome_band, "022", "22")


cnv_coords = data.frame(rbindlist(lapply(str_split(cnv_long$Full.Location, ":|-"), function(x) {data.frame(t(x), stringsAsFactors = FALSE)})), stringsAsFactors = FALSE) #split the full location into a list

names(cnv_coords) = c("chromosome", "startpos", "endpos")
cnv_coords$chromosome = NULL #delete chromosome information since it's repeated
cnv_coords$startpos = as.numeric(cnv_coords$startpos)
cnv_coords$endpos = as.numeric(cnv_coords$endpos)
cnv_long = data.frame(cbind(cnv_long, cnv_coords), stringsAsFactors = FALSE)

#losses plot setup
cnv_long_cnstate_loss = subset(cnv_long, Type == "Loss") 
cnv_long_cnstate_loss = cnv_long_cnstate_loss[order(cnv_long_cnstate_loss$chromosome_band), ]

cnv_long_cnstate_loss_unique = unique(cnv_long_cnstate_loss[, c("TissueID", "chromosome_band")]) #TissueID's can have repeated chromosome + band combinations which overinflates the counts
cnv_long_cnstate_loss_bandfreq = data.frame(table(cnv_long_cnstate_loss_unique$chromosome_band), stringsAsFactors = FALSE)
names(cnv_long_cnstate_loss_bandfreq) = c("chromosome_band", "Freq")
cnv_long_cnstate_loss_bandfreq$chromosome_band = as.character(cnv_long_cnstate_loss_bandfreq$chromosome_band)
cnv_long_cnstate_loss_bandfreq = subset(cnv_long_cnstate_loss_bandfreq, Freq > gt_x_samples) #only include counts greater than x samples defined above


cnv_long_cnstate_loss_bandfreq$Chromosome = as.numeric(str_match(cnv_long_cnstate_loss_bandfreq$chromosome_band, "^[0-9]*"))
cnv_long_cnstate_loss_bandfreq$Genes = ""
cnv_long_cnstate_loss_bandfreq$startpos = 0
cnv_long_cnstate_loss_bandfreq$endpos = 0
for(i in 1:nrow(cnv_long_cnstate_loss_bandfreq)) {
	cnv_long_cnstate_loss_bandfreq$Genes[i] = paste(subset(cnv_long_cnstate_loss, chromosome_band == cnv_long_cnstate_loss_bandfreq$chromosome_band[i])$Genes, collapse = ", ") #populate genes in cnv bandfreq dataframe to use for annotating the plot
	cnv_long_cnstate_loss_bandfreq$startpos[i] = mean(subset(cnv_long_cnstate_loss, chromosome_band == cnv_long_cnstate_loss_bandfreq$chromosome_band[i])$startpos)
	cnv_long_cnstate_loss_bandfreq$endpos[i] = mean(subset(cnv_long_cnstate_loss, chromosome_band == cnv_long_cnstate_loss_bandfreq$chromosome_band[i])$endpos)
}


plot_losses_df = list()
for (i in 1:22) {
	cnv_long_cnstate_loss_bandfreq_sub = subset(cnv_long_cnstate_loss_bandfreq, Chromosome == i)
	if (nrow(cnv_long_cnstate_loss_bandfreq_sub) > 0) {
		chrom_lowerbound = min(cnv_long_cnstate_loss_bandfreq_sub$startpos)
		chrom_upperbound = chromosome_coords$end[i]
		cnv_long_cnstate_loss_bandfreq_sub$y_values = i + cnv_long_cnstate_loss_bandfreq_sub$startpos / chrom_upperbound
		plot_losses_df[[i]] = data.frame(x_values = cnv_long_cnstate_loss_bandfreq_sub$Freq, y_values = cnv_long_cnstate_loss_bandfreq_sub$y_values, chromosome_band = cnv_long_cnstate_loss_bandfreq_sub$chromosome_band, Genes = cnv_long_cnstate_loss_bandfreq_sub$Genes, stringsAsFactors = FALSE)
	}
}
plot_losses_df = as.data.frame(rbindlist(plot_losses_df))
plot_losses_df$x_values = plot_losses_df$x_values / number_of_samples * 100


#gains plot setup
cnv_long_cnstate_gain = subset(cnv_long, Type == "Gain") 
cnv_long_cnstate_gain = cnv_long_cnstate_gain[order(cnv_long_cnstate_gain$chromosome_band), ]

cnv_long_cnstate_gain_unique = unique(cnv_long_cnstate_gain[, c("TissueID", "chromosome_band")]) #TissueID's can have repeated chromosome + band combinations which overinflates the counts
cnv_long_cnstate_gain_bandfreq = data.frame(table(cnv_long_cnstate_gain_unique$chromosome_band), stringsAsFactors = FALSE)
names(cnv_long_cnstate_gain_bandfreq) = c("chromosome_band", "Freq")
cnv_long_cnstate_gain_bandfreq$chromosome_band = as.character(cnv_long_cnstate_gain_bandfreq$chromosome_band)
cnv_long_cnstate_gain_bandfreq = subset(cnv_long_cnstate_gain_bandfreq, Freq > gt_x_samples) #only include counts greater than x samples defined above


cnv_long_cnstate_gain_bandfreq$Chromosome = as.numeric(str_match(cnv_long_cnstate_gain_bandfreq$chromosome_band, "^[0-9]*"))
cnv_long_cnstate_gain_bandfreq$Genes = ""
cnv_long_cnstate_gain_bandfreq$startpos = 0
cnv_long_cnstate_gain_bandfreq$endpos = 0
for(i in 1:nrow(cnv_long_cnstate_gain_bandfreq)) {
	cnv_long_cnstate_gain_bandfreq$Genes[i] = paste(subset(cnv_long_cnstate_gain, chromosome_band == cnv_long_cnstate_gain_bandfreq$chromosome_band[i])$Genes, collapse = ", ") #populate genes in cnv bandfreq dataframe to use for annotating the plot
	cnv_long_cnstate_gain_bandfreq$startpos[i] = mean(subset(cnv_long_cnstate_gain, chromosome_band == cnv_long_cnstate_gain_bandfreq$chromosome_band[i])$startpos)
	cnv_long_cnstate_gain_bandfreq$endpos[i] = mean(subset(cnv_long_cnstate_gain, chromosome_band == cnv_long_cnstate_gain_bandfreq$chromosome_band[i])$endpos)
}


plot_gains_df = list()
for (i in 1:22) {
	cnv_long_cnstate_gain_bandfreq_sub = subset(cnv_long_cnstate_gain_bandfreq, Chromosome == i)
	if (nrow(cnv_long_cnstate_gain_bandfreq_sub) > 0) {
		chrom_upperbound = chromosome_coords$end[i]
		cnv_long_cnstate_gain_bandfreq_sub$y_values = i + cnv_long_cnstate_gain_bandfreq_sub$startpos / chrom_upperbound
		plot_gains_df[[i]] = data.frame(x_values = cnv_long_cnstate_gain_bandfreq_sub$Freq, y_values = cnv_long_cnstate_gain_bandfreq_sub$y_values, chromosome_band = cnv_long_cnstate_gain_bandfreq_sub$chromosome_band, Genes = cnv_long_cnstate_gain_bandfreq_sub$Genes, stringsAsFactors = FALSE)
	}
}
plot_gains_df = as.data.frame(rbindlist(plot_gains_df))
plot_gains_df$x_values = plot_gains_df$x_values / number_of_samples * 100



annotate_regions_in_cnv_plot = TRUE

png(filename = "output/main/cnv_gain_loss_counts_2018-12-11.png", height = 800, width = 1000)
par(mar = c(5.1, 1, 4.1, 1.5))
par(mfrow = c(1,2))
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5, cex.sub = 1.5)
cex_text_parameter = 1.5
plot(0, type = "n", ylim = c(26, 0), bty = "n", xaxt = "n", yaxt = "n", xlab = "Percent Patients with Loss", xlim = c(55, 0))
for (i in 1:nrow(plot_losses_df)) {
	segments(x0 = 0, y0 = plot_losses_df$y_values[i], x1 = plot_losses_df$x_values[i], y1 = plot_losses_df$y_values[i], col = "blue")
}

axis(side = 4, at = c(1:23), lab = c(1:22, ""), las = 2) #y axis with labels
axis(side = 1, at = seq(50, 0, -10), line = -2) #x axis

if (annotate_regions_in_cnv_plot == TRUE) {

	KEAP1_locations = subset(plot_losses_df, str_detect(plot_losses_df$Genes, " KEAP1,|^KEAP1, |, KEAP1$"))
	KEAP1_max_x  = max(KEAP1_locations$x_values)
	KEAP1_min_y = min(KEAP1_locations$y_values)
	KEAP1_max_y = max(KEAP1_locations$y_values)
	
	segments(x0 = KEAP1_max_x + 2.5, y0 = KEAP1_min_y, x1 = KEAP1_max_x + 5, y1 = KEAP1_min_y) #top horizontal line
	segments(x0 = KEAP1_max_x + 5, y0 = KEAP1_min_y, x1 = KEAP1_max_x + 5, y1 = KEAP1_max_y) #vertical line
	segments(x0 = KEAP1_max_x + 5, y0 = KEAP1_max_y, x1 = KEAP1_max_x + 2.5, y1 = KEAP1_max_y) #bottom horizontal line
	segments(x0 = KEAP1_max_x + 5, y0 = (KEAP1_min_y + KEAP1_max_y) / 2, x1 = KEAP1_max_x + 7.5, y1 = (KEAP1_min_y + KEAP1_max_y) / 2) #midpoint horizontal line
	text("KEAP1", x = KEAP1_max_x + 7.5, y = (KEAP1_min_y + KEAP1_max_y) / 2, pos = 2, offset = 0, cex = cex_text_parameter)
	
	
	chrom03_p_locations = subset(plot_losses_df, str_detect(plot_losses_df$chromosome_band, "03\\_p"))
	chrom03_p_max_x  = max(chrom03_p_locations$x_values)
	chrom03_p_min_y = min(chrom03_p_locations$y_values)
	chrom03_p_max_y = max(chrom03_p_locations$y_values)
	
	segments(x0 = chrom03_p_max_x + 2.5, y0 = chrom03_p_min_y, x1 = chrom03_p_max_x + 5, y1 = chrom03_p_min_y) #top horizontal line
	segments(x0 = chrom03_p_max_x + 5, y0 = chrom03_p_min_y, x1 = chrom03_p_max_x + 5, y1 = chrom03_p_max_y) #vertical line
	segments(x0 = chrom03_p_max_x + 5, y0 = chrom03_p_max_y, x1 = chrom03_p_max_x + 2.5, y1 = chrom03_p_max_y) #bottom horizontal line
	segments(x0 = chrom03_p_max_x + 5, y0 = (chrom03_p_min_y + chrom03_p_max_y) / 2, x1 = chrom03_p_max_x + 7.5, y1 = (chrom03_p_min_y + chrom03_p_max_y) / 2) #midpoint horizontal line
	text("3p", x = chrom03_p_max_x + 7.5, y = (chrom03_p_min_y + chrom03_p_max_y) / 2, pos = 2, offset = 0, cex = cex_text_parameter)
	
	chrom05_q_locations = subset(plot_losses_df, str_detect(plot_losses_df$chromosome_band, "05\\_q"))
	chrom05_q_max_x  = max(chrom05_q_locations$x_values)
	chrom05_q_min_y = min(chrom05_q_locations$y_values)
	chrom05_q_max_y = max(chrom05_q_locations$y_values)
	
	segments(x0 = chrom05_q_max_x + 2.5, y0 = chrom05_q_min_y, x1 = chrom05_q_max_x + 5, y1 = chrom05_q_min_y) #top horizontal line
	segments(x0 = chrom05_q_max_x + 5, y0 = chrom05_q_min_y, x1 = chrom05_q_max_x + 5, y1 = chrom05_q_max_y) #vertical line
	segments(x0 = chrom05_q_max_x + 5, y0 = chrom05_q_max_y, x1 = chrom05_q_max_x + 2.5, y1 = chrom05_q_max_y) #bottom horizontal line
	segments(x0 = chrom05_q_max_x + 5, y0 = (chrom05_q_min_y + chrom05_q_max_y) / 2, x1 = chrom05_q_max_x + 7.5, y1 = (chrom05_q_min_y + chrom05_q_max_y) / 2) #midpoint horizontal line
	text("5q", x = chrom05_q_max_x + 7.5, y = (chrom05_q_min_y + chrom05_q_max_y) / 2, pos = 2, offset = 0, cex = cex_text_parameter)
	
	chrom09_p_locations = subset(plot_losses_df, str_detect(plot_losses_df$chromosome_band, "09\\_p"))
	chrom09_p_max_x  = max(chrom09_p_locations$x_values)
	chrom09_p_min_y = min(chrom09_p_locations$y_values)
	chrom09_p_max_y = max(chrom09_p_locations$y_values)
	
	segments(x0 = chrom09_p_max_x + 2.5, y0 = chrom09_p_min_y, x1 = chrom09_p_max_x + 5, y1 = chrom09_p_min_y) #top horizontal line
	segments(x0 = chrom09_p_max_x + 5, y0 = chrom09_p_min_y, x1 = chrom09_p_max_x + 5, y1 = chrom09_p_max_y) #vertical line
	segments(x0 = chrom09_p_max_x + 5, y0 = chrom09_p_max_y, x1 = chrom09_p_max_x + 2.5, y1 = chrom09_p_max_y) #bottom horizontal line
	segments(x0 = chrom09_p_max_x + 5, y0 = (chrom09_p_min_y + chrom09_p_max_y) / 2, x1 = chrom09_p_max_x + 7.5, y1 = (chrom09_p_min_y + chrom09_p_max_y) / 2) #midpoint horizontal line
	text("9p", x = chrom09_p_max_x + 7.5, y = (chrom09_p_min_y + chrom09_p_max_y) / 2, pos = 2, offset = 0, cex = cex_text_parameter)
}


par(mar = c(5.1, 1.5, 4.1, 1))
#plot(plot_gains_df$x_values, plot_gains_df$y_values, ylim = c(26, 0), bty = "n", xaxt = "n", yaxt = "n", xlab = "Percentage Gain", xlim = c(10, 50), type = "l", col = "red")
plot(0, type = "n", ylim = c(26, 0), bty = "n", xaxt = "n", yaxt = "n", xlab = "Percent Patients with Gain", xlim = c(0, 55))
for (i in 1:nrow(plot_gains_df)) {
	segments(x0 = 0, y0 = plot_gains_df$y_values[i], x1 = plot_gains_df$x_values[i], y1 = plot_gains_df$y_values[i], col = "red")	
}


axis(side = 2, at = c(1:23), lab = rep("", 23)) #yaxis without labels
#axis(side = 1, at = seq(10, 70, 10), line = -0.5) #x axis
axis(side = 1, at = seq(0, 50, 10), line = -2) #x axis

if (annotate_regions_in_cnv_plot == TRUE) {
	
	#plot p63 sox2 pik3ca locations
	p63_sox2_pik3ca_locations = subset(plot_gains_df, str_detect(plot_gains_df$Genes, " TP63,| SOX2,| PIK3CA,"))
	p63_sox2_pik3ca_max_x  = max(p63_sox2_pik3ca_locations$x_values) + 10
	p63_sox2_pik3ca_min_y = min(p63_sox2_pik3ca_locations$y_values)
	p63_sox2_pik3ca_max_y = max(p63_sox2_pik3ca_locations$y_values)
	
	segments(x0 = p63_sox2_pik3ca_max_x + 2.5, y0 = p63_sox2_pik3ca_min_y, x1 = p63_sox2_pik3ca_max_x + 5, y1 = p63_sox2_pik3ca_min_y) #top horizontal line
	segments(x0 = p63_sox2_pik3ca_max_x + 5, y0 = p63_sox2_pik3ca_min_y, x1 = p63_sox2_pik3ca_max_x + 5, y1 = p63_sox2_pik3ca_max_y) #vertical line
	segments(x0 = p63_sox2_pik3ca_max_x + 5, y0 = p63_sox2_pik3ca_max_y, x1 = p63_sox2_pik3ca_max_x + 2.5, y1 = p63_sox2_pik3ca_max_y) #bottom horizontal line
	segments(x0 = p63_sox2_pik3ca_max_x + 5, y0 = (p63_sox2_pik3ca_min_y + p63_sox2_pik3ca_max_y) / 2, x1 = p63_sox2_pik3ca_max_x + 7.5, y1 = (p63_sox2_pik3ca_min_y + p63_sox2_pik3ca_max_y) / 2) #midpoint horizontal line
	text("TP63/SOX2/PIK3CA", x = p63_sox2_pik3ca_max_x + 7.5, y = (p63_sox2_pik3ca_min_y + p63_sox2_pik3ca_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	
	nrf2_locations = subset(plot_gains_df, str_detect(plot_gains_df$Genes, " NFE2L2,"))
	nrf2_max_x  = max(nrf2_locations$x_values)
	nrf2_min_y = min(nrf2_locations$y_values)
	nrf2_max_y = max(nrf2_locations$y_values)
	
	segments(x0 = nrf2_max_x + 2.5, y0 = nrf2_min_y, x1 = nrf2_max_x + 5, y1 = nrf2_min_y) #top horizontal line
	segments(x0 = nrf2_max_x + 5, y0 = nrf2_min_y, x1 = nrf2_max_x + 5, y1 = nrf2_max_y) #vertical line
	segments(x0 = nrf2_max_x + 5, y0 = nrf2_max_y, x1 = nrf2_max_x + 2.5, y1 = nrf2_max_y) #bottom horizontal line
	segments(x0 = nrf2_max_x + 5, y0 = (nrf2_min_y + nrf2_max_y) / 2, x1 = nrf2_max_x + 7.5, y1 = (nrf2_min_y + nrf2_max_y) / 2) #midpoint horizontal line
	text("NFE2L2", x = nrf2_max_x + 7.5, y = (nrf2_min_y + nrf2_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	
	cdk6_locations = subset(plot_gains_df, str_detect(plot_gains_df$Genes, " CDK6,"))
	cdk6_max_x  = max(cdk6_locations$x_values)
	cdk6_min_y = min(cdk6_locations$y_values)
	cdk6_max_y = max(cdk6_locations$y_values)
	
	segments(x0 = cdk6_max_x + 2.5, y0 = cdk6_min_y, x1 = cdk6_max_x + 5, y1 = cdk6_min_y) #top horizontal line
	segments(x0 = cdk6_max_x + 5, y0 = cdk6_min_y, x1 = cdk6_max_x + 5, y1 = cdk6_max_y) #vertical line
	segments(x0 = cdk6_max_x + 5, y0 = cdk6_max_y, x1 = cdk6_max_x + 2.5, y1 = cdk6_max_y) #bottom horizontal line
	segments(x0 = cdk6_max_x + 5, y0 = (cdk6_min_y + cdk6_max_y) / 2, x1 = cdk6_max_x + 7.5, y1 = (cdk6_min_y + cdk6_max_y) / 2) #midpoint horizontal line
	text("CDK6", x = cdk6_max_x + 7.5, y = (cdk6_min_y + cdk6_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	
	bcl2l1_locations = subset(plot_gains_df, str_detect(plot_gains_df$Genes, " BCL2L1,"))
	bcl2l1_max_x  = max(bcl2l1_locations$x_values)
	bcl2l1_min_y = min(bcl2l1_locations$y_values)
	bcl2l1_max_y = max(bcl2l1_locations$y_values)
	
	segments(x0 = bcl2l1_max_x + 2.5, y0 = bcl2l1_min_y, x1 = bcl2l1_max_x + 5, y1 = bcl2l1_min_y) #top horizontal line
	segments(x0 = bcl2l1_max_x + 5, y0 = bcl2l1_min_y, x1 = bcl2l1_max_x + 5, y1 = bcl2l1_max_y) #vertical line
	segments(x0 = bcl2l1_max_x + 5, y0 = bcl2l1_max_y, x1 = bcl2l1_max_x + 2.5, y1 = bcl2l1_max_y) #bottom horizontal line
	segments(x0 = bcl2l1_max_x + 5, y0 = (bcl2l1_min_y + bcl2l1_max_y) / 2, x1 = bcl2l1_max_x + 7.5, y1 = (bcl2l1_min_y + bcl2l1_max_y) / 2) #midpoint horizontal line
	text("BCL2L1", x = bcl2l1_max_x + 7.5, y = (bcl2l1_min_y + bcl2l1_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	
	
	ccnd_locations = subset(plot_gains_df, str_detect(plot_gains_df$Genes, " CCND1,"))
	ccnd_max_x  = max(ccnd_locations$x_values)
	ccnd_min_y = min(ccnd_locations$y_values)
	ccnd_max_y = max(ccnd_locations$y_values)
	
	segments(x0 = ccnd_max_x + 2.5, y0 = ccnd_min_y, x1 = ccnd_max_x + 5, y1 = ccnd_min_y) #top horizontal line
	segments(x0 = ccnd_max_x + 5, y0 = ccnd_min_y, x1 = ccnd_max_x + 5, y1 = ccnd_max_y) #vertical line
	segments(x0 = ccnd_max_x + 5, y0 = ccnd_max_y, x1 = ccnd_max_x + 2.5, y1 = ccnd_max_y) #bottom horizontal line
	segments(x0 = ccnd_max_x + 5, y0 = (ccnd_min_y + ccnd_max_y) / 2, x1 = ccnd_max_x + 7.5, y1 = (ccnd_min_y + ccnd_max_y) / 2) #midpoint horizontal line
	text("CCND1", x = ccnd_max_x + 7.5, y = (ccnd_min_y + ccnd_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	PDGFRA_locations = subset(plot_gains_df, str_detect(plot_gains_df$Genes, " PDGFRA,|^PDGFRA, |, PDGFRA$"))
	PDGFRA_max_x  = max(PDGFRA_locations$x_values)
	PDGFRA_min_y = min(PDGFRA_locations$y_values)
	PDGFRA_max_y = max(PDGFRA_locations$y_values)
	
	segments(x0 = PDGFRA_max_x + 2.5, y0 = PDGFRA_min_y, x1 = PDGFRA_max_x + 5, y1 = PDGFRA_min_y) #top horizontal line
	segments(x0 = PDGFRA_max_x + 5, y0 = PDGFRA_min_y, x1 = PDGFRA_max_x + 5, y1 = PDGFRA_max_y) #vertical line
	segments(x0 = PDGFRA_max_x + 5, y0 = PDGFRA_max_y, x1 = PDGFRA_max_x + 2.5, y1 = PDGFRA_max_y) #bottom horizontal line
	segments(x0 = PDGFRA_max_x + 5, y0 = (PDGFRA_min_y + PDGFRA_max_y) / 2, x1 = PDGFRA_max_x + 7.5, y1 = (PDGFRA_min_y + PDGFRA_max_y) / 2) #midpoint horizontal line
	text("PDGFRA", x = PDGFRA_max_x + 7.5, y = (PDGFRA_min_y + PDGFRA_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	
	EGFR_locations = subset(plot_gains_df, str_detect(plot_gains_df$Genes, " EGFR,|^EGFR, |, EGFR$"))
	EGFR_max_x  = max(EGFR_locations$x_values) + 15
	EGFR_min_y = min(EGFR_locations$y_values)
	EGFR_max_y = max(EGFR_locations$y_values)
	
	segments(x0 = EGFR_max_x + 2.5, y0 = EGFR_min_y, x1 = EGFR_max_x + 5, y1 = EGFR_min_y) #top horizontal line
	segments(x0 = EGFR_max_x + 5, y0 = EGFR_min_y, x1 = EGFR_max_x + 5, y1 = EGFR_max_y) #vertical line
	segments(x0 = EGFR_max_x + 5, y0 = EGFR_max_y, x1 = EGFR_max_x + 2.5, y1 = EGFR_max_y) #bottom horizontal line
	segments(x0 = EGFR_max_x + 5, y0 = (EGFR_min_y + EGFR_max_y) / 2, x1 = EGFR_max_x + 7.5, y1 = (EGFR_min_y + EGFR_max_y) / 2) #midpoint horizontal line
	text("EGFR", x = EGFR_max_x + 7.5, y = (EGFR_min_y + EGFR_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	
	FGFR1_locations = subset(plot_gains_df, str_detect(plot_gains_df$Genes, " FGFR1,|^FGFR1, |, FGFR1$"))
	FGFR1_max_x  = max(FGFR1_locations$x_values) + 10
	FGFR1_min_y = min(FGFR1_locations$y_values)
	FGFR1_max_y = max(FGFR1_locations$y_values)
	
	segments(x0 = FGFR1_max_x + 2.5, y0 = FGFR1_min_y, x1 = FGFR1_max_x + 5, y1 = FGFR1_min_y) #top horizontal line
	segments(x0 = FGFR1_max_x + 5, y0 = FGFR1_min_y, x1 = FGFR1_max_x + 5, y1 = FGFR1_max_y) #vertical line
	segments(x0 = FGFR1_max_x + 5, y0 = FGFR1_max_y, x1 = FGFR1_max_x + 2.5, y1 = FGFR1_max_y) #bottom horizontal line
	segments(x0 = FGFR1_max_x + 5, y0 = (FGFR1_min_y + FGFR1_max_y) / 2, x1 = FGFR1_max_x + 7.5, y1 = (FGFR1_min_y + FGFR1_max_y) / 2) #midpoint horizontal line
	text("FGFR1", x = FGFR1_max_x + 7.5, y = (FGFR1_min_y + FGFR1_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	
	
	chrom01_q_locations = subset(plot_gains_df, str_detect(plot_gains_df$chromosome_band, "01\\_q"))
	chrom01_q_max_x  = max(chrom01_q_locations$x_values)
	chrom01_q_min_y = min(chrom01_q_locations$y_values)
	chrom01_q_max_y = max(chrom01_q_locations$y_values)
	
	segments(x0 = chrom01_q_max_x + 2.5, y0 = chrom01_q_min_y, x1 = chrom01_q_max_x + 5, y1 = chrom01_q_min_y) #top horizontal line
	segments(x0 = chrom01_q_max_x + 5, y0 = chrom01_q_min_y, x1 = chrom01_q_max_x + 5, y1 = chrom01_q_max_y) #vertical line
	segments(x0 = chrom01_q_max_x + 5, y0 = chrom01_q_max_y, x1 = chrom01_q_max_x + 2.5, y1 = chrom01_q_max_y) #bottom horizontal line
	segments(x0 = chrom01_q_max_x + 5, y0 = (chrom01_q_min_y + chrom01_q_max_y) / 2, x1 = chrom01_q_max_x + 7.5, y1 = (chrom01_q_min_y + chrom01_q_max_y) / 2) #midpoint horizontal line
	text("1q", x = chrom01_q_max_x + 7.5, y = (chrom01_q_min_y + chrom01_q_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	
	chrom05_p15.33_locations = subset(plot_gains_df, str_detect(plot_gains_df$chromosome_band, fixed("05_p15.33")))
	chrom05_p15.33_max_x  = max(chrom05_p15.33_locations$x_values)
	chrom05_p15.33_min_y = min(chrom05_p15.33_locations$y_values)
	chrom05_p15.33_max_y = max(chrom05_p15.33_locations$y_values)
	
	segments(x0 = chrom05_p15.33_max_x + 2.5, y0 = chrom05_p15.33_min_y, x1 = chrom05_p15.33_max_x + 5, y1 = chrom05_p15.33_min_y) #top horizontal line
	segments(x0 = chrom05_p15.33_max_x + 5, y0 = chrom05_p15.33_min_y, x1 = chrom05_p15.33_max_x + 5, y1 = chrom05_p15.33_max_y) #vertical line
	segments(x0 = chrom05_p15.33_max_x + 5, y0 = chrom05_p15.33_max_y, x1 = chrom05_p15.33_max_x + 2.5, y1 = chrom05_p15.33_max_y) #bottom horizontal line
	segments(x0 = chrom05_p15.33_max_x + 5, y0 = (chrom05_p15.33_min_y + chrom05_p15.33_max_y) / 2, x1 = chrom05_p15.33_max_x + 7.5, y1 = (chrom05_p15.33_min_y + chrom05_p15.33_max_y) / 2) #midpoint horizontal line
	text("5p15.33", x = chrom05_p15.33_max_x + 7.5, y = (chrom05_p15.33_min_y + chrom05_p15.33_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	
	chrom08_q_locations = subset(plot_gains_df, str_detect(plot_gains_df$chromosome_band, "08\\_q"))
	chrom08_q_max_x  = max(chrom08_q_locations$x_values)
	chrom08_q_min_y = min(chrom08_q_locations$y_values)
	chrom08_q_max_y = max(chrom08_q_locations$y_values)
	
	segments(x0 = chrom08_q_max_x + 2.5, y0 = chrom08_q_min_y, x1 = chrom08_q_max_x + 5, y1 = chrom08_q_min_y) #top horizontal line
	segments(x0 = chrom08_q_max_x + 5, y0 = chrom08_q_min_y, x1 = chrom08_q_max_x + 5, y1 = chrom08_q_max_y) #vertical line
	segments(x0 = chrom08_q_max_x + 5, y0 = chrom08_q_max_y, x1 = chrom08_q_max_x + 2.5, y1 = chrom08_q_max_y) #bottom horizontal line
	segments(x0 = chrom08_q_max_x + 5, y0 = (chrom08_q_min_y + chrom08_q_max_y) / 2, x1 = chrom08_q_max_x + 7.5, y1 = (chrom08_q_min_y + chrom08_q_max_y) / 2) #midpoint horizontal line
	text("8q", x = chrom08_q_max_x + 7.5, y = (chrom08_q_min_y + chrom08_q_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	chrom14_q32.33_locations = subset(plot_gains_df, str_detect(plot_gains_df$chromosome_band, "14\\_q32\\.33"))
	chrom14_q32.33_max_x  = max(chrom14_q32.33_locations$x_values)
	chrom14_q32.33_min_y = min(chrom14_q32.33_locations$y_values)
	chrom14_q32.33_max_y = max(chrom14_q32.33_locations$y_values)
	
	segments(x0 = chrom14_q32.33_max_x + 2.5, y0 = chrom14_q32.33_min_y, x1 = chrom14_q32.33_max_x + 5, y1 = chrom14_q32.33_min_y) #top horizontal line
	segments(x0 = chrom14_q32.33_max_x + 5, y0 = chrom14_q32.33_min_y, x1 = chrom14_q32.33_max_x + 5, y1 = chrom14_q32.33_max_y) #vertical line
	segments(x0 = chrom14_q32.33_max_x + 5, y0 = chrom14_q32.33_max_y, x1 = chrom14_q32.33_max_x + 2.5, y1 = chrom14_q32.33_max_y) #bottom horizontal line
	segments(x0 = chrom14_q32.33_max_x + 5, y0 = (chrom14_q32.33_min_y + chrom14_q32.33_max_y) / 2, x1 = chrom14_q32.33_max_x + 7.5, y1 = (chrom14_q32.33_min_y + chrom14_q32.33_max_y) / 2) #midpoint horizontal line
	text("14q32.33", x = chrom14_q32.33_max_x + 7.5, y = (chrom14_q32.33_min_y + chrom14_q32.33_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
	
	chrom17_q21.31_locations = subset(plot_gains_df, str_detect(plot_gains_df$chromosome_band, "17\\_q21\\.31"))
	chrom17_q21.31_max_x  = max(chrom17_q21.31_locations$x_values)
	chrom17_q21.31_min_y = min(chrom17_q21.31_locations$y_values)
	chrom17_q21.31_max_y = max(chrom17_q21.31_locations$y_values)
	
	segments(x0 = chrom17_q21.31_max_x + 2.5, y0 = chrom17_q21.31_min_y, x1 = chrom17_q21.31_max_x + 5, y1 = chrom17_q21.31_min_y) #top horizontal line
	segments(x0 = chrom17_q21.31_max_x + 5, y0 = chrom17_q21.31_min_y, x1 = chrom17_q21.31_max_x + 5, y1 = chrom17_q21.31_max_y) #vertical line
	segments(x0 = chrom17_q21.31_max_x + 5, y0 = chrom17_q21.31_max_y, x1 = chrom17_q21.31_max_x + 2.5, y1 = chrom17_q21.31_max_y) #bottom horizontal line
	segments(x0 = chrom17_q21.31_max_x + 5, y0 = (chrom17_q21.31_min_y + chrom17_q21.31_max_y) / 2, x1 = chrom17_q21.31_max_x + 7.5, y1 = (chrom17_q21.31_min_y + chrom17_q21.31_max_y) / 2) #midpoint horizontal line
	text("17q21.31", x = chrom17_q21.31_max_x + 7.5, y = (chrom17_q21.31_min_y + chrom17_q21.31_max_y) / 2, pos = 4, offset = 0, cex = cex_text_parameter)
}

mtext("Chromosome", side = 3, line = -6, outer = TRUE)
#mtext("Alteration Percentage", side = 3, line = -53, outer = TRUE)
dev.off()


#reviewer comments - summarize cnv
cnv_long_gain_unique = unique(cnv_long_cnstate_gain[, c("TissueID", "chromosome_band")])
cnv_long_gain_unique_counts = data.frame(sort(table(cnv_long_gain_unique$chromosome_band), decreasing = TRUE))
cnv_long_gain_unique_counts$Type = "Gain"
cnv_long_loss_unique = unique(cnv_long_cnstate_loss[, c("TissueID", "chromosome_band")])
cnv_long_loss_unique_counts = data.frame(sort(table(cnv_long_loss_unique$chromosome_band), decreasing = TRUE))
cnv_long_loss_unique_counts$Type = "Loss"

redox_patients = subset(tmt_meta, subtype == "Redox")$TissueID
cnv_long_gain_unique_redox = unique(subset(cnv_long_cnstate_gain, TissueID %in% redox_patients)[, c("TissueID", "chromosome_band")])
cnv_long_gain_unique_redox_counts = data.frame(sort(table(cnv_long_gain_unique_redox$chromosome_band), decreasing = TRUE))
cnv_long_gain_unique_redox_counts$Type = "Gain"
cnv_long_loss_unique_redox = unique(subset(cnv_long_cnstate_loss, TissueID %in% redox_patients)[, c("TissueID", "chromosome_band")])
cnv_long_loss_unique_redox_counts = data.frame(sort(table(cnv_long_loss_unique_redox$chromosome_band), decreasing = TRUE))
cnv_long_loss_unique_redox_counts$Type = "Loss"

inflamed_patients = subset(tmt_meta, subtype == "Inflamed")$TissueID
cnv_long_gain_unique_inflamed = unique(subset(cnv_long_cnstate_gain, TissueID %in% inflamed_patients)[, c("TissueID", "chromosome_band")])
cnv_long_gain_unique_inflamed_counts = data.frame(sort(table(cnv_long_gain_unique_inflamed$chromosome_band), decreasing = TRUE))
cnv_long_gain_unique_inflamed_counts$Type = "Gain"
cnv_long_loss_unique_inflamed = unique(subset(cnv_long_cnstate_loss, TissueID %in% inflamed_patients)[, c("TissueID", "chromosome_band")])
cnv_long_loss_unique_inflamed_counts = data.frame(sort(table(cnv_long_loss_unique_inflamed$chromosome_band), decreasing = TRUE))
cnv_long_loss_unique_inflamed_counts$Type = "Loss"

mixed_patients = subset(tmt_meta, subtype == "Mixed")$TissueID
cnv_long_gain_unique_mixed = unique(subset(cnv_long_cnstate_gain, TissueID %in% mixed_patients)[, c("TissueID", "chromosome_band")])
cnv_long_gain_unique_mixed_counts = data.frame(sort(table(cnv_long_gain_unique_mixed$chromosome_band), decreasing = TRUE))
cnv_long_gain_unique_mixed_counts$Type = "Gain"
cnv_long_loss_unique_mixed = unique(subset(cnv_long_cnstate_loss, TissueID %in% mixed_patients)[, c("TissueID", "chromosome_band")])
cnv_long_loss_unique_mixed_counts = data.frame(sort(table(cnv_long_loss_unique_mixed$chromosome_band), decreasing = TRUE))
cnv_long_loss_unique_mixed_counts$Type = "Loss"

sum(cnv_long_gain_unique_redox_counts$Freq) + sum(cnv_long_loss_unique_redox_counts$Freq)
sum(cnv_long_gain_unique_inflamed_counts$Freq) + sum(cnv_long_loss_unique_inflamed_counts$Freq)
sum(cnv_long_gain_unique_mixed_counts$Freq) + sum(cnv_long_loss_unique_mixed_counts$Freq)

sum(cnv_long_gain_unique_redox_counts$Freq) + sum(cnv_long_gain_unique_inflamed_counts$Freq) + sum(cnv_long_gain_unique_mixed_counts$Freq)
sum(cnv_long_loss_unique_redox_counts$Freq) + sum(cnv_long_loss_unique_inflamed_counts$Freq) + sum(cnv_long_loss_unique_mixed_counts$Freq)


length(cnv_long_gain_unique_counts$Var1)
length(cnv_long_loss_unique_counts$Var1)

#copied from TCGA supplemental 2018-03-06
tcga_gistic_amplifications = c("03_q26.33", "08_p11.23", "011_q13.3", "02_p16.1", "04_q12", "08_q24.21", "07_p11.2", "01_q21.2", "09_p21.1", "015_q26.3", "017_q25.1", "020_q11.21", "02_q31.2", "06_q12", "014_q13.3", "019_q13.13", "07_q21.2", "05_p15.33", "021_q21.1", "017_p11.2", "022_q11.21", "03_p11.1", "014_q32.31", "02_q11.2", "011_p11.2", "06_q22.31", "09_p24.1", "012_q15", "05_p13.1", "012_p13.33", "018_q11.2", "01_p34.2") 
tcga_gistic_losses = c("9_p21.3", "08_p23.2", "02_q22.1", "010_q23.31", "05_q11.2", "04_q22.1", "03_p13", "011_p15.5", "02_q37.3", "04_q35.2", "022_q11.21", "01_p13.1", "011_q25", "021_q21.1", "019_p13.3", "03_p25.3", "019_q13.33", "017_q11.2", "03_p12.2", "04_q32.2", "016_q23.1", "01_p36.12", "09_p24.1", "09_p13.1", "015_q11.2", "06_q27", "018_q23", "016_p13.3", "07_q11.22", "013_q12.11", "014_q21.1", "05_q31.1", "07_q36.1", "05_q35.3", "03_p14.2", "010_p15.3", "04_q25", "05_p12", "03_q13.31", "010_q26.3", "04_p16.3")

table(tcga_gistic_amplifications %in% cnv_long_gain_unique_counts$Var1)
table(tcga_gistic_losses %in% cnv_long_loss_unique_counts$Var1)

#finding genomic locations and counting number of amplifications ----

#NFE2L2
cnv_long_cnstate_gain_bandfreq[str_detect(cnv_long_cnstate_gain_bandfreq$Genes, "NFE2L2"), 1:3] #genomic location of NRF2?
table(subset(tmt_meta, subtype == "Redox")$TissueID %in% subset(cnv_long, chromosome_band == "02_q23.3" | chromosome_band == "02_q24.1" | chromosome_band == "02_q24.2" | chromosome_band == "02_q31.1")$TissueID) #how many redox patients with NRF2 amplified?
redox_table = data.frame(TissueID = subset(tmt_meta, subtype == "Redox")$TissueID, stringsAsFactors = FALSE)
redox_table$NRF2_amplified = redox_table$TissueID %in% subset(cnv_long, chromosome_band == "02_q23.3" | chromosome_band == "02_q24.1" | chromosome_band == "02_q24.2" | chromosome_band == "02_q31.1")$TissueID
redox_table$KEAP1_loss = redox_table$TissueID %in% subset(cnv_long_cnstate_loss, chromosome_band == "19_p13.3")$TissueID
redox_table$NRF2_mutated = redox_table$TissueID %in% subset(tmt_meta, subtype == "Redox" & NFE2L2 != "wt")$TissueID
redox_table$KEAP1_mutated = redox_table$TissueID %in% subset(tmt_meta, subtype == "Redox" & KEAP1 != "wt")$TissueID
redox_table$at_least_one_true = rowSums(redox_table[, c("NRF2_amplified", "NRF2_mutated", "KEAP1_mutated", "KEAP1_loss")])
nrf2_keap1_venn = list("NRF2 MUT" = subset(redox_table, NRF2_mutated == "TRUE")$TissueID, "NRF2 AMP" = subset(redox_table, NRF2_amplified == "TRUE")$TissueID, "KEAP1 MUT" = subset(redox_table, KEAP1_mutated == "TRUE")$TissueID, "KEAP1 LOSS" = subset(redox_table, KEAP1_loss == "TRUE")$TissueID)
venn.diagram(nrf2_keap1_venn, "output/main/nrf2_keap1_venn_2018-12-06.png", cex = 2)

nrow(subset(redox_table, redox_table$NRF2_amplified == TRUE & redox_table$NRF2_mutated == TRUE | redox_table$NRF2_amplified == TRUE & redox_table$KEAP1_mutated == TRUE))
nrow(subset(redox_table, redox_table$NRF2_amplified == TRUE & redox_table$NRF2_mutated != TRUE | redox_table$NRF2_amplified == TRUE & redox_table$KEAP1_mutated != TRUE))
nrow(subset(redox_table, redox_table$NRF2_mutated != TRUE & redox_table$KEAP1_mutated != TRUE))


#redox group KEAP1 loss 
cnv_long_cnstate_loss_bandfreq[str_detect(cnv_long_cnstate_loss_bandfreq$Genes, "KEAP1"), 1:3] #genomic location of KEAP1?
table(subset(tmt_meta, subtype == "Redox")$TissueID %in% subset(cnv_long_cnstate_loss, chromosome_band == "19_p13.3")$TissueID)

#BCL2L1
cnv_long_cnstate_gain_bandfreq[str_detect(cnv_long_cnstate_gain_bandfreq$Genes, " BCL2L1, "), 1:3] #genomic location of BCL2L1?
table(tmt_meta$TissueID %in% subset(cnv_long, chromosome_band == "20_p11.21" | chromosome_band == "20_p12.1" | chromosome_band == "20_p12.2" | chromosome_band == "20_p13" | chromosome_band == "20_q11.21")$TissueID) #how many redox patients with BCL2L1 amplified?



#for manuscript heatmap

tmt_meta$amp_3q2 = tmt_meta$TissueID %in% cnv_long[str_detect(cnv_long$chromosome_band, "^03\\_q2[0-9]*"), ]$TissueID
tmt_meta$amp_nfe2l2 = tmt_meta$TissueID %in% subset(cnv_long, chromosome_band == "02_q23.3" | chromosome_band == "02_q24.1" | chromosome_band == "02_q24.2" | chromosome_band == "02_q31.1")$TissueID
tmt_meta$loss_keap1 = tmt_meta$TissueID %in% subset(cnv_long_cnstate_loss, chromosome_band == "19_p13.3")$TissueID
tmt_meta$amp_3q2 = str_replace(tmt_meta$amp_3q2, "TRUE", "#800026")
tmt_meta$amp_3q2 = str_replace(tmt_meta$amp_3q2, "FALSE", "lightgrey") 
tmt_meta$amp_nfe2l2 = str_replace(tmt_meta$amp_nfe2l2, "TRUE", "#800026")
tmt_meta$amp_nfe2l2 = str_replace(tmt_meta$amp_nfe2l2, "FALSE", "lightgrey")
tmt_meta$loss_keap1 = str_replace(tmt_meta$loss_keap1, "TRUE", "#08306B")
tmt_meta$loss_keap1 = str_replace(tmt_meta$loss_keap1, "FALSE", "lightgrey")
copy_number_heatmap_coloring = as.matrix(tmt_meta[, c("loss_keap1", "amp_nfe2l2", "amp_3q2")])
colnames(copy_number_heatmap_coloring) = c("KEAP1 Loss", "NFE2L2 Amp.","3q2 Amp.")


#DRIVE results  ----
drive_results = read.delim("data/main/drive_rsa_lusq_vs_all_short_list_clustered_remove_header_2018-02-14.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
drive_meta = read.delim("data/main/drive_rsa_lusq_vs_all_short_list_clustered_meta_2018-02-14.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
drive_colors = brewer.pal(3, "Set3")

drive_sqlc_cells = subset(drive_meta, subtype == "LUSQ")$column_name
drive_adc_cells = subset(drive_meta, subtype == "LUAD")$column_name

drive_sqlc_tp63 = data.frame(cell_line = drive_sqlc_cells, tp63 = as.numeric(drive_results[drive_results$Symbol == "TP63", drive_sqlc_cells]), color = drive_colors[1], stringsAsFactors = FALSE)
drive_sqlc_tp63 = drive_sqlc_tp63[order(drive_sqlc_tp63$tp63), ]
drive_adc_tp63 = data.frame(cell_line = drive_adc_cells, tp63 = as.numeric(drive_results[drive_results$Symbol == "TP63", drive_adc_cells]), color = drive_colors[2], stringsAsFactors = FALSE)
drive_adc_tp63 = drive_adc_tp63[order(drive_adc_tp63$tp63), ]


#tp63 adc and sqlc
par(mar = c(5.1, 4.1, 7.1, 2.1))
drive_plot_tp63 = rbind(drive_sqlc_tp63, drive_adc_tp63)
drive_plot_tp63$cell_line = str_replace(drive_plot_tp63$cell_line, "\\_lung", "")

tp63_barplot = barplot(drive_plot_tp63$tp63, xaxt = "n", ylim = c(round(min(drive_plot_tp63$tp63)), 0), col = drive_plot_tp63$color, ylab = "ATARiS Score")
axis(3, at = tp63_barplot, labels = drive_plot_tp63$cell_line, las = 2)
abline(h = -3, lty = 2)
legend("bottomright", legend = c("SQLC", "ADC"), fill = c(drive_colors[1], drive_colors[2]))

#tp63 sqlc alone
par(mar = c(5.1, 4.1, 10, 2.1))
par(mfrow = c(1,3))
drive_plot_tp63 = drive_sqlc_tp63
drive_plot_tp63$cell_line = str_replace(drive_plot_tp63$cell_line, "\\_lung", "")

tp63_barplot = barplot(drive_plot_tp63$tp63, xaxt = "n", ylim = c(-6, 0), col = "red", ylab = "ATARiS Score", main = "TP63")
axis(3, at = tp63_barplot, labels = drive_plot_tp63$cell_line, las = 2)
abline(h = -3, lty = 2)


#PSAT1 SQLC alone
drive_sqlc_PSAT1 = data.frame(cell_line = drive_sqlc_cells, PSAT1 = as.numeric(drive_results[drive_results$Symbol == "PSAT1", drive_sqlc_cells]), color = drive_colors[1], stringsAsFactors = FALSE)
drive_sqlc_PSAT1 = drive_sqlc_PSAT1[order(drive_sqlc_PSAT1$PSAT1), ]
drive_plot_PSAT1 = drive_sqlc_PSAT1
drive_plot_PSAT1$cell_line = str_replace(drive_plot_PSAT1$cell_line, "\\_lung", "")

PSAT1_barplot = barplot(drive_plot_PSAT1$PSAT1, xaxt = "n", ylim = c(-6, 0), col = "red", ylab = "ATARiS Score", main = "PSAT1")
axis(3, at = PSAT1_barplot, labels = drive_plot_PSAT1$cell_line, las = 2)
abline(h = -3, lty = 2)

#AKR1C3 SQLC alone
drive_sqlc_AKR1C3 = data.frame(cell_line = drive_sqlc_cells, AKR1C3 = as.numeric(drive_results[drive_results$Symbol == "AKR1C3", drive_sqlc_cells]), color = drive_colors[1], stringsAsFactors = FALSE)
drive_sqlc_AKR1C3 = drive_sqlc_AKR1C3[order(drive_sqlc_AKR1C3$AKR1C3), ]
drive_plot_AKR1C3 = drive_sqlc_AKR1C3
drive_plot_AKR1C3$cell_line = str_replace(drive_plot_AKR1C3$cell_line, "\\_lung", "")

AKR1C3_barplot = barplot(drive_plot_AKR1C3$AKR1C3, xaxt = "n", ylim = c(-6, 0), col = "red", ylab = "ATARiS Score", main = "AKR1C3")
axis(3, at = AKR1C3_barplot, labels = drive_plot_AKR1C3$cell_line, las = 2)
abline(h = -3, lty = 2)



#poster & manuscript
drive_sqlc_cells
#999 solid, 135 hashed
drive_sqlc_nrf2_keap1_dens_for_plot = c(999, 999, 25, 25, 999, 25, 25, 25, 25)
png(filename = "output/main/scc_redox_vulnerabilities_2018-11-28.png", width = 1200, height = 1000)
par(mar = c(5.1, 6.1, 12, 2.1))
par(mfrow = c(2,3))

drive_sqlc_TP63 = data.frame(cell_line = drive_sqlc_cells, TP63 = as.numeric(drive_results[drive_results$Symbol == "TP63", drive_sqlc_cells]), color = drive_colors[1], density = drive_sqlc_nrf2_keap1_dens_for_plot, stringsAsFactors = FALSE)
drive_sqlc_TP63 = drive_sqlc_TP63[order(drive_sqlc_TP63$TP63), ]
drive_plot_TP63 = drive_sqlc_TP63
drive_plot_TP63$cell_line = str_replace(drive_plot_TP63$cell_line, "\\_lung", "")


TP63_barplot = barplot(drive_plot_TP63$TP63, xaxt = "n", ylim = c(-12, 0), col = "#FF000A", ylab = "RSA", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = drive_plot_TP63$density)
axis(3, at = TP63_barplot, labels = drive_plot_TP63$cell_line, las = 2, cex.axis = 2)
abline(h = -3, lty = 2)

drive_sqlc_PSAT1 = data.frame(cell_line = drive_sqlc_cells, PSAT1 = as.numeric(drive_results[drive_results$Symbol == "PSAT1", drive_sqlc_cells]), color = drive_colors[1], density = drive_sqlc_nrf2_keap1_dens_for_plot, stringsAsFactors = FALSE)

drive_sqlc_PSAT1 = drive_sqlc_PSAT1[order(drive_sqlc_PSAT1$PSAT1), ]
drive_plot_PSAT1 = drive_sqlc_PSAT1
drive_plot_PSAT1$cell_line = str_replace(drive_plot_PSAT1$cell_line, "\\_lung", "")

PSAT1_barplot = barplot(drive_plot_PSAT1$PSAT1, xaxt = "n", ylim = c(-12, 0), col = "#FF000A", ylab = "RSA", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = drive_plot_PSAT1$density)

axis(3, at = PSAT1_barplot, labels = drive_plot_PSAT1$cell_line, las = 2, cex.axis = 2)
abline(h = -3, lty = 2)

drive_sqlc_AKR1C3 = data.frame(cell_line = drive_sqlc_cells, AKR1C3 = as.numeric(drive_results[drive_results$Symbol == "AKR1C3", drive_sqlc_cells]), color = drive_colors[1], density = drive_sqlc_nrf2_keap1_dens_for_plot, stringsAsFactors = FALSE)

drive_sqlc_AKR1C3 = drive_sqlc_AKR1C3[order(drive_sqlc_AKR1C3$AKR1C3), ]
drive_plot_AKR1C3 = drive_sqlc_AKR1C3
drive_plot_AKR1C3$cell_line = str_replace(drive_plot_AKR1C3$cell_line, "\\_lung", "")

AKR1C3_barplot = barplot(drive_plot_AKR1C3$AKR1C3, xaxt = "n", ylim = c(-12, 0), col = "#FF000A", ylab = "RSA", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = drive_plot_AKR1C3$density)

axis(3, at = AKR1C3_barplot, labels = drive_plot_AKR1C3$cell_line, las = 2, cex.axis = 2)
abline(h = -3, lty = 2)


drive_sqlc_GSR = data.frame(cell_line = drive_sqlc_cells, GSR = as.numeric(drive_results[drive_results$Symbol == "GSR", drive_sqlc_cells]), color = drive_colors[1], density = drive_sqlc_nrf2_keap1_dens_for_plot, stringsAsFactors = FALSE)

drive_sqlc_GSR = drive_sqlc_GSR[order(drive_sqlc_GSR$GSR), ]
drive_plot_GSR = drive_sqlc_GSR
drive_plot_GSR$cell_line = str_replace(drive_plot_GSR$cell_line, "\\_lung", "")

GSR_barplot = barplot(drive_plot_GSR$GSR, xaxt = "n", ylim = c(-12, 0), col = "#CB6527", ylab = "RSA", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = drive_plot_GSR$density)

axis(3, at = GSR_barplot, labels = drive_plot_GSR$cell_line, las = 2, cex.axis = 2)
abline(h = -3, lty = 2)


drive_sqlc_TFRC = data.frame(cell_line = drive_sqlc_cells, TFRC = as.numeric(drive_results[drive_results$Symbol == "TFRC", drive_sqlc_cells]), color = drive_colors[1], density = drive_sqlc_nrf2_keap1_dens_for_plot, stringsAsFactors = FALSE)

drive_sqlc_TFRC = drive_sqlc_TFRC[order(drive_sqlc_TFRC$TFRC), ]
drive_plot_TFRC = drive_sqlc_TFRC
drive_plot_TFRC$cell_line = str_replace(drive_plot_TFRC$cell_line, "\\_lung", "")

TFRC_barplot = barplot(drive_plot_TFRC$TFRC, xaxt = "n", ylim = c(-12, 0), col = "#CB6527", ylab = "RSA", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = drive_plot_TFRC$density)

axis(3, at = TFRC_barplot, labels = drive_plot_TFRC$cell_line, las = 2, cex.axis = 2)
abline(h = -3, lty = 2)


drive_sqlc_SERPINB5 = data.frame(cell_line = drive_sqlc_cells, SERPINB5 = as.numeric(drive_results[drive_results$Symbol == "SERPINB5", drive_sqlc_cells]), color = drive_colors[1], density = drive_sqlc_nrf2_keap1_dens_for_plot, stringsAsFactors = FALSE)

drive_sqlc_SERPINB5 = drive_sqlc_SERPINB5[order(drive_sqlc_SERPINB5$SERPINB5), ]
drive_plot_SERPINB5 = drive_sqlc_SERPINB5
drive_plot_SERPINB5$cell_line = str_replace(drive_plot_SERPINB5$cell_line, "\\_lung", "")

SERPINB5_barplot = barplot(drive_plot_SERPINB5$SERPINB5, xaxt = "n", ylim = c(-12, 0), col = "#CB6527", ylab = "RSA", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = drive_plot_SERPINB5$density)

axis(3, at = SERPINB5_barplot, labels = drive_plot_SERPINB5$cell_line, las = 2, cex.axis = 2)
abline(h = -3, lty = 2)
dev.off()



drive_top_genes_enrichr = enrichr(genes = c("SMARCA4","DHX30","FABP4","TK1","SMC1A","ITGB1","RACGAP1","PI4KA,PI4KAP2,PI4KAP1","F13B","UBA3","SOD1","UBE2L3","KIF18A","TP63","CTHRC1","PRR5L","CDC7","ANGPTL4"), "GO_Biological_Process_2017")

png(filename = "output/main/pathway_enrichment/top_genes_drive_enrichr.png", width = 800, height = 250)
par(mar =  c(5.1, 38, 2, 2.1))
plot_enrichr_results(enrichr_in = drive_top_genes_enrichr$GO_Biological_Process_2017, number_of_pathways_to_plot = 10, xlim_in = c(0, 10), ylim_in = c(0, 10), plot_title = "")
dev.off()

drive_frac_pvals = read.delim("data/main/drive_ataris_sqlc_fractionsig_pvals_2018-03-13.txt", header = TRUE, stringsAsFactors = FALSE)
drive_frac_pvals$sqlc_adc_ratio = drive_frac_pvals$sqlc_fraction_sig / drive_frac_pvals$adc_fraction_sig
drive_frac_pvals$sqlc_other_ratio = drive_frac_pvals$sqlc_fraction_sig / drive_frac_pvals$other_fraction_sig
drive_frac_pvals = drive_frac_pvals[order(drive_frac_pvals$sqlc_adc_ratio, decreasing = TRUE), ]
drive_frac_pvals = subset(drive_frac_pvals, sqlc_fraction_sig > 0.12)
drive_pvals_sig = subset(drive_pvals, drive_pvals$sqlc_vs_adc_pval <= 0.001) #| drive_pvals$sqlc_vs_other_pval <= 0.05)
genes_drive_sig = drive_pvals_sig$Symbol
drive_pvals_sig_enrichr = enrichr(genes_drive_sig, "GO_Biological_Process_2017")

drive_frac_pvals = read.delim("data/main/drive_ataris_sqlc_fractionsig_pvals_2018-03-13.txt", header = TRUE, stringsAsFactors = FALSE)

drive_pval_cutoff = 0.05
drive_cutoff_genes = subset(drive_frac_pvals, sqlc_vs_adc_pval < drive_pval_cutoff | sqlc_vs_other_pval < drive_pval_cutoff)$Symbol



#CIBERSORT relative ----
#unlog and replace NA with 0 for uploading to CIBERSORT website
# rnaseq_for_cibersort = rnaseq_expression[, c("Symbol", tmt_meta$GenomicsSample)]
# rnaseq_for_cibersort[, tmt_meta$GenomicsSample] = 2^rnaseq_for_cibersort[, tmt_meta$GenomicsSample]
# rnaseq_for_cibersort[is.na(rnaseq_for_cibersort)] = 0
# write.table(rnaseq_for_cibersort, "output/main/rnaseq_for_cibersort_2018-05-14.txt", sep = "\t", quote = F, row.names = FALSE)

cibersort_rnaseq = read.delim("data/main/original_data/CIBERSORT_output_rna_for_cibersort_2018-05-15.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

pg_path = read.delim("data/main/SQLC116-HE_scores_20180111.txt", header = TRUE, stringsAsFactors = FALSE, row.name = 1)

cibersort_rnaseq = as.data.frame(t(subset(cibersort_rnaseq, P.value <= 0.05)))
cibersort_meta = subset(tmt_meta, GenomicsSample %in% names(cibersort_rnaseq))
cibersort_meta_path = merge(cibersort_meta, pg_path, by.x = "Deidentified.Patient.ID.x", by.y = "row.names")
cibersort_rnaseq = cibersort_rnaseq[!row.names(cibersort_rnaseq) %in% c("P.value", "Pearson.Correlation", "RMSE"), ]
cibersort_rnaseq = cibersort_rnaseq[, cibersort_meta$GenomicsSample]
#cibersort_rnaseq[cibersort_rnaseq == 0] = NA
cibersort_rnaseq$tln_correlation = apply(cibersort_rnaseq[, cibersort_meta_path$GenomicsSample], 1, function(x) {
	cor.test(x, cibersort_meta_path$TOL, method = "spearman", use = "pairwise.complete.obs")$estimate
})

cibersort_rnaseq$na_count = rowSums(is.na(cibersort_rnaseq[, cibersort_meta_path$GenomicsSample]))
#protein_heatmap_coloring_proteinonly = data.frame("Proteomics Group" = k5_protein_coloring[, 19], stringsAsFactors = FALSE, check.names = FALSE)
# 
# heatmap_obj = heatmap.4(cibersort_rnaseq[, cibersort_meta$GenomicsSample], Colv = FALSE, col = bluered, trace = "none", density.info = "none", main = "", scale = "row", dendrogram = "none", labCol = FALSE, labRow = row.names(cibersort_rnaseq), ColSideColors = as.matrix(protein_heatmap_coloring_proteinonly), ColSideColorsSize = 2, cexCol = 1.5, margins = c(2,12), na.color = "darkgrey", breaks = c(-2, seq(-2, 0, length = 50), seq(0, 2, length = 50)), xlab = "Patients", ylab = "",  distfun = function(c){as.dist(1 - cor(t(c), use = "pairwise.complete.obs"))},colsep = c(23, 43, 55, 94), sepcolor = "black")
# legend("bottomleft", legend = c("Inflamed", "Redox", "Mixed"), fill = c(darkcolors[1], darkcolors[2], darkcolors[3]), cex = 1)
# 



inflamed = subset(cibersort_meta, consensusClassK5_protein == 1 | consensusClassK5_protein == 3)$GenomicsSample
redoxgroup = subset(cibersort_meta, consensusClassK5_protein == 2 | consensusClassK5_protein == 4)$GenomicsSample
Mixedgroup = group3_proteincc_samples = subset(cibersort_meta, consensusClassK5_protein == 5)$GenomicsSample
inflamed_a = subset(cibersort_meta, consensusClassK5_protein == 1)$GenomicsSample
inflamed_b = subset(cibersort_meta, consensusClassK5_protein == 3)$GenomicsSample

par(mfrow = c(3,3))
for (i in 1:nrow(cibersort_rnaseq)) {
	boxplot(list("Inflamed" = as.numeric(cibersort_rnaseq[row.names(cibersort_rnaseq)[i], inflamed]), "Redox" = as.numeric(cibersort_rnaseq[row.names(cibersort_rnaseq)[i], redoxgroup]), "Double Neg" = as.numeric(cibersort_rnaseq[row.names(cibersort_rnaseq)[i], Mixedgroup])), main = row.names(cibersort_rnaseq)[i], ylab = "CIBERSORT Score")
}


cibersort_meta$subtype_numeric4 = NA
for (i in 1:nrow(cibersort_meta)) {
	if (cibersort_meta$consensusClassK5_protein[i] == 1) {
		cibersort_meta$subtype_numeric4[i] = 1 #inflamed A lump group
	}
	if (cibersort_meta$consensusClassK5_protein[i] == 3) {
		cibersort_meta$subtype_numeric4[i] = 2 #inflamed B lump group
	} 
	if (cibersort_meta$consensusClassK5_protein[i] == 2 | cibersort_meta$consensusClassK5_protein[i] == 4) {
		cibersort_meta$subtype_numeric4[i] = 3 #metabolism lump group
	} 
	if (cibersort_meta$consensusClassK5_protein[i] == 5) {
		cibersort_meta$subtype_numeric4[i] = 4 #double negative lump group
	} 
}


#wilcoxon test a given group vs the rest
cibersort_score_testing = data.frame(row.names = row.names(cibersort_rnaseq), Inflamed = rep(NA, nrow(cibersort_rnaseq)), Redox = rep(NA, nrow(cibersort_rnaseq)), Mixed = rep(NA, nrow(cibersort_rnaseq)), stringsAsFactors = FALSE)
for (i in 1:nrow(cibersort_rnaseq)) {
	for(j in 1:3) {

		tryCatch({
			current_subtype = subset(cibersort_meta, subtype_numeric == j)$GenomicsSample
			other_subtype = subset(cibersort_meta, subtype_numeric != j)$GenomicsSample
			cibersort_score_testing[i, j] = wilcox.test(as.numeric(cibersort_rnaseq[i, current_subtype]), as.numeric(cibersort_rnaseq[i, other_subtype]))$p.value
		}, error = function(e) {
			NA
		})
	}
}

cibersort_sig = cibersort_score_testing[rowSums(cibersort_score_testing < 0.05) > 0, ]
cibersort_rnaseq_sig = subset(cibersort_rnaseq, row.names(cibersort_rnaseq) %in% row.names(cibersort_sig))
par(mfrow = c(3,3))

for (i in 1:nrow(cibersort_rnaseq_sig)) {
	boxplot(list("Inflamed" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], inflamed]), "Redox" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], redoxgroup]), "Double Neg" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], Mixedgroup])), main = row.names(cibersort_rnaseq_sig)[i], ylab = "CIBERSORT Score")
}


#wilcoxon test inflamed a vs inflamed b
cibersort_score_testing_inflamed = data.frame(row.names = row.names(cibersort_rnaseq), InflamedAvsInflamedB = rep(NA, nrow(cibersort_rnaseq)), stringsAsFactors = FALSE)
for (i in 1:nrow(cibersort_rnaseq)) {
		tryCatch({
			current_subtype = subset(cibersort_meta, subtype_numeric4 == 1)$GenomicsSample
			other_subtype = subset(cibersort_meta, subtype_numeric4 == 2)$GenomicsSample
			cibersort_score_testing_inflamed[i, 1] = wilcox.test(as.numeric(cibersort_rnaseq[i, current_subtype]), as.numeric(cibersort_rnaseq[i, other_subtype]))$p.value
		}, error = function(e) {
			NA
		})
}

cibersort_sig_inflamed = row.names(cibersort_score_testing_inflamed)[rowSums(cibersort_score_testing_inflamed < 0.05) > 0]
cibersort_rnaseq_sig_inflamed = subset(cibersort_rnaseq, row.names(cibersort_rnaseq) %in% cibersort_sig_inflamed)
#cibersort_rnaseq_sig_inflamed = cibersort_rnaseq #all samples
par(mfrow = c(3,3))
for (i in 1:nrow(cibersort_rnaseq_sig_inflamed)) {
	boxplot(list("Inflamed A" = as.numeric(cibersort_rnaseq_sig_inflamed[row.names(cibersort_rnaseq_sig_inflamed)[i], inflamed_a]), "Inflamed B" = as.numeric(cibersort_rnaseq_sig_inflamed[row.names(cibersort_rnaseq_sig_inflamed)[i], inflamed_b])), main = row.names(cibersort_rnaseq_sig_inflamed)[i], ylab = "CIBERSORT Score")
}


#wilcoxon test a given group vs the rest (combined inflamed)
cibersort_score_testing = data.frame(row.names = row.names(cibersort_rnaseq), Inflamed = rep(NA, nrow(cibersort_rnaseq)), Redox = rep(NA, nrow(cibersort_rnaseq)), Mixed = rep(NA, nrow(cibersort_rnaseq)), stringsAsFactors = FALSE)
for (i in 1:nrow(cibersort_rnaseq)) {
	for(j in 1:3) {
		tryCatch({
			current_subtype = subset(cibersort_meta, subtype_numeric == j)$GenomicsSample
			other_subtype = subset(cibersort_meta, subtype_numeric != j)$GenomicsSample
			cibersort_score_testing[i, j] = wilcox.test(as.numeric(cibersort_rnaseq[i, current_subtype]), as.numeric(cibersort_rnaseq[i, other_subtype]))$p.value
		}, error = function(e) {
			NA
		})
	}
}

cibersort_sig = cibersort_score_testing[rowSums(cibersort_score_testing < 0.05) > 0, ]
cibersort_sig_for_radar = row.names(cibersort_score_testing[rowSums(cibersort_score_testing < 0.1) > 0, ])
#cibersort_sig = cibersort_score_testing
cibersort_rnaseq_sig = subset(cibersort_rnaseq, row.names(cibersort_rnaseq) %in% row.names(cibersort_sig))
par(mfrow = c(3,3))

for (i in 1:nrow(cibersort_rnaseq_sig)) {
	boxplot(list("Inflamed" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], c(inflamed_a, inflamed_b)]), "Redox" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], redoxgroup]), "Double Neg" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], Mixedgroup])), main = row.names(cibersort_rnaseq_sig)[i], ylab = "CIBERSORT Score")
}

#paper plots
png(file = "output/main/cibersort_all_2018-06-26.png", width = 900, height = 800)
par(mfrow = c(2,3))
par(mar = c(5.1, 5.1, 4.1, 2.1))
cibersort_sig = cibersort_score_testing[rowSums(cibersort_score_testing < 0.05) > 0, ]
#cibersort_sig = cibersort_score_testing
cibersort_rnaseq_sig = subset(cibersort_rnaseq, row.names(cibersort_rnaseq) %in% row.names(cibersort_sig))
cibersort_rnaseq_plot_mains = c("Memory B-cells", "Plasma Cells", "Resting NK Cells", "Monocytes", "M2 Macrophages", "Neutrophils", "Memory B-cells", "Regulatory T-cells")
for (i in 1:nrow(cibersort_rnaseq_sig)) {
	#boxplot(list("Inflamed" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], c(inflamed_a, inflamed_b)]), "Redox" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], redoxgroup]), "Mixed" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], Mixedgroup])), main = row.names(cibersort_rnaseq_sig)[i], ylab = "CIBERSORT Score", cex.axis = 2, cex.lab = 2, col = c(darkcolors[1], darkcolors[2], darkcolors[3]))
	#no main title
	boxplot(list("Inflamed" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], c(inflamed_a, inflamed_b)]), "Redox" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], redoxgroup]), "Mixed" = as.numeric(cibersort_rnaseq_sig[row.names(cibersort_rnaseq_sig)[i], Mixedgroup])), ylab = "CIBERSORT Score", cex.axis = 2, cex.lab = 2, cex.main = 2, col = c(darkcolors[1], darkcolors[2], darkcolors[3]), main = cibersort_rnaseq_plot_mains[i])
}
dev.off()

png(file = "output/main/cibersort_inflamedab_2018-06-26.png", width = 900, height = 800)
par(mar = c(5.1, 5.1, 4.1, 2.1))
par(mfrow = c(2,3))
cibersort_sig_inflamed = row.names(cibersort_score_testing_inflamed)[rowSums(cibersort_score_testing_inflamed < 0.1) > 0]
cibersort_rnaseq_sig_inflamed = subset(cibersort_rnaseq, row.names(cibersort_rnaseq) %in% cibersort_sig_inflamed)
#cibersort_rnaseq_sig_inflamed = cibersort_rnaseq #all samples
cibersort_rnaseq_inflamed_plot_mains = c("Memory B-cells", "Regulatory T-cells", "Resting NK Cells", "Activated Mast Cells", "Neutrophils")
for (i in 1:nrow(cibersort_rnaseq_sig_inflamed)) {
	#boxplot(list("Inflamed A" = as.numeric(cibersort_rnaseq_sig_inflamed[row.names(cibersort_rnaseq_sig_inflamed)[i], inflamed_a]), "Inflamed B" = as.numeric(cibersort_rnaseq_sig_inflamed[row.names(cibersort_rnaseq_sig_inflamed)[i], inflamed_b])), main = row.names(cibersort_rnaseq_sig_inflamed)[i], ylab = "CIBERSORT Score", cex.axis = 2, cex.lab = 2, col = darkcolors[1])
	#no main title
	boxplot(list("Inflamed A" = as.numeric(cibersort_rnaseq_sig_inflamed[row.names(cibersort_rnaseq_sig_inflamed)[i], inflamed_a]), "Inflamed B" = as.numeric(cibersort_rnaseq_sig_inflamed[row.names(cibersort_rnaseq_sig_inflamed)[i], inflamed_b])), ylab = "CIBERSORT Score", cex.axis = 2, cex.lab = 2, cex.main = 2, col = darkcolors[1], main = cibersort_rnaseq_inflamed_plot_mains[i])
}
dev.off()



#number of wt in each group ----
inflamed_muts = subset(mutation_status_prelim, tumor_06s %in% subset(tmt_meta, subtype == "Inflamed")$TissueID)
redox_muts = subset(mutation_status_prelim, tumor_06s %in% subset(tmt_meta, subtype == "Redox")$TissueID)
Mixed_muts = subset(mutation_status_prelim, tumor_06s %in% subset(tmt_meta, subtype == "Mixed")$TissueID)
table(inflamed_muts[, 3:ncol(inflamed_muts)] == "wt")
table(redox_muts[, 3:ncol(redox_muts)] == "wt")
table(Mixed_muts[, 3:ncol(Mixed_muts)] == "wt")
5617 / nrow(inflamed_muts)
6593 / nrow(redox_muts)
1816 / nrow(Mixed_muts)

#DRIVE Venn ----
#needs RNA-protein correlations, CNV-protein correlations, DRIVE data, and differential expression R code run first ----

drive_venn = list("High Corr" = correlation_table_cnv_protein_positive_matched$Accession, "Redox" = metabolismhigher_de_genes, "DRIVE" = drive_results$Symbol)
venn.diagram(drive_venn, "output/main/drive_correlation.5_fdr.25_2018-12-14.png", cex = 3)


#Proteomics PCA for manuscript ----

proteomics_10pct_pca_df = protein_expression_10pct[, tmt_meta$TissueID]
proteomics_10pct_pca_df[is.na(proteomics_10pct_pca_df)] = 0
proteomics_10pct_pca = prcomp(t(proteomics_10pct_pca_df), center = TRUE, scale = TRUE)

ggbiplot(proteomics_10pct_pca, obs.scale = 1, var.scale = 1, 
				 ellipse = FALSE, circle = TRUE, var.axes = FALSE, groups = as.factor(tmt_meta$subtype), size = 4) +
	#scale_color_discrete(name = '') +
	geom_point(aes(color = as.factor(tmt_meta$subtype)), size = 5) +
	scale_color_manual(name = "Proteomic Subtype", values = c(darkcolors[1], darkcolors[3], darkcolors[2])) +
	theme(legend.direction = 'horizontal', legend.position = 'top')
ggsave("output/main/proteomics_pca_manuscript_2018-05-31.png")

#APC lollipop plot ----
library(trackViewer)
library(GenomicRanges)

mutation_apc = read.delim("output/main/apc_mutation_type_position_expanded_2018-06-12.txt", header = TRUE, stringsAsFactors = FALSE)

apc_domains = read.delim("output/main/apc_cbioportal_pfam_mutations_2018-06-12.txt", header = TRUE, stringsAsFactors = FALSE)
apc_domains$seqnames = "chr1"
apc_domains[21, ] = c("start", "start", 1, 1, "chr1")

APC_mutations <- GRanges("chr1", IRanges(mutation_apc$numeric_positions, width=1, names = mutation_apc$NM_000038))
APC_mutations$color = mutation_apc$color
APC_mutations$border = mutation_apc$border
gene_regions = makeGRangesFromDataFrame(apc_domains[, c("seqnames","start", "end")])
gene_regions = setNames(gene_regions, apc_domains$PFAM_short)
gene_regions$fill = as.factor(apc_domains$PFAM_short)

legend = c(darkcolors[2], darkcolors[3])
names(legend) = c("Redox", "Mixed")
lolliplot(APC_mutations, gene_regions, legend = legend, xaxis = c(1, 500, 1000, 1500, 2000, 2500, 2843))

#image-related for manuscript reviews ----
neutrophil_proteins = protein_expression_10pct$Symbol[str_detect(protein_expression_10pct$Symbol, "DEFA")]
neutrophil_proteins = c(neutrophil_proteins, "ELANE")
inflamed_a_tissueid = subset(tmt_meta, tmt_meta$consensusClassK5_protein == 1)$TissueID
inflamed_b_tissueid = subset(tmt_meta, tmt_meta$consensusClassK5_protein == 3)$TissueID

boxplot(list(InflamedA = as.numeric(protein_expression_10pct[protein_expression_10pct$Symbol == neutrophil_proteins[2], inflamed_a_tissueid]), InflamedB = as.numeric(protein_expression_10pct[protein_expression_10pct$Symbol == neutrophil_proteins[2], inflamed_b_tissueid])))
#inflamed_a_neutrophil_highest
inflamed_a_neutrophil_highest = sort(protein_expression_10pct[protein_expression_10pct$Symbol == neutrophil_proteins[2], inflamed_a_tissueid], decreasing = TRUE)[1:4]
#inflamed_b_neutrophil_lowest
inflamed_b_neutrophil_lowest = sort(protein_expression_10pct[protein_expression_10pct$Symbol == neutrophil_proteins[2], inflamed_b_tissueid], decreasing = TRUE)[(length(inflamed_b_tissueid) - 2):length(inflamed_b_tissueid)]

ffpe = read.delim("data/main/ffpe_slide_mapping_table_08-22-17.txt", header = TRUE, stringsAsFactors = FALSE)
tmt_meta = merge(tmt_meta, ffpe[, c("Corresponding.Tissue.Core.Aliquot.ID", "FFPE.H.E.Slide.Label")], by.x = "TissueID", by.y = "Corresponding.Tissue.Core.Aliquot.ID")
subset(tmt_meta, TissueID %in% names(inflamed_a_neutrophil_highest) | TissueID %in% names(inflamed_b_neutrophil_lowest))[, c("TissueID", "FFPE.H.E.Slide.Label", "consensusClassK5_protein") ]

#pickles for manuscript reviews ----
pickles_meta = read.delim("output/manuscript_revisions/pickles_avana_lung_cell_lines_squamous_annotated_2018-11-14.txt", stringsAsFactors = FALSE)
pickles_meta_sqlc = subset(pickles_meta, histology == "squamous")
pickles_TP63 = read.csv("output/manuscript_revisions/pickles_lung_TP63_2018-11-14.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
names(pickles_TP63)[names(pickles_TP63) == "BF"] = "TP63"
pickles_AKR1C3 = read.csv("output/manuscript_revisions/pickles_lung_AKR1C3_2018-11-14.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
names(pickles_AKR1C3)[names(pickles_AKR1C3) == "BF"] = "AKR1C3"
pickles_PSAT1 = read.csv("output/manuscript_revisions/pickles_lung_PSAT1_2018-11-14.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
names(pickles_PSAT1)[names(pickles_PSAT1) == "BF"] = "PSAT1"
pickles_TFRC = read.csv("output/manuscript_revisions/pickles_lung_TFRC_2018-11-14.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
names(pickles_TFRC)[names(pickles_TFRC) == "BF"] = "TFRC"
pickles_SERPINB5 = read.csv("output/manuscript_revisions/pickles_lung_SERPINB5_2018-11-14.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
names(pickles_SERPINB5)[names(pickles_SERPINB5) == "BF"] = "SERPINB5"
pickles_GSR = read.csv("output/manuscript_revisions/pickles_lung_GSR_2018-11-14.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
names(pickles_GSR)[names(pickles_GSR) == "BF"] = "GSR"
pickles_bf_df = data.frame(row.names = pickles_meta$Cell_Line, TP63 = pickles_TP63[, "TP63"], AKR1C3 = pickles_AKR1C3[, "AKR1C3"], PSAT1 = pickles_PSAT1[, "PSAT1"], TFRC = pickles_TFRC[, "TFRC"], SERPINB5 = pickles_SERPINB5[, "SERPINB5"], stringsAsFactors = FALSE)
pickles_bf_df_sqlc = pickles_bf_df[row.names(pickles_bf_df) %in% pickles_meta_sqlc$Cell_Line, ]
pickles_bf_df_sqlc = pickles_bf_df_sqlc[pickles_meta_sqlc$Cell_Line, ]
pickles_bf_df_sqlc$density = pickles_meta_sqlc$nrf2_gina_paper == "yes" | pickles_meta_sqlc$meyerson_nrf2_keap1_mutation == "yes"
for (i in 1:nrow(pickles_bf_df_sqlc)) {
	if (pickles_bf_df_sqlc$density[i] == TRUE) {
		pickles_bf_df_sqlc$density[i] = 25
	} else {
		pickles_bf_df_sqlc$density[i] = 999
	}
}  

png(filename = "output/main/pickles_scc_drive_matched_2018-12-21.png", width = 1200, height = 1000)
par(mar = c(5.1, 6.1, 12, 2.1))
par(mfrow = c(2,3))

pickles_bf_df_sqlc = pickles_bf_df_sqlc[order(pickles_bf_df_sqlc$TP63, decreasing = T), ]
TP63_pickles_barplot = barplot(pickles_bf_df_sqlc$TP63, xaxt = "n", ylim = c(-50, 50), col = "#FF000A", ylab = "BF", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = pickles_bf_df_sqlc$density)
axis(3, at = TP63_pickles_barplot, labels = row.names(pickles_bf_df_sqlc), las = 2, cex.axis = 2)
abline(h = 0, lty = 1)
abline(h = 3, lty = 2)

pickles_bf_df_sqlc = pickles_bf_df_sqlc[order(pickles_bf_df_sqlc$PSAT1, decreasing = T), ]
PSAT1_pickles_barplot = barplot(pickles_bf_df_sqlc$PSAT1, xaxt = "n", ylim = c(-50, 50), col = "#FF000A", ylab = "BF", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = pickles_bf_df_sqlc$density)
axis(3, at = PSAT1_pickles_barplot, labels = row.names(pickles_bf_df_sqlc), las = 2, cex.axis = 2)
abline(h = 0, lty = 1)
abline(h = 3, lty = 2)


pickles_bf_df_sqlc = pickles_bf_df_sqlc[order(pickles_bf_df_sqlc$AKR1C3, decreasing = T), ]
AKR1C3_pickles_barplot = barplot(pickles_bf_df_sqlc$AKR1C3, xaxt = "n", ylim = c(-50, 50), col = "#FF000A", ylab = "BF", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = pickles_bf_df_sqlc$density)
axis(3, at = AKR1C3_pickles_barplot, labels = row.names(pickles_bf_df_sqlc), las = 2, cex.axis = 2)
abline(h = 0, lty = 1)
abline(h = 3, lty = 2)


pickles_bf_df_sqlc = pickles_bf_df_sqlc[order(pickles_bf_df_sqlc$TFRC, decreasing = T), ]
TFRC_pickles_barplot = barplot(pickles_bf_df_sqlc$TFRC, xaxt = "n", ylim = c(-50, 50), col = "#CB6527", ylab = "BF", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = pickles_bf_df_sqlc$density)
axis(3, at = TFRC_pickles_barplot, labels = row.names(pickles_bf_df_sqlc), las = 2, cex.axis = 2)
abline(h = 0, lty = 1)
abline(h = 3, lty = 2)

pickles_bf_df_sqlc = pickles_bf_df_sqlc[order(pickles_bf_df_sqlc$SERPINB5, decreasing = T), ]
SERPINB5_pickles_barplot = barplot(pickles_bf_df_sqlc$SERPINB5, xaxt = "n", ylim = c(-50, 50), col = "#CB6527", ylab = "BF", main = "", cex.lab = 2, cex.axis = 2, angle = 45, density = pickles_bf_df_sqlc$density)
axis(3, at = SERPINB5_pickles_barplot, labels = row.names(pickles_bf_df_sqlc), las = 2, cex.axis = 2)
abline(h = 0, lty = 1)
abline(h = 3, lty = 2)
dev.off()


#revisions: transcript-protein correlations in redox and inflamed ----
#need section correlation: protein/rnaseq run
inflamed_proteomics_names = subset(tmt_meta, subtype == "Inflamed")$TissueID
redox_proteomics_names = subset(tmt_meta, subtype == "Redox")$TissueID
inflamed_rnaseq_names = subset(tmt_meta, subtype == "Inflamed")$GenomicsSample
redox_rnaseq_names = subset(tmt_meta, subtype == "Redox")$GenomicsSample
mixed_proteomics_names = subset(tmt_meta, subtype == "Mixed")$TissueID
mixed_rnaseq_names = subset(tmt_meta, subtype == "Mixed")$GenomicsSample


correlation_table_rnaseq_by_gene_inflamed = data.frame(Accession = protein_rnaseq_expression_for_correlations$Symbol, spearman_corr = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(protein_rnaseq_expression_for_correlations)) {
	corr_result = tryCatch({
		proteins_for_corr = as.numeric(protein_rnaseq_expression_for_correlations[i, inflamed_proteomics_names])
		genes_for_corr = as.numeric(protein_rnaseq_expression_for_correlations[i, inflamed_rnaseq_names])
		if (sum(!is.na(proteins_for_corr)) >= number_of_samples * 0.1 & sum(!is.na(genes_for_corr)) >= number_of_samples * 0.1) {
			cor.test(proteins_for_corr, genes_for_corr, method = "spearman", use = "pairwise.complete.obs")
		} else {
			NA
		}
	}, error = function(e) {
		NA
	})
	if (is.na(corr_result)) {
		correlation_table_rnaseq_by_gene_inflamed$spearman_corr[i] = NA
		correlation_table_rnaseq_by_gene_inflamed$pval[i] = NA
	} else {
		correlation_table_rnaseq_by_gene_inflamed$spearman_corr[i] = corr_result$estimate
		correlation_table_rnaseq_by_gene_inflamed$pval[i] = corr_result$p.value
	}
}
correlation_table_rnaseq_by_gene_inflamed$padj = p.adjust(correlation_table_rnaseq_by_gene_inflamed$pval, method = "BH")


correlation_table_rnaseq_by_gene_redox = data.frame(Accession = protein_rnaseq_expression_for_correlations$Symbol, spearman_corr = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(protein_rnaseq_expression_for_correlations)) {
	corr_result = tryCatch({
		proteins_for_corr = as.numeric(protein_rnaseq_expression_for_correlations[i, redox_proteomics_names])
		genes_for_corr = as.numeric(protein_rnaseq_expression_for_correlations[i, redox_rnaseq_names])
		if (sum(!is.na(proteins_for_corr)) >= number_of_samples * 0.1 & sum(!is.na(genes_for_corr)) >= number_of_samples * 0.1) {
			cor.test(proteins_for_corr, genes_for_corr, method = "spearman", use = "pairwise.complete.obs")
		} else {
			NA
		}
	}, error = function(e) {
		NA
	})
	if (is.na(corr_result)) {
		correlation_table_rnaseq_by_gene_redox$spearman_corr[i] = NA
		correlation_table_rnaseq_by_gene_redox$pval[i] = NA
	} else {
		correlation_table_rnaseq_by_gene_redox$spearman_corr[i] = corr_result$estimate
		correlation_table_rnaseq_by_gene_redox$pval[i] = corr_result$p.value
	}
}
correlation_table_rnaseq_by_gene_redox$padj = p.adjust(correlation_table_rnaseq_by_gene_redox$pval, method = "BH")



correlation_table_rnaseq_by_gene_mixed = data.frame(Accession = protein_rnaseq_expression_for_correlations$Symbol, spearman_corr = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(protein_rnaseq_expression_for_correlations)) {
	corr_result = tryCatch({
		proteins_for_corr = as.numeric(protein_rnaseq_expression_for_correlations[i, mixed_proteomics_names])
		genes_for_corr = as.numeric(protein_rnaseq_expression_for_correlations[i, mixed_rnaseq_names])
		if (sum(!is.na(proteins_for_corr)) >= number_of_samples * 0.1 & sum(!is.na(genes_for_corr)) >= number_of_samples * 0.1) {
			cor.test(proteins_for_corr, genes_for_corr, method = "spearman", use = "pairwise.complete.obs")
		} else {
			NA
		}
	}, error = function(e) {
		NA
	})
	if (is.na(corr_result)) {
		correlation_table_rnaseq_by_gene_mixed$spearman_corr[i] = NA
		correlation_table_rnaseq_by_gene_mixed$pval[i] = NA
	} else {
		correlation_table_rnaseq_by_gene_mixed$spearman_corr[i] = corr_result$estimate
		correlation_table_rnaseq_by_gene_mixed$pval[i] = corr_result$p.value
	}
}
correlation_table_rnaseq_by_gene_mixed$padj = p.adjust(correlation_table_rnaseq_by_gene_mixed$pval, method = "BH")


#analyses suggested from 2018-12-17 meeting
SOX2_targets = read.delim("data/manuscript_revisions/download/BENPORATH_SOX2_TARGETS_msigdb_2018-12-03.txt", header = FALSE, stringsAsFactors = FALSE)
NFE2L2_targets = read.delim("data/manuscript_revisions/download/nrf2_01_2018-12-20.txt", header = FALSE, stringsAsFactors = FALSE)
TP63_targets = read.delim("data/manuscript_revisions/download/perez_TP63_targets_2018-12-19.txt", header = FALSE, stringsAsFactors = FALSE)
SOX2_targets$SOX2_target = "TRUE"
names(SOX2_targets)[names(SOX2_targets) == "V1"] = "Accession"
NFE2L2_targets$NFE2L2_target = "TRUE"
names(NFE2L2_targets)[names(NFE2L2_targets) == "V1"] = "Accession"
TP63_targets$TP63_target = "TRUE"
names(TP63_targets)[names(TP63_targets) == "V1"] = "Accession"

SOX2_list = list(Overall = subset(correlation_table_rnaseq_by_gene, Accession %in% SOX2_targets$Accession)$spearman_corr, Inflamed = subset(correlation_table_rnaseq_by_gene_inflamed, Accession %in% SOX2_targets$Accession)$spearman_corr, Redox = subset(correlation_table_rnaseq_by_gene_redox, Accession %in% SOX2_targets$Accession)$spearman_corr, Mixed = subset(correlation_table_rnaseq_by_gene_mixed, Accession %in% SOX2_targets$Accession)$spearman_corr)
lapply(SOX2_list, function(x) {mean(na.omit(x))})
wilcox.test(SOX2_list$Inflamed, SOX2_list$Redox)
NFE2L2_list = list(Overall = subset(correlation_table_rnaseq_by_gene, Accession %in% NFE2L2_targets$Accession)$spearman_corr, Inflamed = subset(correlation_table_rnaseq_by_gene_inflamed, Accession %in% NFE2L2_targets$Accession)$spearman_corr, Redox = subset(correlation_table_rnaseq_by_gene_redox, Accession %in% NFE2L2_targets$Accession)$spearman_corr, Mixed = subset(correlation_table_rnaseq_by_gene_mixed, Accession %in% NFE2L2_targets$Accession)$spearman_corr)
lapply(NFE2L2_list, function(x) {mean(na.omit(x))})
wilcox.test(NFE2L2_list$Inflamed, NFE2L2_list$Redox)
TP63_list = list(Overall = subset(correlation_table_rnaseq_by_gene, Accession %in% TP63_targets$Accession)$spearman_corr, Inflamed = subset(correlation_table_rnaseq_by_gene_inflamed, Accession %in% TP63_targets$Accession)$spearman_corr, Redox = subset(correlation_table_rnaseq_by_gene_redox, Accession %in% TP63_targets$Accession)$spearman_corr, Mixed = subset(correlation_table_rnaseq_by_gene_mixed, Accession %in% TP63_targets$Accession)$spearman_corr)
lapply(TP63_list, function(x) {mean(na.omit(x))})
wilcox.test(TP63_list$Inflamed, TP63_list$Redox)


all(correlation_table_rnaseq_by_gene_inflamed$Accession == correlation_table_rnaseq_by_gene_redox$Accession)
all(correlation_table_rnaseq_by_gene_redox$Accession == correlation_table_rnaseq_by_gene$Accession)
head(correlation_table_rnaseq_by_gene_inflamed)

names(correlation_table_rnaseq_by_gene_inflamed) = c("Accession", "spearman_corr_inflamed", "pval_inflamed", "padj_inflamed")
names(correlation_table_rnaseq_by_gene_redox) = c("Accession", "spearman_corr_redox", "pval_redox", "padj_redox")

combined_correlation_table = cbind(correlation_table_rnaseq_by_gene[, c("Accession", "spearman_corr", "pval", "padj")], correlation_table_rnaseq_by_gene_inflamed[, -1], correlation_table_rnaseq_by_gene_redox[, -1])
combined_correlation_table = merge(combined_correlation_table, SOX2_targets, by = "Accession", all = TRUE)
combined_correlation_table = merge(combined_correlation_table, NFE2L2_targets, by = "Accession", all = TRUE)
combined_correlation_table = merge(combined_correlation_table, TP63_targets, by = "Accession", all = TRUE)
combined_correlation_table = subset(combined_correlation_table, SOX2_target == TRUE | NFE2L2_target == TRUE | TP63_target == TRUE)
names(combined_correlation_table)[1:4] = c("Gene Symbol", "spearman_corr_all", "pval_all", "padj_all")
combined_correlation_table = subset(combined_correlation_table, !is.na(spearman_corr_all) | !is.na(spearman_corr_inflamed) | !is.na(spearman_corr_redox))


write.table(combined_correlation_table, "output/main/redox_inflamed_correlation_tf_targets_2019-01-07.txt", sep = "\t", quote = F, row.names = F)

#keap1 mutant vs keap1 wt protein-protein correlations ----
keap1_mut_samples = subset(tmt_meta, KEAP1 != "wt")$TissueID
keap1_wt_samples = subset(tmt_meta, KEAP1 == "wt")$TissueID

keap1_mut_df = protein_expression_10pct[, c("Symbol", keap1_mut_samples)]
keap1_wt_df = protein_expression_10pct[, c("Symbol", keap1_wt_samples)]


keap1_mut_correlations_df = data.frame(Accession = protein_expression_10pct$Symbol, spearman_corr = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(keap1_mut_df)) {
	corr_result = tryCatch({
		keap1_expression = as.numeric(keap1_mut_df[keap1_mut_df$Symbol == "KEAP1", keap1_mut_samples])
		keap1_protein_expression = as.numeric(keap1_mut_df[i, keap1_mut_samples])
		if (sum(!is.na(keap1_expression)) >= number_of_samples * 0.1 & sum(!is.na(keap1_protein_expression)) >= number_of_samples * 0.1) {
			cor.test(keap1_expression, keap1_protein_expression, method = "spearman", use = "pairwise.complete.obs")
		} else {
			NA
		}
	}, error = function(e) {
		NA
	})
	if (is.na(corr_result)) {
		keap1_mut_correlations_df$spearman_corr[i] = NA
		keap1_mut_correlations_df$pval[i] = NA
	} else {
		keap1_mut_correlations_df$spearman_corr[i] = corr_result$estimate
		keap1_mut_correlations_df$pval[i] = corr_result$p.value
	}
}
keap1_mut_correlations_df$padj = p.adjust(keap1_mut_correlations_df$pval, method = "BH")



keap1_wt_correlations_df = data.frame(Accession = protein_expression_10pct$Symbol, spearman_corr = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(keap1_wt_df)) {
	corr_result = tryCatch({
		keap1_expression = as.numeric(keap1_wt_df[keap1_wt_df$Symbol == "KEAP1", keap1_wt_samples])
		keap1_protein_expression = as.numeric(keap1_wt_df[i, keap1_wt_samples])
		if (sum(!is.na(keap1_expression)) >= number_of_samples * 0.1 & sum(!is.na(keap1_protein_expression)) >= number_of_samples * 0.1) {
			cor.test(keap1_expression, keap1_protein_expression, method = "spearman", use = "pairwise.complete.obs")
		} else {
			NA
		}
	}, error = function(e) {
		NA
	})
	if (is.na(corr_result)) {
		keap1_wt_correlations_df$spearman_corr[i] = NA
		keap1_wt_correlations_df$pval[i] = NA
	} else {
		keap1_wt_correlations_df$spearman_corr[i] = corr_result$estimate
		keap1_wt_correlations_df$pval[i] = corr_result$p.value
	}
}
keap1_wt_correlations_df$padj = p.adjust(keap1_wt_correlations_df$pval, method = "BH")

keap1_mut_correlations_df$mut_wt_difference = keap1_mut_correlations_df$spearman_corr - keap1_wt_correlations_df$spearman_corr
keap1_mut_correlations_df_0.5 = subset(keap1_mut_correlations_df, mut_wt_difference >= 0.5)

subset(keap1_mut_correlations_df, padj <= 0.25 & Accession %in% keap1_mut_correlations_df_0.5)
write.table(keap1_mut_correlations_df_0.5, "output/main/keap1_correlations_2018-20-18.txt", sep = "\t", quote = F, row.names = F)

# revisions: heatmap protein correlations ----
protein1kheatmap_rnaseq_expression_for_correlations = merge(protein_expression_heatmap, rnaseq_expression, by.x = "Symbol", by.y = "Symbol")
correlation_table_rnaseq_by_gene = data.frame(Accession = protein1kheatmap_rnaseq_expression_for_correlations$Symbol, spearman_corr = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(protein1kheatmap_rnaseq_expression_for_correlations)) {
	corr_result = tryCatch({
		proteins_for_corr = as.numeric(protein1kheatmap_rnaseq_expression_for_correlations[i, tmt_meta$TissueID])
		genes_for_corr = as.numeric(protein1kheatmap_rnaseq_expression_for_correlations[i, tmt_meta$GenomicsSample])
		if (sum(!is.na(proteins_for_corr)) >= number_of_samples * 0.1 & sum(!is.na(genes_for_corr)) >= number_of_samples * 0.1) {
			cor.test(proteins_for_corr, genes_for_corr, method = "spearman", use = "pairwise.complete.obs")
		} else {
			NA
		}
	}, error = function(e) {
		NA
	})
	if (is.na(corr_result)) {
		correlation_table_rnaseq_by_gene$spearman_corr[i] = NA
		correlation_table_rnaseq_by_gene$pval[i] = NA
	} else {
		correlation_table_rnaseq_by_gene$spearman_corr[i] = corr_result$estimate
		correlation_table_rnaseq_by_gene$pval[i] = corr_result$p.value
	}
}
correlation_table_rnaseq_by_gene$padj = p.adjust(correlation_table_rnaseq_by_gene$pval, method = "BH")
dim(correlation_table_rnaseq_by_gene)
mean(na.omit(correlation_table_rnaseq_by_gene$spearman_corr))

#4 bicarbonate transport (GO:0015701) padj 0.0003474
#3 glutathione metabolic process (GO:0006749) padj 4.053e-10
#2 platelet degranulation (GO:0002576) padj 1.453e-19
#5 extracellular matrix organization (GO:0030198) padj 9.660e-22
#1 neutrophil degranulation (GO:0043312) padj 4.944e-21


mean(na.omit(subset(correlation_table_rnaseq_by_gene, correlation_table_rnaseq_by_gene$Accession %in% pathway_list$`neutrophil degranulation (GO:0043312)`)$spearman_corr))
mean(na.omit(subset(correlation_table_rnaseq_by_gene, correlation_table_rnaseq_by_gene$Accession %in% pathway_list$`extracellular matrix organization (GO:0030198)`)$spearman_corr))
mean(na.omit(subset(correlation_table_rnaseq_by_gene, correlation_table_rnaseq_by_gene$Accession %in% pathway_list$`platelet degranulation (GO:0002576)`)$spearman_corr))
mean(na.omit(subset(correlation_table_rnaseq_by_gene, correlation_table_rnaseq_by_gene$Accession %in% pathway_list$`glutathione metabolic process (GO:0006749)`)$spearman_corr))
mean(na.omit(subset(correlation_table_rnaseq_by_gene, correlation_table_rnaseq_by_gene$Accession %in% pathway_list$`bicarbonate transport (GO:0015701)`)$spearman_corr))

#tumor cellularity ----
cellularity_list = list(Inflamed = subset(tmt_meta, subtype == "Inflamed")$Tumor.Cellularity.Percentage...PQC, Redox = subset(tmt_meta, subtype == "Redox")$Tumor.Cellularity.Percentage...PQC, Mixed = subset(tmt_meta, subtype == "Mixed")$Tumor.Cellularity.Percentage...PQC)
												
wilcox.test(subset(tmt_meta, subtype == "Mixed")$Tumor.Cellularity.Percentage...PQC, subset(tmt_meta, subtype != "Mixed")$Tumor.Cellularity.Percentage...PQC)

#how many cnv in each group? ----
cnv_long = read.delim("data/main/Standard_400kbp.50counts_updated_processed_20171220.txt", header = TRUE, stringsAsFactors = FALSE)

#overall
length(unique(cnv_long_cnstate_gain$chromosome_band))
length(unique(cnv_long_cnstate_loss$chromosome_band))

redox_patients = subset(tmt_meta, subtype == "Redox")$TissueID
inflamed_patients = subset(tmt_meta, subtype == "Inflamed")$TissueID
mixed_patients = subset(tmt_meta, subtype == "Mixed")$TissueID

length(unique(subset(cnv_long_cnstate_gain, TissueID %in% redox_patients)$chromosome_band))
length(unique(subset(cnv_long_cnstate_loss, TissueID %in% redox_patients)$chromosome_band))

length(unique(subset(cnv_long_cnstate_gain, TissueID %in% inflamed_patients)$chromosome_band))
length(unique(subset(cnv_long_cnstate_loss, TissueID %in% inflamed_patients)$chromosome_band))

length(unique(subset(cnv_long_cnstate_gain, TissueID %in% mixed_patients)$chromosome_band))
length(unique(subset(cnv_long_cnstate_loss, TissueID %in% mixed_patients)$chromosome_band))


dim(subset(cnv_long, TissueID %in% subset(tmt_meta, subtype == "Inflamed")$TissueID & Type == "Gain"))
dim(subset(cnv_long, TissueID %in% subset(tmt_meta, subtype == "Redox")$TissueID))
dim(subset(cnv_long, TissueID %in% subset(tmt_meta, subtype == "Mixed")$TissueID))


#cd33 from james ----
cd33_table = read.delim("data/manuscript_revisions/james_cd33_ihc_scoring_2019-01-11.txt", header = TRUE, stringsAsFactors = FALSE)
tmt_meta_for_cd33 = subset(tmt_meta, consensusClassK5_protein == 1 | consensusClassK5_protein == 3)
cd_33_scoring = merge(tmt_meta_for_cd33[, c("consensusClassK5_protein", "Deidentified.Patient.ID.x")], cd33_table, by.x = "Deidentified.Patient.ID.x", by.y = "De.identified.Patient.ID")

#png("output/main/spine_mosaic_TOL_plot_2018-06-11.png", width = 400, height = 600)
yellow_to_red = brewer.pal(9, "YlOrRd")
cd33_factor = factor(cd_33_scoring$S.CD33)
subtype_factor = factor(cd_33_scoring$consensusClassK5_protein)
spineplot(cd33_factor ~ subtype_factor, col = yellow_to_red[c(5, 7, 9)], ylab = "CD33-S Score", xlab = "", off = 1, xaxlabels	
 = c("InflamedA", "InflamedB"))
legend("bottomleft", legend = c("0", "1", "2", "3"), fill = c(yellow_to_red[3], yellow_to_red[5], yellow_to_red[7], yellow_to_red[9]), cex = 1)

spineplot(cd33_factor ~ subtype_factor, col = yellow_to_red[c(3, 5, 7, 9)], ylab = "CD33-S Score", xlab = "", off = 1, xaxlabels	
					= c("InflamedA", "InflamedB"))
legend("bottomleft", legend = c("0", "1", "2", "3"), fill = c(yellow_to_red[3], yellow_to_red[5], yellow_to_red[7], yellow_to_red[9]), cex = 1)

#CIBERSORT xCell comparison ----
xcell_cibersort_mapping = read.delim("data/manuscript_revisions/cibersort_xcell_name_mapping_2019-04-08.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
xcell_cibersort_mapping = subset(xcell_cibersort_mapping, xCell != "")
xcell_cibersort_mapping$spearman_estimate = NA
cibersort_rnaseq = read.delim("data/main/original_data/CIBERSORT_output_rna_for_cibersort_2018-05-15.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

cibersort_rnaseq = as.data.frame(t(cibersort_rnaseq))

xcell_rnaseq = read.delim("output/main/xCell_rnaseq_for_cibersort_2018-05-14_xCell_1434092518.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

for (i in 1:nrow(xcell_cibersort_mapping)) {
	xcell_cibersort_mapping$spearman_estimate[i] = cor.test(as.numeric(xcell_rnaseq[xcell_cibersort_mapping[i, "xCell"], tmt_meta$GenomicsSample]), as.numeric(cibersort_rnaseq[xcell_cibersort_mapping[i, "CIBERSORT"], tmt_meta$GenomicsSample]), method = "spearman")$estimate
	
}

write.table(xcell_cibersort_mapping, "output/manuscript_revisions/xcell_cibersort_correlations_2019-04-08.txt", sep = "\t", quote = F)

#tSNE ----
set.seed(1234)
tmt_meta$tsne_coloring = tmt_meta$subtype
tmt_meta$tsne_coloring = str_replace(tmt_meta$tsne_coloring, "Inflamed", darkcolors[1])
tmt_meta$tsne_coloring = str_replace(tmt_meta$tsne_coloring, "Redox", darkcolors[2])
tmt_meta$tsne_coloring = str_replace(tmt_meta$tsne_coloring, "Mixed", darkcolors[3])

head(protein_expression_10pct[, tmt_meta$TissueID])
d = dist(t(protein_expression_10pct[, tmt_meta$TissueID]))
tsne_out = Rtsne(d, is_distance=TRUE, perplexity=10, verbose = TRUE)
plot(tsne_out$Y, col = tmt_meta$tsne_coloring, pch = 16, main='Proteogenomics protein tSNE (10pct)', cex = 2, yaxt = "n", xaxt = "n", bty = "n", xlab = "", ylab = "")
legend("topleft", legend = c("Inflamed", "Redox", "Mixed"), col = c(darkcolors[1], darkcolors[2], darkcolors[3]), pch = 16, bty = "n")


proteomics_10pct_pca_df = protein_expression_10pct[, tmt_meta$TissueID]
proteomics_10pct_pca_df[is.na(proteomics_10pct_pca_df)] = 0
proteomics_10pct_pca = prcomp(t(proteomics_10pct_pca_df), center = TRUE, scale = TRUE)
summary(proteomics_10pct_pca)
proteomics_10pct_pca_50pc = proteomics_10pct_pca$x[, 1:50]

par(mfrow = c(2,3))
for (i in c(5, 10, 15, 25, 30, 35)) {
	tsne_pca_out = Rtsne(proteomics_10pct_pca_50pc, is_distance = FALSE, perplexity = i, verbose = TRUE, pca = FALSE)
	plot(tsne_pca_out$Y, col = tmt_meta$tsne_coloring, pch = 16, main = paste("tSNE first 50 pcs, perplexity = ", i), cex = 2, yaxt = "n", xaxt = "n", bty = "n", xlab = "", ylab = "")
	
}
