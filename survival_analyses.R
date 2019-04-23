
library(readxl)
library(plyr)
library(tidyverse)
library(broom)
library(compareGroups)
library(Hmisc)
library(lubridate)
library(janitor)
library(survival)
library(KMsurv)
library(survminer)
library(rio)
library(qvalue)
library(parallel)
library(ggpubr)
library(vcdExtra)
library(reshape2)
library(meta)
library(metafor)
library(glue)

# Function to calculate row summary
rowSummary <- function(x, ...) apply(X = x, MARGIN = 1, ...)

# Clinical data 
pgx_followup <- read_excel("PGX follow-up data final_jewel_2.5.18.xlsx", 
                           sheet = "Data", skip = 1) %>% 
  clean_names() 

# Data for smoking status
smoking_status <- read_excel("PGX follow-up data final_jewel_2.5.18.xlsx", 
                             sheet = "Smoking Status", skip = 1) %>% 
  clean_names() 

# Data for clusters/subtypes info
protein_subtypes <- read_delim(
  "proteogenomics_sample_class_immune_subtypes_2018-02-06.txt",
  delim = "\t") %>% clean_names() 

# Data for pathology features
pathology <- read_delim("SQLC116-HE_scores_20180111.txt", delim = "\t") %>% 
  clean_names() 

# Merge all four datasets above
clinical_data <- list(pgx_followup, smoking_status, 
                      protein_subtypes, pathology) %>% 
  reduce2(list(c("deidentified_patient_id", "aliquot_sample_id" = "sample_id"), 
               c("aliquot_sample_id" = "tissue_id"), 
               "deidentified_patient_id"), 
          left_join) 

# Format stage and lymph node status
# regional_lymph_nodes_examined value 97 means unknown/not stated
clinical_data <- clinical_data %>% 
  mutate(
    gender_cancer_registry = str_to_title(gender_cancer_registry),
    ethnicity_cancer_registry = recode_factor(
      ethnicity_cancer_registry, "NON-SPANISH" = "Non-hispanic",
      "PUERTO RICAN" = "Hispanic", "UNKNOWN" = "Unknown"),
    race_cancer_registry_1 = recode_factor(
      race_cancer_registry_1, "WHITE" = "White", "BLACK" = "Black"),
    ajcc_7_stage_grp2 = str_replace_all(ajcc_7_stage_grp, "A|B", ""),
    ajcc_7_stage_grp3 = recode_factor(
      ajcc_7_stage_grp2, "I" = "I", "II" = "II/III","III" = "II/III"),
    grade_differentiation = recode_factor(
      grade_differentiation, "POORLY DIFFEREN." = "Poor", 
      "MODERATELY DIFFEREN." = "Moderate",
      "WELL DIFFERENTIATED" = "Well", "NOT DETERMINED OR NA" = "N/A"),
    grade_differentiation2 = recode_factor(
      grade_differentiation, "Poor" = "Poor", "Moderate" = "Moderate/Well",
      "Well" = "Moderate/Well", "N/A" = NA_character_),
    recently_reported_cigarette_smoking_status = recode_factor(
      recently_reported_cigarette_smoking_status, 
      "NEVER" = "Never-smoker", "EVER" = "Ever/missing", 
      "MISSING" = "Ever/missing"),
    lymph_nodes_status = if_else(
      regional_lymph_nodes_positive == 0, "Negative", "Positive"),
    protein_group2 = case_when(
      protein_group == "Inflamed" & subtype == "A" ~ "Inflamed A",
      protein_group == "Inflamed" & subtype == "B" ~ "Inflamed B", 
      protein_group == "Redox" ~ "Redox", 
      protein_group == "DoubleNeg" ~ "Double negative"),
    protein_group2 = factor(
      protein_group2, levels = c("Inflamed A", "Inflamed B", "Redox", 
                                 "Double negative")),
    protein_group = recode_factor(
      protein_group, "Inflamed" = "Inflamed", "Redox" = "Redox",
      "DoubleNeg" = "Double negative"),
    # tln2 = if_else(tln %in% 2:3, "2/3", as.character(tln)),
    tln2 = if_else(tln == 0, "0", "Others"),
    subtype = recode_factor(subtype, "A" = "Inflamed A", "B" = "Inflamed B"))

# Summary variables
clinical_vars <- c(
  "age_at_diagnosis", "gender_cancer_registry", "ethnicity_cancer_registry", 
  "race_cancer_registry_1", "ajcc_7_stage_grp2", "grade_differentiation", 
  "grade_differentiation2", "recently_reported_cigarette_smoking_status", 
  "lymph_nodes_status", "tln", "tln2")

# Clinical variable labels
clinical_var_labels <- c(
  "Age at Diagnosis", "Gender", "Ethnicity", "Race", "AJCC-7 Stage", 
  "Grade/Differentiation", "Grade/Differentiation", "Smoking", "Lymph nodes", 
  "TLN", "TLN")

# Label for summary variables
clinical_data <- Hmisc::upData(clinical_data, 
        labels = set_names(clinical_var_labels, clinical_vars))

# Summary table using compareGroups package
summary_clinical_overall <- compareGroups(
  ~ age_at_diagnosis + gender_cancer_registry + ethnicity_cancer_registry + 
    race_cancer_registry_1 + ajcc_7_stage_grp2 + grade_differentiation + 
    recently_reported_cigarette_smoking_status + lymph_nodes_status + tln + tln2, 
  data = as.data.frame(clinical_data), method = c(tln = 3)) %>% 
  createTable(digits = 1)


# Summary by three protein groups -----------------------------------------
summary_three_protein <- compareGroups(
  protein_group ~ age_at_diagnosis + gender_cancer_registry + ethnicity_cancer_registry + 
    race_cancer_registry_1 + ajcc_7_stage_grp2 + grade_differentiation + grade_differentiation2 +
    recently_reported_cigarette_smoking_status + lymph_nodes_status + tln + tln2, 
  data = as.data.frame(clinical_data), method = c(tln = 3)) %>% 
  createTable(digits = 1, show.all = TRUE)

# Mantel-Haenszel test for linear trend in tln by three protein groups 
mh_tln <- CMHtest(
  ~ protein_group + tln, data = clinical_data, types = "cmeans")$table 

# Summary by four protein groups ------------------------------------------
summary_four_protein <- compareGroups(
  protein_group2 ~ age_at_diagnosis + gender_cancer_registry + ethnicity_cancer_registry + 
    race_cancer_registry_1 + ajcc_7_stage_grp2 + grade_differentiation + grade_differentiation2 +
    recently_reported_cigarette_smoking_status + lymph_nodes_status + tln + tln2, 
  data = as.data.frame(clinical_data), method = c(tln = 3)) %>% 
  createTable(digits = 1, show.all = TRUE)

# Mantel-Haenszel test for linear trend in tln by four protein groups
mh_tln2 <- CMHtest(
  ~ protein_group2 + tln, data = clinical_data, types = "cmeans")$table 


# Subtype analysis --------------------------------------------------------
summary_subtype <- filter(clinical_data, !is.na(subtype)) %>% 
  # droplevels() %>% 
  as.data.frame() %>% 
  compareGroups(
    subtype ~ ajcc_7_stage_grp2 + grade_differentiation + 
      grade_differentiation2 + tln + tln2, data = ., method = c(tln = 3)) %>% 
  createTable(digits = 1, show.all = TRUE)

# Mantel-Haenszel test for linear trend in tln by subtypes
mh_tln_subtype <- CMHtest(
  ~ subtype + tln, data = filter(clinical_data, !is.na(subtype)), 
  types = "cmeans")$table 


## ---- survival
# OS and RFS ----------------------------------------------------------------------------
# Variables related to OS and RFS
clinical_data <- clinical_data %>% 
  mutate(rfs_status = if_else(recurrence_final == "1" | vital_status_final == "Dead", 1, 0),
         os_status = if_else(vital_status_final == "Dead", 1, 0),
         overall_survival_2_years = overall_survival_2_months / 12,
         recurrence_free_survival_2_years = recurrence_free_survival_2_months_specimen_collection_date_to_date_of_first_clinically_confirmed_recurrence_date_of_last_contact_death_whichever_occurred_first / 12)

# ggplot theme to be used for all OS and RFS plots
mytheme <- theme(plot.title = element_text(hjust = 0.5),
                 axis.text = element_text(size = 12),
                 axis.title = element_text(size = 14),
                 # legend.position = c(0.8, 0.7),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 11),
                 legend.key = element_blank(),
                 legend.key.width = unit(1.3,"cm"),
                 legend.background = element_blank())

## KM plots of OS and RFS by different variables
# Strata variables for fitting OS and PRF
strata_vars <- c(
  "protein_group", "protein_group2", "ajcc_7_stage_grp2", 
  "gender_cancer_registry", "grade_differentiation2", "tln2", "subtype")

strata_var_labels <- c(
  "Protein group", "Protein group", "AJCC-7 Stage", "Gender", 
  "Grade/differentiation", "TLN", "Subtype")

# Survival time variables
surv_time_vars <- c("overall_survival_2_years", "recurrence_free_survival_2_years")

# Survival time variable labs
surv_time_var_labels <- c("OS from date of sample collection", "RFS from date of sample collection")

# Survival status variables
status_vars <- c("os_status", "rfs_status")

# Plot legend labs
legend_labs <- list(
  c("Inflamed", "Redox", "Double negative"), 
  c("Inflamed A", "Inflamed B","Redox", "Double negative"),
  c("Stage I", "Stage II", "Stage III"), c("Female", "Male"), c("Poor", "Moderate/Well"),
  c("No TLN", "With TLN"), c("Inflamed A", "Inflamed B"))

# Y-axis label for OS and PFS plots
ylab <- c("Overall survival", "Recurrence free survival")

# KM curves
km_os_rfs <- lapply(seq_along(strata_vars), function(x){
  lapply(seq_along(surv_time_vars), function(y){
    # Surv object
    surv_obj = as.formula(paste0("Surv(", surv_time_vars[y], ", ", status_vars[y], ") ~ ", 
                                 strata_vars[x]))
    # KM fit
    fit = survfit(surv_obj, data = clinical_data)
    # Set formula object to surv_obj prevent error " Error: object of type 'symbol' is not subsettable"
    fit$call$formula = surv_obj
    # KM plot
    kmPlot = ggsurvplot(fit, data = clinical_data, size = 0.5, pval = TRUE, ggtheme = theme_bw(),
                        linetype = seq_along(unique(clinical_data[[strata_vars[x]]])), 
                        palette = c("red", "blue", "black", "bisque4"), title = paste0(strata_var_labels[x], " (", surv_time_var_labels[y], ")") ,
                        pval.size = 4, xlab = "Time (years)", legend = c(0.7, 0.8),
                        legend.labs = legend_labs[[x]],
                        ylab = ylab[y], break.x.by = 2, break.y.by = 0.2)
    # Customize plot
    kmPlot$plot = kmPlot$plot + mytheme
  }) %>% 
    set_names(surv_time_var_labels)
}) %>% 
  set_names(strata_var_labels)



## Proteomic data ----------------------------------------------------------
# Proteomic data with annotation data
proteomic <- import_list(
  c(expression = "tmt29_processed_data_uncorrected.txt", 
    protein_group_annotation = "tmt29_protein_group_annotation.txt", 
    sample_annotation = "tmt29_sample_annotation.txt"))

names(proteomic) <- c("expression", "protein_group_annotation", "sample_annotation")

# Select previously analized 108 patients from expression data
proteomic$expression <- proteomic$expression %>% 
  dplyr::select(protein_group = ProteinGroup, one_of(clinical_data$aliquot_sample_id))

# RNA-seq data ------------------------------------------------------------
# # RNA-seq data with metadata
# rna_seq <- import_list(c(expression = "v2.1_rnaseq_iron_untilt_filt3_debatched.txt", 
#                          metadata = "v2.1_lusc_proteogenomics_rnaseq_metadata.xlsx"))

# New RNA-seq data with metadata
rna_seq <- import_list(
  c(expression = "v2.2_lusq116_rnaseq_iron_untilt_filt3_debatched.txt", 
    metadata = "v2.2_lusq_proteogenomics_rnaseq_metadata_with_batches_2018-05-09.xlsx"))

names(rna_seq) <- c("expression", "metadata")

# Filter metadata to include previously analyzed 108 patients and redone sample
rna_seq$metadata <- rna_seq$metadata %>% 
  clean_names() %>% 
  filter(tissue_id %in% clinical_data$aliquot_sample_id, 
         (redo == "Redo" | is.na(redo)))

# Select sample id from RNA-seq data that are in metadata only
rna_seq$expression <- rna_seq$expression %>% 
  dplyr::select(Symbol, one_of(rna_seq$metadata[["genomics_sample"]]))


# Targeted exome sequencing -----------------------------------------------
mutation <- read_delim("mutations_minGQ15_mapped_sample.txt", delim = "\t") %>% 
  filter(tumor_06s %in% clinical_data$aliquot_sample_id)

# Copy number  ------------------------------------------------------------
gene_cnv <- read_delim("Gene_CNV_400kbp.50counts_noXYgene_2018-02-12.txt", delim = "\t")

# Select previously analized 108 patients from expression data
gene_cnv <- gene_cnv %>% 
  dplyr::select(Gene.Symbol, one_of(str_c("X", clinical_data$aliquot_sample_id)))

# QC ----------------------------------------------------------------------
# Missing proportion and mean expression by gene in proteomic, RNA-seq and gene CNV data
missing_proportion <- list(proteomic$expression, rna_seq$expression, gene_cnv) %>% 
  map(~ .x %>% 
        mutate_if(is.numeric, funs(if_else(is.na(.), 0, as.numeric(.)))) %>%  
        mutate(missing_prop = rowMeans(dplyr::select(., -1) == 0),
               exprs_mean = rowMeans(dplyr::select(., -1))) %>% 
        dplyr::select(1, missing_prop, exprs_mean)) %>% 
  set_names(c("Proteomic", "RNA-seq", "Gene CNV"))

# export(missing_proportion, "missing_proportion.xlsx")

# Distribution of mutation status by gene
mutation_status_dist <- dplyr::select(mutation, -tumor_06s, -dna_sample) %>% 
  map_df(~ table(.x) %>% 
           enframe() %>% 
           spread(name, value), .id = "gene")

# export(mutation_status_dist, "missing_proportion.xlsx",  
#        which = "Mutation status", overwrite = FALSE)

# Data filtering ----------------------------------------------------------
# Min of proteomic expression
proteomic_exprs_min <- min(proteomic$expression[, -1], na.rm = TRUE)

# Keep the proteins with <10% missing across all samples in proteomic data and 
# impute the missing values using overall minimum detected positive protein level 
proteomic$expression <- proteomic$expression %>% 
  filter(rowMeans(is.na(select(., -protein_group))) < 0.1) %>% 
  mutate_at(vars(-protein_group), funs(coalesce(., proteomic_exprs_min)))

# Filter out genes in RNA-seq data with mean expression less than or 
# equal to 5 in normalized log2 scale
rna_seq$expression <- rna_seq$expression %>% 
  filter(rowMeans(select(., -Symbol), na.rm = TRUE) > 5)

# Mutation status for mutated vs wild type
mutation <- mutation %>% 
  mutate_at(
    vars(-tumor_06s), 
    funs(recode_factor(., "wt" = "Wild type", "mut" = "Mutated", "trunc" = "Mutated")))

# List of suggested gene (sent by Dr. Chen on 3/22/2018) to analyze in mutation data
mutation_genes <- c(
  "MLL2", "NFE2L2", "CDKN2A", "NOTCH1", "KEAP1", "PIK3CA", "PTEN",
  "RB1", "TP53", "LAMA2", "ERBB4", "CYP2B6", "APC", 'GNAS')

# Select genes in mutation_genes only
mutation <- mutation %>% 
  dplyr::select(symbol = tumor_06s, one_of(mutation_genes))

# Select only the genes with at least 10% samples have some gain/loss in 
# gene CNV data. Replace missing, if any, with 0
gene_cnv <- gene_cnv %>% 
  mutate_at(vars(-Gene.Symbol), funs(if_else(is.na(.), 0, as.numeric(.)))) %>% 
  filter(rowMeans(select(., -Gene.Symbol) == 0) <= 0.9) 

# Association between TLN and mutation for NRF2 (aka NFE2L2) and KEAP1 genes
# Clinical data with mutation status for those two genes
clinical_data <- left_join(clinical_data, 
                           select(mutation, symbol, NFE2L2, KEAP1),
                           by = c("aliquot_sample_id" = "symbol"))

tln_gene_association <- dplyr::select(clinical_data, NFE2L2, KEAP1) %>% 
  map_df(~ table(.x,  clinical_data$tln2) %>% 
           fisher.test() %>% 
           tidy(),
         .id = "Gene") %>% 
  dplyr::select(Gene, OR = estimate, `CI (low)` = conf.low, 
                `CI (high)` = conf.high, P = p.value)

tln_gene_table <- c(NFE2L2 = "NFE2L2", KEAP1 = "KEAP1") %>% 
  map_df(function(x) tabyl(clinical_data, !!rlang::sym(x), tln2) %>% 
           rename(Mutation = x),
         .id = "Gene")

# export(list(`Fisher test` = tln_gene_association, 
#             `Cross table` = tln_gene_table), "TLN & gene association.xlsx")

# Univariable survival analysis -------------------------------------------
# Proteomic expression data with survival informoation
proteomic_exprs_surv <- proteomic$expression %>% 
  gather(sample, exprs, -protein_group) %>% 
  spread(protein_group, exprs) %>% 
  left_join(select(clinical_data, aliquot_sample_id, overall_survival_2_years, os_status, 
                   recurrence_free_survival_2_years, rfs_status, tln2),
            by = c("sample" = "aliquot_sample_id"))

# RNA-seq expression data with survival informoation
rna_exprs_surv <- rna_seq$expression %>% 
  gather(sample, exprs, -Symbol) %>% 
  spread(Symbol, exprs) %>% 
  left_join(select(rna_seq$metadata, genomics_sample, tissue_id), 
            by = c("sample" = "genomics_sample")) %>% 
  left_join(select(clinical_data, aliquot_sample_id, overall_survival_2_years, 
                   os_status, recurrence_free_survival_2_years, rfs_status, tln2),
            by = c("tissue_id" = "aliquot_sample_id"))

# Mutation data with survival information
mutation_surv <- mutation %>% 
  left_join(select(clinical_data, aliquot_sample_id, overall_survival_2_years, 
                   os_status, recurrence_free_survival_2_years, rfs_status, tln2),
            by = c("symbol" = "aliquot_sample_id"))

# Gene CNV data with survival information
gene_cnv_surv <- gene_cnv  %>% 
  gather(sample, exprs, -Gene.Symbol) %>% 
  spread(Gene.Symbol, exprs) %>% 
  mutate(sample = str_replace(sample, "^X", "")) %>% 
  left_join(select(clinical_data, aliquot_sample_id, overall_survival_2_years, 
                   os_status, recurrence_free_survival_2_years, rfs_status, tln2),
            by = c("sample" = "aliquot_sample_id"))


# ##############################################################################
# # Fitting Cox models w/ and w/o adjusting by TLN for all four datasets
##############################################################################
# Parallelizing the codes for Cox fits -----------------------------------------
# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
cl <- makeCluster(no_cores)

# Load packages and export data to global env
clusterExport(cl, list("proteomic", "proteomic_exprs_surv",
                       "rna_exprs_surv", "mutation_surv", "gene_cnv_surv"))

clusterEvalQ(cl, list(library(survival), library(broom), library(dplyr)))

# Fit Cox model to proteomic data ----------------------------------------------
coxfit_exprs_proteomic <- parLapply(
  cl, proteomic$expression[["protein_group"]], function(x){
    tryCatch({coxph(
      as.formula(
        paste0("Surv(overall_survival_2_years, os_status) ~ ", "`", x, "`")),
      data = proteomic_exprs_surv) %>%
        tidy(exponentiate = TRUE)}, error = function(e) NULL)
  }) %>%
  ldply() %>%
  dplyr::select(protein_group = term, HR = estimate, `HR (low)` = conf.low,
                `HR (high)` = conf.high, z = statistic, P = p.value) %>%
  # adjusted p-value (FDR)
  mutate(beta = log(HR),
         P_fdr = p.adjust(P, method = "fdr"),
         Qvalue = qvalue(P)$qvalues,
         protein_group = as.integer(str_replace_all(protein_group, "`", "")))

# saveRDS(coxfit_exprs_proteomic, file = "coxfit_exprs_proteomic.RDS")

stopCluster(cl)

# Do the same for other datasets

##############################################################################


# # Load clinical and proteomic datasets
# load("clinical_protemoic_data.RData")

# Saved Cox fit RDS filenames of all 4 datasets
coxfit_filenames <- list.files(pattern = "coxfit.*\\.RDS")

# Read the .RDS files and assign to global env
map2(str_replace(coxfit_filenames, ".RDS", ""), coxfit_filenames, ~ assign(.x, readRDS(.y), env = .GlobalEnv))

# # Mean and SD of proteomic expression data
# proteomic_exprs_mean <- proteomic$expression  %>% 
#   mutate(mean = rowSummary(dplyr::select(., -protein_group), FUN = mean, na.rm = TRUE),
#          sd = rowSummary(dplyr::select(., -protein_group), FUN = sd, na.rm = TRUE)) %>% 
#   dplyr::select(protein_group, mean, sd)

# Add annotation data to fitted Cox models
coxfit_exprs_proteomic <- coxfit_exprs_proteomic %>% 
  map(~ .x %>% 
        left_join(select(missing_proportion$Proteomic, -exprs_mean), by = "protein_group") %>% 
        # left_join(proteomic_exprs_mean, by = "protein_group") %>% 
        left_join(proteomic$protein_group_annotation, by = c("protein_group" = "ProteinGroup")))

# Adjusted by tln
coxfit_exprs_proteomic_tln <- coxfit_exprs_proteomic_tln %>% 
  map(~ .x %>% 
        left_join(select(missing_proportion$Proteomic, -exprs_mean), by = "protein_group") %>% 
        # left_join(proteomic_exprs_mean, by = "protein_group") %>% 
        left_join(proteomic$protein_group_annotation, by = c("protein_group" = "ProteinGroup"))) 

# # Mean and SD of RNA-seq expression data
# rna_exprs_mean <- rna_seq$expression  %>% 
#   mutate(mean = rowSummary(dplyr::select(., -Symbol), FUN = mean, na.rm = TRUE),
#          sd = rowSummary(dplyr::select(., -Symbol), FUN = sd, na.rm = TRUE)) %>% 
#   dplyr::select(Symbol, mean, sd)

# Add mean and missing proportion to coxfit_exprs_rna
coxfit_exprs_rna <- coxfit_exprs_rna %>% 
  map(~ .x %>% 
        left_join(select(missing_proportion$`RNA-seq`, -exprs_mean), by = c("gene" = "Symbol")))

# Add mean and missing proportion to coxfit_exprs_rna_tln
coxfit_exprs_rna_tln <- coxfit_exprs_rna_tln %>% 
  map(~ .x %>% 
        left_join(select(missing_proportion$`RNA-seq`, -exprs_mean), by = c("gene" = "Symbol")))

# Add number of mutated/Wild type to coxfit_exprs_mutation
coxfit_exprs_mutation <- coxfit_exprs_mutation %>% 
  map(~ .x %>% 
        mutate(gene2 = word(gene, sep = "- ")) %>% 
        left_join(mutation_status_dist, -exprs_mean, by = c("gene2" = "gene")) %>% 
        dplyr::select(-gene2))

# Add number of mutated/Wild type to coxfit_exprs_mutation_tln
coxfit_exprs_mutation_tln <- coxfit_exprs_mutation_tln %>% 
  map(~ .x %>% 
        mutate(gene2 = word(gene, sep = "- ")) %>% 
        left_join(mutation_status_dist, -exprs_mean, by = c("gene2" = "gene")) %>% 
        dplyr::select(-gene2))

# Add missing proportion to coxfit_exprs_cnv
coxfit_exprs_cnv <- coxfit_exprs_cnv %>% 
  map(~ .x %>% 
        left_join(select(missing_proportion$`Gene CNV`, -exprs_mean), by = c("gene" = "Gene.Symbol")))

# Add missing proportion to coxfit_exprs_cnv_tln
coxfit_exprs_cnv_tln <- coxfit_exprs_cnv_tln %>% 
  map(~ .x %>% 
        left_join(select(missing_proportion$`Gene CNV`, -exprs_mean), by = c("gene" = "Gene.Symbol")))

#######################################################################################################
## Meta-analysis 
# Overall survival --------------------------------------------------------
# Proteomic and rna-seq
proteomic_rna_os <- full_join(
  coxfit_exprs_proteomic$OS, coxfit_exprs_rna$OS, 
  by = c("Symbol" = "gene"), suffix = c("_proteomic", "_rna")) 

## Proteomic, rna-seq and CNV
proteomic_rna_cnv_os <- reduce(
  list(coxfit_exprs_proteomic$OS, coxfit_exprs_rna$OS, coxfit_exprs_cnv$OS),
  full_join, by = c("Symbol" = "gene"))


# RFS ---------------------------------------------------------------------
# Proteomic and rna-seq
proteomic_rna_rfs <- full_join(
  coxfit_exprs_proteomic$RFS, coxfit_exprs_rna$RFS, 
  by = c("Symbol" = "gene"), suffix = c("_proteomic", "_rna")) 

## Proteomic, rna-seq and CNV
proteomic_rna_cnv_rfs <- reduce(
  list(coxfit_exprs_proteomic$RFS, coxfit_exprs_rna$RFS, coxfit_exprs_cnv$RFS),
  full_join, by = c("Symbol" = "gene"))

### Reuse Zach's codes for meta-analysis
proteomic_rna_cnv_os<-filter(proteomic_rna_cnv_os,!is.na(beta.x)&!is.na(beta.y)&!is.na(beta))

proteomic_rna_cnv_os$se.beta.x<-proteomic_rna_cnv_os$beta.x/proteomic_rna_cnv_os$z.x

proteomic_rna_cnv_os$se.beta.y<-proteomic_rna_cnv_os$beta.y/proteomic_rna_cnv_os$z.y

proteomic_rna_cnv_os$se.beta<-proteomic_rna_cnv_os$beta/proteomic_rna_cnv_os$z

metafunc<-function(x){
  as.vector(as.numeric(x[c( "beta.x" ,"beta.y","beta")]))->testy;
  print(testy);
  testy<-testy[!is.na(testy)];
  as.vector(as.numeric(x[c( "se.beta.x","se.beta.y","se.beta")]))->testse;
  testse<-testse[!is.na(testse)];
  random<-rma( yi = testy, sei= testse, method="EB");
  protein_group<-x["protein_group"];
  Symbol<-x["Symbol"];
  data.frame(protein_group=protein_group,Symbol=Symbol,MetaBeta=random$b,MetaBeta_Pvalue=random$pval,
             datav=collapse(as.character(testy),sep="|"))->output;
  return(output);
}#

meta0<-apply(proteomic_rna_cnv_os,1,metafunc);

metaBetaandPvalue<-do.call(rbind.data.frame, meta0 )#rbind.data.frame

# FDR adjusted Fisher's pvalue
metaBetaandPvalue <- metaBetaandPvalue %>% 
  mutate(MetaBeta_Fisher_p_fdr = as.numeric(p.adjust( MetaBeta_Pvalue, method = "fdr")),
         MetaBeta_Qvalue = as.numeric(qvalue(MetaBeta_Pvalue)$qvalues))

proteomic_rna_cnv_os_meta<-metaBetaandPvalue;

test<-cbind(proteomic_rna_cnv_os,proteomic_rna_cnv_os_meta )
# test<-left_join(proteomic_rna_cnv_os,proteomic_rna_cnv_os_meta, by = c("protein_group", "Symbol"))

proteomic_rna_cnv_os_meta_unadjusted_overlapping<-test


# Forest plot of HR (with 95% CI) by genes for OS (unadjusted)
proteomic_rna_cnv_os_forest <- proteomic_rna_cnv_os_meta_unadjusted_overlapping %>% 
  set_tidy_names() %>% 
  select(protein_group = protein_group..1, Symbol = Symbol..18, 
         MetaBeta_Qvalue, matches("HR$|HR\\.|HR\\s")) %>% 
  filter(MetaBeta_Qvalue <= 0.3) %>% 
  distinct(Symbol, .keep_all = TRUE) %>% 
  gather(hr_key, hr_value, HR.x:`HR (high)`) %>% 
  mutate(hr_value = if_else(hr_value > 3.5, 3.5, hr_value),
         data = case_when(str_detect(hr_key, ".x$") ~ "Proteomic",
                          str_detect(hr_key, ".y$") ~ "RNA-seq",
                          TRUE ~ "CNV"),
         hr_key = str_remove(hr_key, ".x$|.y$")) %>% 
  spread(hr_key, hr_value) %>% 
  ggplot(aes(x = data, y = HR, ymin = `HR (low)`,
             ymax = `HR (high)`, color = data)) +
  geom_pointrange(size = 0.2) +
  coord_flip() +
  scale_color_manual(values = c("blue", "red", "gray40"),
                     labels = c("CNV", "Proteomic", "RNA-seq")) +
  ylab("Hazard ratio (95% CI)") +
  ggtitle("OS (unadjusted)") +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ Symbol, strip.position = "left", nrow = 16, scales = "free_y") +
  # theme_bw() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        strip.text.y  =  element_text(hjust = 0, vjust  =  0.5, angle = 180, face = "bold"))

# ggsave(filename = "forest_hr_meta_OS.pdf", plot = proteomic_rna_cnv_os_forest, height = 10, width = 8)

