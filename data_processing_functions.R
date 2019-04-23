library(stringr)
library(data.table)

#Parses the PTM position in MaxQuant output with PTMs. Uses a "probability_threshold" to arbitrarily flag high probability modifications.
get_ptm_position = function(sequence_in,protein_positions_in,peptide_position_in,probability_threshold) {
    sequence_start = as.numeric(str_split(protein_positions_in,";")[[1]][1]) - as.numeric(str_split(peptide_position_in,";")[[1]][1]) 
    sequence_location = 0
    ptm_site_location_peptide = numeric()
    sequence_in_split = unlist(str_split(sequence_in,""))
    for (i in 1:length(sequence_in_split)) {
        if (str_detect(sequence_in_split[i],"[A-Z]")) {
            sequence_location = sequence_location + 1
        }
        else if (str_detect(sequence_in_split[i],"\\(")) {
            ptm_site_location_peptide = c(ptm_site_location_peptide,sequence_location)
        }
    }
    aa_sequence = unlist(str_split(paste(unlist(str_split(sequence_in,"\\(([0-9]*\\.[0-9]*)\\)|\\(([0-9])\\)")),collapse=""),"")) #split by the probabilities in parenthesis and extract out the amino acids
    aa = aa_sequence[ptm_site_location_peptide] #get the modified amino acid
    modification_probs = as.numeric(unlist(str_extract_all(sequence_in,"([0-9]*\\.[0-9]*)|([0-9])")))
    ptm_site_location_protein = sequence_start + ptm_site_location_peptide
    modification_df = data.frame(aa=aa,ptm_site_location_peptide=ptm_site_location_peptide,ptm_site_location_protein=ptm_site_location_protein,modification_probs = modification_probs,stringsAsFactors=FALSE)
    modification_string = character()
    for (i in 1:nrow(modification_df)) {
        if (as.numeric(modification_df$modification_probs[i]) >= probability_threshold) {
            modification_string = c(modification_string,paste(c(modification_df[i,1],modification_df[i,3],"*"),collapse=""))
        } else {
            modification_string = c(modification_string,paste(c(modification_df[i,1],modification_df[i,3]),collapse=""))
        }
    }
    return(paste(modification_string,collapse=","))
}

#Average technical replicates together, but don't consider missing values in the average.
average_replicates = function(df_in, metadata_in, missing_value, log2_output = TRUE) {
    sample_column = "sample"
    length_unique_sample_column = length(unique(metadata_in[, sample_column]))
    averaged_expression_matrix = matrix(data = NA, nrow = dim(df_in)[1], ncol = length_unique_sample_column)
    for (i in 1:length_unique_sample_column) {
        current_sample = unique(metadata_in[, sample_column])[i]
        current_sample_full_text = as.character(metadata_in[metadata_in[, sample_column] == current_sample, "expression_column"])
        current_sample_df = df_in[, current_sample_full_text]
        for (j in 1:nrow(current_sample_df)) {
            if (is.na(missing_value)) {
                current_peptide_missing = sum(is.na(current_sample_df[j, ]))
            } else {
                current_peptide_missing = sum(current_sample_df[j, ] == missing_value)
            }
            current_peptide_sum = sum(na.omit(as.numeric(current_sample_df[j, ])))
            current_sample_df_ncol = ncol(current_sample_df)
            if (current_peptide_missing < current_sample_df_ncol) {
                averaged_expression_matrix[j, i] = current_peptide_sum/(current_sample_df_ncol - current_peptide_missing) # "average" the nonzero expression
            } else if (current_peptide_missing == current_sample_df_ncol) {
                averaged_expression_matrix[j, i] = NA
            }
        }
    }
    colnames(averaged_expression_matrix) = unique(metadata_in[, sample_column])
    if (log2_output == TRUE) {
        return(log2(data.frame(averaged_expression_matrix, stringsAsFactors = FALSE)))   
    }
    else {
        return(data.frame(averaged_expression_matrix, stringsAsFactors = FALSE))
    }
}


process_mq_data = function(df_in, protein_id_col, metadata_in, is_ptm_data = FALSE, ptm_threshold = 0.50, IRON = TRUE, find_median_col = "", missing_values = 0, avg_replicates = TRUE, run_peptide_qc = TRUE) {
    if (sum(c("expression_column", "sample") %in% names(metadata_in)) != 2) {
        stop("Metadata is missing the column 'expression_column' or 'sample'")
    }
    if (avg_replicates == TRUE) {
        if (any("replicate" %in% names(metadata_in)) == FALSE) {
            stop("Metadata is missing the column 'replicate'")
        }
    }
    if (is_ptm_data == TRUE) {
        modification_prob_column = names(df_in)[str_detect(names(df_in), "Probabilities")]
    }
    
    #normalize data with IRON ----
    if (IRON == TRUE) {
        if (find_median_col == "") {
            ironed_df = iron_proteomics(df_in[, metadata_in$expression_column])
            metadata_in$pre_iron_expression = metadata_in$expression_column
            metadata_in$expression_column = names(ironed_df)
            df_in = cbind(df_in, ironed_df)
        } else {
        ironed_df = iron_proteomics(df_in[, metadata_in$expression_column], median_sample = find_median_col)
        metadata_in$pre_iron_expression = metadata_in$expression_column
        metadata_in$expression_column = names(ironed_df)
        df_in = cbind(df_in, ironed_df)
        }
    }
    
    #peptide_qc ----
    if (run_peptide_qc == TRUE) {
        #flag non-human entries
        uniprot_entries = read.delim("~/Google Drive/data/SwissProt_Human_2018_05_list_of_uniprot_2018-07-11.txt", header = T, stringsAsFactors = FALSE)
        df_in$qc_nothuman = sapply(df_in[, protein_id_col], function(x) {
            sum(unlist(str_extract_all(x, uniprot_regex(isoforms = "none"))) %in% uniprot_entries$Entry) == 0
        })
        
        #flag reverse only entries
        df_in$qc_reverse_only = sapply(df_in[, protein_id_col], function(x) {
            x_split = unlist(str_split(x, ";"))
            if (sum(str_detect(x_split, "REV__(CON__){0,1}|(REV__){0,1}CON__ENSEMBL|(REV__){0,1}CON__REFSEQ|(REV__){0,1}CON__H-INV")) == length(x_split)) {
                TRUE
            } else {
                FALSE
            }
        })
        
        #flag completely missing rows; na.omit() works for SILAC ratios because sum(na.omit(c(NA,NA,NA))) is 0.
        df_in$qc_missing_all_expression = apply(df_in, 1, function(x) {
            sum(na.omit(as.numeric(x[metadata_in$expression_column]))) == 0
        })
        
        #flag PEP score greater than threshold
        df_in$qc_pep = sapply(df_in$PEP, function(x) {
            if (x > 0.1) {
                TRUE
            } else {
                FALSE
            }
        })
        
        if (is_ptm_data == TRUE) {
            #flag missing modification count
            df_in$qc_mod_missing = sapply(df_in[, modification_prob_column], function(x) {
                if (x == "" | x == 0) {
                    TRUE
                } else {
                    FALSE
                }
            })
        }
        
        #get rid of rows that don't pass all qc
        if (is_ptm_data == TRUE) {
            qc_columns = c("qc_nothuman", "qc_reverse_only", "qc_missing_all_expression", "qc_pep", "qc_mod_missing")
            
        } else if (is_ptm_data == FALSE) {
            qc_columns = c("qc_nothuman", "qc_reverse_only", "qc_missing_all_expression", "qc_pep")
            
        }
        df_in$qc_pass = apply(df_in[, qc_columns], 1, function(x) { #weird things occurring with apply coersion; specify qc_columns in apply so that the booleans are read in as booleans and not characters
            sum(as.numeric(x[qc_columns]))
        })
        df_in_keep = subset(df_in, qc_pass == 0)
    } else {
        df_in_keep = df_in
    }
    
    #calculate number of PTMs for each site ----
    if (is_ptm_data == TRUE) {
        df_in_keep$PTM_position = ""
        for (i in 1:nrow(df_in_keep)) {
            df_in_keep$PTM_position[i] = get_ptm_position(df_in_keep[i, modification_prob_column], df_in_keep$Positions[i], df_in_keep$Position.in.peptide[i], ptm_threshold)
        }
    }
    
    #"average" expression by replicate but ignore missing values
    if (avg_replicates == TRUE) {
        combined_expression = average_replicates(df_in_keep, metadata_in, missing_values) #returns log2 expression
        return(list(df_out_processed = cbind(df_in_keep, combined_expression), metadata_out_processed = metadata_in))
    } else {
        df_in_keep[, metadata_in$expression_column] = log2(df_in_keep[, metadata_in$expression_column])
        df_in_keep[df_in_keep == -Inf] = NA
        return(list(df_out_processed = df_in_keep, metadata_out_processed = metadata_in))
    }
}

#MaxQuant no longer provides PEP scores in protein-level output. These need to be calculated as a product of the PEP scores of all peptides belonging to a particular protein group.
get_pep_from_peptides = function(peptide_ids_to_map, peptide_df) {
    pepscores = matrix(peptide_df[, "PEP"])
    row.names(pepscores) = peptide_df[, "id"]
    calculated_pep_score = vapply(peptide_ids_to_map, FUN.VALUE = numeric(1), function(x) {
        prod(pepscores[unlist(str_split(x, ";")),])
    })
    return(calculated_pep_score)
}

#Function for processing microarray data and finding the probe to a particular gene with the highest average intensity.
get_best_probe = function(df_in, gene_name_col, metadata_in, metadata_in_expression_colnames) {
    if (!metadata_in_expression_colnames %in% names(metadata_in)) {
        stop("Metadata is missing the specified expression column")
    }
    all_best_genes = list()
    unique_gene_names = unique(df_in[, gene_name_col])
    for (i in 1:length(unique_gene_names)) {
        current_gene = unique_gene_names[i]
        current_gene_df = df_in[df_in[, gene_name_col] == current_gene, ] 
        best_gene = -Inf
        for (j in 1:nrow(current_gene_df)) {
            current_gene_df_rowmean = mean(as.numeric(current_gene_df[j, metadata_in[, metadata_in_expression_colnames]]))
            if (current_gene_df_rowmean > best_gene) {
                best_gene = current_gene_df_rowmean
                best_gene_row = j
            }
        }
        all_best_genes[[i]] = data.frame(current_gene_df[best_gene_row, ])
    }
    return(as.data.frame(rbindlist(all_best_genes)))
}

#Function for binning a vector and coloring according to a blue-white-red gradient
center_bin_color = function(vector_in) {
    vector_in_standardized = scale(vector_in, center = TRUE, scale = TRUE)
    #floor and ceiling values at +/- 2
    vector_in_standardized[vector_in_standardized > 2] = 2
    vector_in_standardized[vector_in_standardized < -2] = -2
    brks = seq(-2, 2, length.out = length(vector_in_standardized))
    groups = cut(vector_in_standardized, breaks = brks, include.lowest = TRUE)
    vector_out = bluered(length(vector_in_standardized))[groups]
    vector_out[is.na(vector_out)] = "darkgrey"
    return(vector_out)
}

