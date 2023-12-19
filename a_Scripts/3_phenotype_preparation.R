#---------------------------------------- PHEWAS --------------------------------------------------#

# UK BIobank is a repository with ~ 500000 samples (genomic and general data from healthy and diseases 
# participants of all kind). 

# PheWAS analysis takes SNPs and associated them with different traits or diseases.

# This scripts is to prepare the table that will be the phenotype table for the PheWAS analysis. The 
# table format is as follow: ID, vocabulary ID, code, count (how many times is reported - different
# dates).

# ----> LIBRARY LOAD:

if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')

# ---> SET DIRECTORY:

# Set directory to be in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ----> INPUT:

# Since the files are too big, I prefered to keep them in the output table and call them from there. 
eu_table <- read.delim("../../output/2_sample_selection/eu_table_pca.tab") # Save table with only codes. Less heavy.
sa_table <- read.delim("../../output/2_sample_selection/sa_table_pca.tab")

table_list <- list(eu_table, sa_table)
name_list <- c("eu_phenotype", "sa_phenotype")

# ---> DATA MANIPULATION:

for (i in 1:length(table_list)) {

# Get the table for each population: 
pop_table <- table_list[[i]]

# Identify codes for the columns wit want to put together in a single column
codes <- c("41203", "41202", "20001", "20002", "40001", "40002", "40006", "40013")
vocabulary <- c("ICD9CM", "ICD10CM", "Cancer", "Non_cancer", "ICD10D", "ICD10D2", "ICD10Cancer", "ICD9Cancer")

# Get a table with only the columns of interest: 
code_table <- pop_table[, grepl(paste(codes, collapse="|"), colnames(pop_table)) | colnames(pop_table) == "f.eid"]
code_table <- code_table %>% mutate_if(is.numeric, as.character)

# Get a table with the covariates: 

cov_table <- pop_table[, c(1,533:545)]

# Output cov table:

write.table(cov_table, 
            file = paste("../../output/3_phenotype_file_preparation/cov_", name_list[[i]], "_pca.tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

  
# Get all the values for all the columns except the first one in a single column with the identifier of
# the first column. 
id_voc_count <- code_table %>% pivot_longer(cols = colnames(code_table)[-1],
                                                  names_to = "vocabulary_id", # Where the column names will go.
                                                  values_to = "code") # Where the values will go. 

# pivot longer makes sure of taking all the values from the different columns into one. 

# Add a column with 1, so we can do the count later: 
id_voc_count$count <- 1

#Create control table for all the vocabularies: 
controls_total <- data.frame()

# This loop replace all the values in vocabulary id columns with the new names from vocabulary vector.
for (k in seq_along(codes)) {
  
  id_voc_count$vocabulary_id[grepl(codes[[k]], id_voc_count$vocabulary_id)] <- vocabulary[[k]]
  
  # Save the columns with NA for each group as possible controls: 
  
  controls <- id_voc_count[id_voc_count$vocabulary_id == vocabulary[[k]], ] # We have to do it for every vocabulary.
  controls_no_dupl <- controls[!duplicated(controls), ] # Delete duplicated rows
  
  # Get a list of the samples with none medical record (we want controls without any condition for that
  # specific code): 
  
  cases <- unique(controls_no_dupl[!is.na(controls_no_dupl$code), ]$f.eid)
  
  # Get the samples without any condition
  final_controls <- controls_no_dupl[!(controls_no_dupl$f.eid %in% cases), ]
  final_controls$count <- NA
  
  controls_total <- rbind(controls_total, final_controls)
  
}

# Write the output table for the controls:

write.table(controls_total, 
            file = paste("../../output/3_phenotype_file_preparation/Controls_", name_list[[i]], ".tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Delete all the columns with NA or 99999 values. 
id_voc_count <- id_voc_count[!is.na(id_voc_count$code), ]
id_voc_count <- id_voc_count[id_voc_count$code != "99999", ]

# We need to aggregate all the codes that are repeat for the same patient: 

id_voc_count <- aggregate(.~f.eid + vocabulary_id + code, data = id_voc_count, FUN = sum)

# Export results: 

write.table(id_voc_count, 
            file = paste("../../output/3_phenotype_file_preparation/", name_list[[i]], ".tab", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

}


