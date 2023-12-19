#---------------------------------------- PHEWAS --------------------------------------------------#

# UK BIobank is a repository with ~ 500000 samples (genomic and general data from healthy and diseases 
# participants of all kind). 

# PheWAS analysis takes SNPs and associated them with different traits or diseases.

# This  script is used to get the final table, adding the assesment centre information (field 54) and
# subsetting the table by populations. An specific one for Europeans and another one for South Asians.

# ----> LIBRARY LOAD:

if (!require('stringr')) install.packages('stringr'); library('stringr')

# ---> SET DIRECTORY:

# Set directory to be in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ----> INPUT:

main_table <- read.delim("column_subset_phewas.tab") 
centre_table <- read.delim("ukb43914.tab")
pca <- read.delim("column_subset_pca.tab") # PCA info
pca_NA <- pca[!is.na(pca$f.22009.0.1), ]
colnames(pca_NA)[c(2:11)] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
withdrawal <- read.csv("w55718_20220222.csv", header = FALSE)
# exclude_eu <- read.delim("exclude_samples_eu.tab", header = FALSE) # Relative samples
# exclude_sa <- read.delim("exclude_samples_sa.tab", header = FALSE) # Relative samples 

# ---> DATA MANIPULATION:

# Merge the main table with the centre so we have everything together:

final_table <- merge(main_table, centre_table)
final_table <- merge(final_table, pca_NA) # To add the PC values

# Exclude the samples that withdrawn from UK Biobank:

final_table <- final_table[!(final_table$f.eid %in% withdrawal$V1), ]

# Get the age:

final_table$age <- 2023 - final_table$f.34.0.0

# Get the center: 

final_table$centre <- paste(final_table$f.54.0.0, final_table$f.54.1.0, 
                            final_table$f.54.2.0, final_table$f.54.3.0)

# Sex:

final_table$sex <- "M"
final_table$sex <- ifelse(!is.na(final_table$f.31.0.0) & final_table$f.31.0.0 == 0, "F", final_table$sex)


final_table$centre <- gsub("NA", "", final_table$centre)
final_table$centre <- gsub(" ", "", final_table$centre) # Remove the white spaces
final_table$centre <- strtrim(final_table$centre, 5)

# Separate by populations:

eu_table <- subset(final_table, 
                     (final_table$f.21000.0.0 == 1001 | final_table$f.21000.1.0 == 1001 | final_table$f.21000.2.0 == 1001) & 
                     final_table$f.22101.0.0 == "#ukbgene") 

sa_table <- subset(final_table, 
                   (final_table$f.21000.0.0 == 3001 | final_table$f.21000.0.0 == 3002 | final_table$f.21000.0.0 == 3003 | 
                      final_table$f.21000.1.0 == 3001 | final_table$f.21000.1.0 == 3002 | final_table$f.21000.1.0 == 3003 |
                      final_table$f.21000.2.0 == 3001 | final_table$f.21000.2.0 == 3002 | final_table$f.21000.2.0 == 3003) & 
                     final_table$f.22101.0.0 == "#ukbgene")

# Remove duplication and related samples: 

# eu_table_no_r <- eu_table[!(eu_table$f.eid %in% exclude_eu$V1), ]
# sa_table_no_r <- sa_table[!(sa_table$f.eid %in% exclude_sa$V1), ]

# ---> OUTPUT FILES:

write.table(eu_table, 
            file = "../../output/2_sample_selection/eu_table_pca.tab",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

writeLines(as.character(sa_table[, 1]), "../../output/2_sample_selection/eu_list_pca.txt", sep = "\n")

write.table(sa_table, 
            file = "../../output/2_sample_selection/sa_table_pca.tab",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

writeLines(as.character(sa_table[, 1]), "../../output/2_sample_selection/sa_list_pca.txt", sep = "\n")

