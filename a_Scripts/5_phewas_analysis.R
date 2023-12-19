#---------------------------------------- PHEWAS --------------------------------------------------#

# UK BIobank is a repository with ~ 500000 samples (genomic and general data from healthy and diseases 
# participants of all kind). 

# PheWAS analysis takes SNPs and associated them with different traits or diseases.

# This scripts used the pheWAS library to run the PheWAS analysis. 

# ----> LIBRARY LOAD:

# To install PheWAS
# devtools::install_github("PheWAS/PheWAS") 
library("PheWAS")
library('stringr')
library('ukbtools')

# ---> SET DIRECTORY:

# Set directory to be in the source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ----> INPUT:

# Coding ICD10: 

code_icd10 <- read.delim("coding_ICD10.tsv")
names(code_icd10)[1] <- "code"
code_icd10[c("code_2", "meaning")] <- str_split_fixed(code_icd10$meaning, ' ', 2)
code_ic_d10_men <- code_icd10[ , c(1,2,6)]
code_icd10 <- code_icd10[ ,c(1,6)]


# Coding ICD9:

code_icd9 <- icd9codes
names(code_icd9)[1] <- "code"
code_icd9[c("code_2", "meaning")] <- str_split_fixed(code_icd9$meaning, ' ', 2)
code_icd9_men <- code_icd9
code_icd9 <- code_icd9[ ,c(1,3)]

# European data:

eu_phe <- read.delim("eu_phenotype.tab")
names(eu_phe)[1]="id"
eu_gen <- read.table("eu_genotypes.raw",header=TRUE)[,c(-2:-6)]
names(eu_gen)[1]="id"
eu_con <- read.delim("Controls_eu_phenotype.tab")
names(eu_con)[1]="id"
eu_cov <- read.delim("cov_eu_phenotype.tab")
names(eu_cov)[1]="id"
eu_cov[,c(3:4)] <- lapply(eu_cov[,c(3:4)], factor)

# Asian data:

#sa_phe <- read.delim("sa_phenotype.tab")
#names(sa_phe)[1]="id"
#sa_gen <- read.table("sa_genotypes.raw",header=TRUE)[,c(-2:-6)]
#names(sa_gen)[1]="id"
#sa_con <- read.delim("Controls_sa_phenotype.tab")
#names(sa_con)[1]="id"
#sa_cov <- read.delim("cov_sa_phenotype.tab")
#names(sa_cov)[1]="id"
#sa_cov[,c(3:4)] <- lapply(sa_cov[,c(3:4)], factor)

# ----> DATA TRANSFORMATION:

# Total Samples:
eu_phe_total <- rbind(eu_phe, eu_con)


icds <- c("ICD9CM", "ICD10CM", "ICD10D", "ICD10D2", "ICD10Cancer", "ICD9Cancer")
icds_voc <- c("ICD9CM", "ICD10CM", "ICD10CM", "ICD10CM", "ICD10CM", "ICD9CM")
snps <- c("rs61823694_T", "rs4266947_T" )

for (i in 1:length(icds)) {
  
  # PHENOTYPE:
  
  # Get the vocabulary for the codes that I will analyse:
  phe <- eu_phe_total[eu_phe_total$vocabulary_id == icds[i], ]
  
  # Get the correct coding since UK Biobank deletes the decimals in the ICD codes:
  if (grepl(9, icds[i]) == TRUE) {
    
    phe_code <- merge(phe, code_icd9, all.x = TRUE) # Merge with the right coding table
    
  } else {
    
    phe_code <- merge(phe, code_icd10, all.x = TRUE)
    
  }
  
  phe_code <- phe_code[, c(2,3,5,4)] # Re-order table
  names(phe_code)[3] <- "code" # Rename the coding column with the right name
  phe_code$vocabulary_id <- icds_voc[i]
  phe_code$vocabulary_id <- factor(phe_code$vocabulary_id) # Get vocabulary column as factor
    
  # GENOTYPE: 
  
  for (j in snps) {
    
    # Get the genotype table only for one SNP at the time: 
    
    gen <- eu_gen[, c("id", j)]
    
    # Get the right phenotype format to run PheWAS:
    
    # A little confused about the min.code.count but it looks like the package consider that a person
    # has a condition if the diagnose has been done more than one. Here I put 1 since I want to have
    # all the conditions. id.sex gets the phenotype by sex. 
    phenotype_count <- createPhenotypes(phe_code, min.code.count = 1, id.sex = eu_cov[, c(1,4)],
                                        aggregate.fun = sum) # Dont think aggregate fun do something
    
    # Here, instead of considering the phenotypes as cases and controls, they became continuous variables,
    # indicating number of times that the phenotype is present is relevant to do the association analysis.
    # I don't think is the right way to do this, but still, good to have it. 
    phenotype_na <- createPhenotypes(phe_code, min.code.count = NA, id.sex = eu_cov[, c(1,4)],
                                     aggregate.fun = sum)
    
    # Run the PheWAS analysis:
    
    # Results for the the two types of created phenotypes. Covariates include age and assesment centre.
    # Sex is already add in the phenotype file. Min.records is the minimum number of cases and controls 
    # to calculate the association. 
    results_count <- phewas(phenotypes=phenotype_count, genotypes=gen, covariates = eu_cov[, c(1:3)],
                            significance.threshold = c("p-value", "bonferroni", "fdr"), min.records = 15)
    
    results_na <- phewas(phenotypes = phenotype_na, genotypes = gen, covariates = eu_cov[, c(1:3)],
                         significance.threshold = c("p-value", "bonferroni", "fdr"), min.records = 15)
    
    # Add descriptions: 
    
    results_count_d <- addPhecodeInfo(results_count)
    results_na_d <- addPhecodeInfo(results_na)
    
    # Output tables:
    
    write.table(results_count_d, 
                file = paste("../../output/5_phewas/EU_Count_", icds[i], "_", j ,".tab", sep = ""),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)

    write.table(results_na_d, 
                file = paste("../../output/5_phewas/EU_Cont_", icds[i], "_", j ,".tab", sep = ""),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
    
    # Generate figures: 
    
    if (all(is.na(results_count$p.value)) == FALSE) {
    
    # Categorical variable:
    
    man_count <- phewasManhattan(results_count, OR.direction = TRUE, suggestive.line = -log10(0.05),
                                 significant.line = -log10(0.05))
    
    png(file = paste("../../output/5_phewas/EU_Count_", icds[i], "_", j ,".png", sep = ""),
        width = 1200, height = 900)
    
    print(man_count)
    
    dev.off()
    
    }
    
    # Continuous variable:
    
    if (all(is.na(results_na$p.value)) == FALSE) {
    
    man_na <- phewasManhattan(results_na, suggestive.line = -log10(0.05), significant.line = -log10(0.05))
    
    png(file = paste("../../output/5_phewas/EU_Cont_", icds[i], "_", j ,".png", sep = ""),
        width = 1200, height = 900)
    
    print(man_na)
    
    dev.off()
    
    }
    
  }
}




