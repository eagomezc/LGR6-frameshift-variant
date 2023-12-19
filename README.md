# Loss of LGR6 expression as a result of a frameshift mutation alters pro-resolving responses in human phagocytes

# Overview: 

This repository contains the **Bash and R scripts** used to study the association between the genetic variant in LGR6, **rs74355478**, with a large number of phenotypes in a population level. Genotype and phenotype data from participant was obtained from the [**UK Biobank dataset**](https://www.ukbiobank.ac.uk/), variant identification were performed using [**Plink 1.9 tools**](https://www.cog-genomics.org/plink/) and phenome-wide association analysis was performed using the R script [**PheWAS**](https://github.com/PheWAS/PheWAS).

To identify genetic variants in **linkage of disequilibrium (LD)** with rs74355478, I used [**ensembl Grch37 website tool**](https://grch37.ensembl.org/Homo_sapiens/Tools/LD?db=core;tl=Gc01pnpSuyUdvsEo-8926382). More details of the parameters used can be found in the thesis.

**NOTE:** **PhenoScanner** (more information [here](http://www.phenoscanner.medschl.cam.ac.uk/) was used to identify previously reported human genotype-phenotype associations of rs74355478 and all the genetic variants in LD (r2 > 0.8).

# System Requirements: 

## Hardware requirements: 

All the scripts and software used for the **PhD Thesis** were run in a standard computer (RAM: 8GB, CP$: 4 cores, 3.60 GHZ/core). 

Big data analysis (including UK Biobank data manipulation and association analysis) were performed in the high-performance computing cluster of QMUL (more information [here](https://docs.hpc.qmul.ac.uk/)). Running time in the cluster depends of memory availability. 

## System requirements:

All the scripts were created and used on **Windows 10**:

**R version**: 4.0.4 

**R Studio version**: 1.1.524

The scripts should be compatible with Mac and Linux operating systems. 

For installing R and R Studio, follows the installation instructions [here](https://www.stats.bris.ac.uk/R/) and [here](https://www.rstudio.com/products/rstudio/download/). With typical installation times on a normal computer not exceeding 3h.

Plink 1.9 tools are already installed in high-performance computing cluster and it was called using the function **module load plink/1.9-170906**.

## Required R packages (libraries): 

The requiered packates to run all scripts should be installed automatically when you run each scritp; however, in case you need to install them manually, you can do it as follow:

```
# Packages ggplot2, rlist, stringr:
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('rlist')) install.packages('rlist'); library('rlist')
if (!require('stringr')) install.packages('stringr'); library('stringr')

# Packages pheWAS:
install.packages("devtools")
#It may be necessary to install required as not all package dependencies are installed by devtools:
install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
devtools::install_github("PheWAS/PheWAS")
library(PheWAS)

```
# Content: 

The repository contains two folders (UK Biobank genotype and phenotype data is not publicly available, so I didn't add a Data folder): 

## [a_Scripts]()

This folder contains the scripts used for association analyzis and data visualization. 

The scripts are: 



## [b_Expected_Output]()

This folder contains the different expected outputs that can be obtain after running the above scripts. 

The subfolders are:



More details about how this files are generated can be seen in the comments of each script. 

# Publication:

Part of the results from this section of my thesis are described in a publication that is being peer-reviewed. 





