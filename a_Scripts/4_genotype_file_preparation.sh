#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 5      # Request 5 core
#$ -l h_rt=72:0:0 # Request 20 hour runtime
#$ -l h_vmem=10G   # Request 50GB RAM
#$ -l highmem      # High Memory flag
#$ -t 1-2	# Number of lines in txt file	

# Create a list with all the populations:

pop=$(sed -n "${SGE_TASK_ID}p" populations.txt)

# Load PlINK module:

module load plink/1.9-170906

# Create the genotype file using recodeA: 
plink --recodeA --bfile ukb_cal_chr1_v2 --extract SNPs.txt --keep ${pop} --out ../../output/4_genotype_file/${pop}_genotypes