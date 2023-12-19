#--------------------------------- LD Plot  -------------------------------------------------------#

# This script is to create LD plots to highlight all the SNPs that are in link of disequilibrium with
# the SNP from LGR6: rs74355478. 

# It's going to generate a plot to highlight LD SNPs. 

# ---> SET DIRECTORY: 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---> LIBRARY LOAD:

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('splitstackshape')) install.packages('splitstackshape'); library('splitstackshape')


# ---> INPUT FILES:

ld_file <- read.delim("LD_grch37_esmbl_ELGH.txt") # LD file
pop <- read.delim("Populations.txt") 
genes_location <- read.delim("../../../../Uk Biobank Project/Results_new_approach/input/7_Plots/LocusZooms-master/Gencode_GRCh37_Genes_UniqueList2021.txt", 
                             stringsAsFactors = FALSE, header = TRUE)
genes_LGR6 <- genes_location[genes_location$Gene == "LGR6", ]
genes_LGR6$Start <- genes_LGR6$Start/1000
genes_LGR6$End <- genes_LGR6$End/1000

# ---> FILE MANIPULATIONS: 

# Get the position of the SNPs in the LD_File by separating the Location column:

ld_split <- cSplit(ld_file, "LD_location",":")
ld_split <- cSplit(ld_split, "Population",":")
colnames(ld_split)[c(11,12,15)] <- c("Chr", "Position", "Population")
ld_split <- merge(ld_split, pop, all.x = T)
ld_final <- ld_split[, c(6,12,13,7,9,10,16,11)]

# We need to create a new dataframe for the Lead SNP itself:

lead_snp <- data.frame(LD_SNP = "rs74355478",
                      Chr = 1,
                      Position = 202183358/1000,
                      LD_consequence = "frameshift_variant",
                      r2 = 1,
                      D = 1,
                      Population = "ALL",
                      Ancestry = "ALL")

# ---> PLOTING ALL POPULATION TOGEGHER:

# For plotting all population together we have to delete repeated SNPs, since analysis from different
# populations will give different r2 values: 

# First we have to order the data based on the population. Thankfully, South Asian, alphabetically, can
# be the first in decreasing order. 

ld_all <- ld_final[order(ld_final$Ancestry, ld_final$r2, decreasing = TRUE), ]

# Here it will delete duplicates from the same population. 
# ld_all <- ld_all[!duplicated(ld_all[, c(1,8)]), ]

# Get the position in Kb:
ld_all$Position <- ld_all$Position/1000

# From previous deletion of duplicates, the file are order by SNPs and population, decreasing. Now we 
# want to delete duplicates considering the r2 value (highest first) but also, when the r2 values are
# equal, we have to make sure that we are choosing the one from the South Asian population. Thanks to
# the previous order you can order by SNP and r2 and have the confidence of choosing South Asian in
# tie cases. 
ld_all <- ld_all[order(ld_all$LD_SNP, ld_all$r2, decreasing = TRUE), ]

# Then you can delete duplicates based on SNP:
# ld_all <- ld_all[!duplicated(ld_all$LD_SNP), ]

# Plot the Gene Figure:

# Get x-asis limits: 
xmin <- min(ld_all$Position) - 1
xmax <- max(ld_all$Position) + 1

gene_fig <- 
  ggplot(genes_LGR6, aes(x = Start, y = 1, xend = End, yend = 1, label = Gene, colour = Gene)) +
  geom_segment(lineend = "round", linejoin = "round", size = 4, arrow = arrow(length = unit(0.3, "inches"))) +
  geom_text(x = (genes_LGR6$Start+genes_LGR6$End)/2, y = 0.9, color = "black", size = 8) +
  scale_x_continuous(name = "Position on Chromosome 1 (Kb)", limits = c(xmin, xmax)) +
  scale_color_manual(values = c("LGR6" = "royalblue2")) +
  guides(color = guide_legend(override.aes = list(color = NA))) +
  labs(color = "Populationsll") +
  ylim(c(0.7,1.1)) +
  theme(axis.title.x = element_text(size = 30, colour = "black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
        axis.text.x  =  element_text(size = 30, colour = "black"),
        axis.text.y  = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.text = element_text(size = 25, vjust = 0.5, colour = "white"), 
        legend.title = element_text(size = 30, vjust = 0.5, colour = "white"),
        legend.key = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 1, 1, 1), "cm"))

# Plot the LD Figure:

ld_fig <- 
  ggplot(ld_all, aes(x = Position, y = r2)) +   
  geom_point(colour = "gray48", size = 5) +
  geom_point(data = ld_all[ld_all$r2 >= 0.8, ], aes(color = as.factor(Ancestry)), size = 6) +
  geom_point(data = lead_snp, aes(Position, r2), size = 9, shape = 18, colour = "royalblue2") + 
  geom_line(y = 0.8, color = "black", size = 1.5, linetype = "dashed") + 
  ylab(expression(r^{2}~with~rs74355478)) +
  scale_x_continuous(name = "Position on Chromosome 1 (Kb)", limits = c(xmin, xmax),
                     breaks = seq(xmin,xmax,50000)) +
  scale_color_manual(values = c("African" = "seagreen4", "American" = "orangered3", "East Asian" = "purple3",
                                "European" = "cyan3", "South Asian" = "darkgoldenrod1")) +
  geom_text(data = lead_snp, aes(x = Position, y = 1.04, label = LD_SNP), color = "royalblue2", size = 9) +
  labs(colour="Population") +
  guides(color = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.y = element_text(size = 35, colour = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
        axis.text.x  =  element_blank(),
        axis.text.y  = element_text(size = 30, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.text = element_text(size = 25, vjust = 0.5, colour = "black"), 
        legend.title = element_text(size = 30, vjust = 0.5, colour = "black"),
        legend.key = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1), 
        axis.ticks.length.y = unit(0.4, "cm"),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 1, 0.5, 1), "cm")) 

# Creating a grid where the two plots are together and they look like one: 

pp <- list(ld_fig, gene_fig)

# Output the Final Figure: 

png(file = "../../output/2_LD_ensembl/All_population_LD_sp_0.8.png",
    width = 1200, height = 900)

plot_grid(plotlist=pp, ncol=1, align='v', axis = "l", rel_heights = c(2,0.5))

dev.off() 

# Output results table: 

write.table(ld_all[ld_all$r2 >= 0.8, ], 
            file = "../../output/2_LD_ensembl/All_populations_LD_high_0.8_sp.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE) 

# ---> PLOTING LD PER POPULATION:

# For plotting by population we have to delete repeated SNPs in the same population, since analysis 
# from different populations will give different r2 values: 

# First we have to order based on the r2 value. 

ld_pop <- ld_final[order(ld_final$r2, decreasing = TRUE), ]

# Get the bases in Kb:
ld_pop$Position <- ld_pop$Position/1000

ld_pop[(ld_pop$Position < genes_LGR6$Start | ld_pop$Position > genes_LGR6$End) 
       & ld_pop$LD_consequence == "intron_variant", ]$LD_consequence <- "intron_variant_other_gene"

# Here it will delete duplicates from the same population. 
ld_pop <- ld_pop[!duplicated(ld_pop[, c(1,7)]), ]
ld_pop$LD_consequence <- gsub("_", " ", ld_pop$LD_consequence)
ld_pop$LD_consequence <- str_to_sentence(ld_pop$LD_consequence)

# Colors for the different consequences: 
colors <- c("darkgoldenrod", "steelblue2", "violet", "olivedrab3",
            "palevioletred4", "plum2", "bisque4")
names(colors) <- c("5 prime utr variant", "Intergenic variant", "Intron variant", "Intron variant other gene",
                   "Non coding transcript exon variant", "Splice donor variant", "Splice region variant")

# Colors for the different subpopulations: 
color_pop <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "#CAB2D6", "blue1", "brown", 
  "darkorange4", "gold1", "#FDBF6F", "gray70", "khaki2","maroon", "orchid1", "deeppink1", "skyblue2", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3","palegreen2", "#FB9A99", "olivedrab3")
names(color_pop) <- unique(ld_pop$Populations)

for (i in unique(ld_pop$Ancestry)) {
  
  # Get the table for single population
  ld_table <- ld_pop[ld_pop$Ancestry == i, ]
  #ld_table_uk <- ld_table[ld_table$LD_SNP %in% c("rs61823694", "rs4266947") & ld_table$r2 >= 0.5, ]
  ld_table_uk <- ld_table[ld_table$LD_SNP %in% c("rs4266947") & ld_table$r2 >= 0.8, ]
  
  # Plot the Gene Figure:
  
  # Get the highest string for the consequence column:
  nochar <- max(nchar(ld_table[ld_table$r2 >= 0.5, ]$LD_consequence))
  
  # Create table with only values with highest string for the consequences:
  chartable <- ld_table[ld_table$r2 >= 0.5 & nchar(ld_table$LD_consequence) == nochar, ]
  
  # Create a new column with the label value:
  genes_LGR6$label <- unique(chartable$LD_consequence)
  
  # Get x-asis limits: 
  xpmin <- min(ld_table$Position) - 1
  xpmax <- max(ld_table$Position) + 1
  
  # Colors: 
  col_sel <- colors[names(colors) %in% ld_table[ld_table$r2 >= 0.5, ]$LD_consequence]
  
  gene_pop_fig <- 
    ggplot(genes_LGR6, aes(x = Start, y = 1, xend = End, yend = 1, label = Gene, colour = label)) +
    geom_segment(lineend = "round", linejoin = "round", size = 4, arrow = arrow(length = unit(0.3, "inches"))) +
    geom_text(x = (genes_LGR6$Start+genes_LGR6$End)/2, y = 0.9, color = "black", size = 8) +
    scale_x_continuous(name = "Position on Chromosome 1 (Kb)", limits = c(xpmin, xpmax)) +
    scale_color_manual(values = ("royalblue2")) +
    guides(color = guide_legend(override.aes = list(colour=NA, size=5))) +
    labs(color = "Consequence") +
    ylim(c(0.7,1.1)) +
    theme(axis.title.x = element_text(size = 30, colour = "black"),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
          axis.text.x  =  element_text(size = 30, colour = "black"),
          axis.text.y  = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.text = element_text(size = 25, vjust = 0.5, colour = "white"), 
          legend.title = element_text(size = 30, vjust = 0.5, colour = "white"),
          legend.key = element_blank(), 
          axis.line = element_line(colour = "black", size = 1),
          axis.ticks.x = element_line(colour = "black", size = 1),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0, 1, 1, 1), "cm"))
  
  # Plot the LD Figure:
  
  ld_pop_fig <- 
    ggplot(ld_table, aes(x = Position, y = r2)) +   
    geom_point(colour = "gray48", size = 5) +
    geom_point(data = ld_table[ld_table$r2 >= 0.8, ], aes(color = as.factor(LD_consequence)), size = 6) +
    geom_point(data = lead_snp, aes(Position, r2), size = 9, shape = 18, colour = "royalblue2") + 
    geom_line(y = 0.8, color = "black", size = 1.5, linetype = "dashed") + 
    ylab(expression(r^{2}~with~rs74355478)) +
    scale_x_continuous(name = "Position on Chromosome 1 (Kb)", limits = c(xmin, xmax),
                       breaks = seq(xpmin,xpmax,50000)) +
    scale_color_manual(values = col_sel) +
    geom_text(data = lead_snp, aes(x = Position, y = 1.03, label = LD_SNP), color = "royalblue2", size = 8) +
    labs(colour="Consequence") +
    guides(color = guide_legend(override.aes = list(size=7))) +
    theme(axis.title.y = element_text(size = 35, colour = "black"),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
          axis.text.x  =  element_blank(),
          axis.text.y  = element_text(size = 30, colour = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.text = element_text(size = 25, vjust = 0.5, colour = "black"), 
          legend.title = element_text(size = 30, vjust = 0.5, colour = "black"),
          legend.key = element_rect(fill = "white"),
          axis.line = element_line(colour = "black", size = 1),
          axis.ticks.y = element_line(colour = "black", size = 1), 
          axis.ticks.length.y = unit(0.4, "cm"),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1, 1, 0.5, 1), "cm")) 
  
  # Creating a grid where the two plots are together and they look like one: 
  
  ppp <- list(ld_pop_fig, gene_pop_fig)

  png(file = paste("../../output/2_LD_ensembl/", i , "_LD_sp_0.8.png", sep = ""),
      width = 1200, height = 900)
  
  print(plot_grid(plotlist=ppp, ncol=1, align='v', axis = "l", rel_heights = c(2,0.5)))
  
  dev.off() 
  
  # Plot the gene for the real populations:
  
  # Get the highest string for the consequence column:
  nochar_p <- max(nchar(ld_table[ld_table$r2 >= 0.5, ]$Populations))
  
  # Create table with only values with highest string for the consequences:
  chartable_p <- ld_table[ld_table$r2 >= 0.5 & nchar(ld_table$Populations) == nochar_p, ]
  
  # Create a new column with the label value:
  genes_LGR6$label_pop <- unique(chartable_p$Populations)
  
  gene_pops_fig <- 
    ggplot(genes_LGR6, aes(x = Start, y = 1, xend = End, yend = 1, label = Gene, colour = label_pop)) +
    geom_segment(lineend = "round", linejoin = "round", size = 4, arrow = arrow(length = unit(0.3, "inches"))) +
    geom_text(x = (genes_LGR6$Start+genes_LGR6$End)/2, y = 0.9, color = "black", size = 8) +
    scale_x_continuous(name = "Position on Chromosome 1 (Kb)", limits = c(xpmin, xpmax)) +
    scale_color_manual(values = ("royalblue2")) +
    guides(color = guide_legend(override.aes = list(colour=NA, size=5))) +
    labs(color = "Populations") +
    ylim(c(0.7,1.1)) +
    theme(axis.title.x = element_text(size = 30, colour = "black"),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
          axis.text.x  =  element_text(size = 30, colour = "black"),
          axis.text.y  = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.text = element_text(size = 25, vjust = 0.5, colour = "white"), 
          legend.title = element_text(size = 30, vjust = 0.5, colour = "white"),
          legend.key = element_blank(), 
          axis.line = element_line(colour = "black", size = 1),
          axis.ticks.x = element_line(colour = "black", size = 1),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0, 1, 1, 1), "cm"))
  
  # Plot the LD Figure for real populations:
  
  col_sel_pop <- color_pop[names(color_pop) %in% ld_table[ld_table$r2 >= 0.5, ]$Populations]
  uk_col <- color_pop[ld_table_uk$Populations]
  
  ld_pops_fig <- 
    ggplot(ld_table, aes(x = Position, y = r2)) +   
    geom_point(colour = "gray48", size = 5) +
    geom_point(data = ld_table[ld_table$r2 >= 0.8, ], aes(color = as.factor(Populations)), size = 6) +
    geom_point(data = lead_snp, aes(Position, r2), size = 10, shape = 18, colour = "royalblue2") + 
   # geom_line(y = 0.5, color = "black", size = 1.5, linetype = "dashed") + 
    geom_line(y = 0.8, color = "black", size = 1.5, linetype = "dashed") + 
    ylab(expression(r^{2}~with~rs74355478)) +
    scale_x_continuous(name = "Position on Chromosome 1 (Kb)", limits = c(xmin, xmax),
                       breaks = seq(xpmin,xpmax,50000)) +
    scale_color_manual(values = col_sel_pop) +
    geom_text(data = lead_snp, aes(x = Position, y = 1.03, label = LD_SNP), color = "royalblue2", size = 8) +
    labs(colour="Populations") +
    geom_text_repel(data = ld_table_uk,  
                    position = "identity",
                    box.padding = 4,
                    label = ld_table_uk$LD_SNP, 
                    size = 8, color = uk_col) +
    guides(color = guide_legend(override.aes = list(size=7))) +
    theme(axis.title.y = element_text(size = 35, colour = "black"),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
          axis.text.x  =  element_blank(),
          axis.text.y  = element_text(size = 30, colour = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          #legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.text = element_text(size = 25, vjust = 0.5, colour = "black"), 
          legend.title = element_text(size = 30, vjust = 0.5, colour = "black"),
          legend.key = element_rect(fill = "white"),
          axis.line = element_line(colour = "black", size = 1),
          axis.ticks.y = element_line(colour = "black", size = 1), 
          axis.ticks.length.y = unit(0.4, "cm"),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1, 1, 0.5, 1), "cm")) 
  
  # Creating a grid where the two plots are together and they look like one: 
  
  ppps <- list(ld_pops_fig, gene_pops_fig)
  
  png(file = paste("../../output/2_LD_ensembl/", i , "_LD_subpop_sp_lab_only47.png", sep = ""),
      width = 1200, height = 900)
  
  print(plot_grid(plotlist=ppps, ncol=1, align='v', axis = "l", rel_heights = c(2,0.5)))
  
  dev.off() 
  
  # Output results table: 
  
  write.table(ld_table[ld_table$r2 >= 0.6, ], 
              file = paste("../../output/2_LD_ensembl/", i , "_LD_high_sp.txt", sep = ""),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE) 
  
}

