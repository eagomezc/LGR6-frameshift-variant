#---------------------------------------- PHEWAS --------------------------------------------------#

# UK BIobank is a repository with ~ 500000 samples (genomic and general data from healthy and diseases 
# participants of all kind). 

# PheWAS analysis takes SNPs and associated them with different traits or diseases.

# This scripts takes the PheWAS results and plot them in a more presentable way. 

# ---> SET DIRECTORY: 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---> LIBRARY LOAD:
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('stringr')) install.packages('stringr'); library('stringr')

# Color squeme for variants: 
trait_color <- function(traits) {
  colors <- c("blue", "cyan3", "brown", "chocolate2",
              "mediumorchid2", "dodgerblue4", "darkolivegreen4", "firebrick2", "chocolate4", 
              "mediumseagreen", "black", "midnightblue", "violetred1", "forestgreen", 
              "aquamarine", "darkorchid", "lightskyblue4")  
  names(colors) <- c("infectious diseases", "neoplasms", "endocrine/metabolic", "hematopoietic",
                     "mental disorders", "neurological", "sense organs", "circulatory system",      
                     "respiratory", "digestive", "genitourinary", "pregnancy complications", "dermatologic",
                     "musculoskeletal", "congenital anomalies", "symptoms", "injuries & poisonings")
  col_sel <- colors[names(colors) %in% traits]
  col_sel
}

# ---> INPUT FILES:

# PheWAS files: 

filenames <- list.files("../6_phewas figures/", pattern="_pca.tab", full.names=TRUE) 
ldf <- lapply(filenames, read.delim) # Open the files and put them in a list. 
names(ldf) <- substr(filenames, 21,100) # Get the names of the tables 
names(ldf) <- gsub("_T_pca.tab", "", names(ldf)) # Delete innecesary strings in filenames

# Create loop to do the analysis for all the tables:

for (i in 1:length(ldf)) {
  
  # Open the phewas table: 
  phewas_table <- ldf[[i]]
  
  # Filter results:
  phewas_table <- phewas_table[!is.na(phewas_table$p), ] # Delete not enough samples row
  phewas_table$adjustment <- NULL # Delete empty column
  
  # Print number of phenotypes:
  no_phe <- paste(names(ldf)[[i]], " phenotypes = ", nrow(phewas_table), sep = "")
  print(no_phe)
  
  # Calculate FDR:
  phewas_table$BH <- p.adjust(phewas_table$p, method = "BH")

  # Create Manhattan Plot: 
  
  # Get columns I care
  table_plot <- phewas_table[, c(2,3,6,9, 18,19)]
  table_plot$snp <- gsub("_T", "", table_plot$snp) # Delete consequence in snp name
  table_plot$group <- factor(table_plot$group, levels = unique(table_plot$group)) # Get order of traits
  table_plot <- table_plot[order(table_plot$group), ]
  table_plot$index <- seq(1, nrow(table_plot)) # Add index value for every data point
  table_plot$shape <- "positive"
  table_plot[table_plot$beta < 0, ]$shape <- "negative"
  
  # Indicates the center for each trait name:
  
  axis_set <- table_plot %>% 
    group_by(group) %>%  
    summarize(center = mean(index))
  
  # Get Uppercase for first letter in group:
  axis_set$group <- str_to_sentence(axis_set$group)
  
  # Established the colors for the different traits:
  color_traits <- trait_color(unique(table_plot$group))
  
  # Get table with p value order:
  table_plot_or <- table_plot[order(table_plot$p), ]
  
  # Create the plot:
  
  phewas_plot <- 
  ggplot(table_plot, aes(x=index, y=-log10(p))) +
    geom_point(aes(color=group, shape=shape), size=5) + # All points
    scale_color_manual(values = color_traits) +
    scale_shape_manual(values = c("negative"="\u25BC", "positive"="\u25B2")) +
    # Location of the labels for every gene:
    scale_x_continuous(label = axis_set$group, breaks= axis_set$center) +
    # Horizontal line representing significant SNPs:
    geom_hline(yintercept = -log10(0.05/nrow(table_plot)), color="red", linetype="dashed", linewidth = 1) +
    geom_hline(yintercept = -log10(0.05), color="mediumblue", linetype="dashed", linewidth = 1) +
    # Add trait name to significant SNPs: 
    geom_text_repel(data = table_plot_or[table_plot_or$p <= (0.05/nrow(table_plot_or)), ], 
                    position = "identity",
                    label = table_plot_or[table_plot_or$p <= (0.05/nrow(table_plot_or)), ]$description, 
                    size = 8, color="red") +
    geom_text_repel(data = head(table_plot_or[table_plot_or$p <= 0.05 & table_plot_or$p > (0.05/nrow(table_plot_or)), ],  n = 10), 
                    position = "identity", 
                    label = head(table_plot_or[table_plot_or$p <= 0.05 & table_plot_or$p > (0.05/nrow(table_plot_or)), ],  n = 10)$description, 
                    size = 8, color="mediumblue", max.overlaps = 6) +
    # x and y labels:
    xlab("Phenotypes") + 
    ylab(expression(-log[10](p))) +
    ggtitle(unique(table_plot$snp)) + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 25),
      axis.title = element_text(size = 30, colour = "black"),
      axis.text.x  =  element_text(size = 20, vjust = 1, hjust=1, angle =45, colour = "black"), 
      axis.text.y  = element_text(size = 25, hjust = 1, colour = "black"),
      axis.line = element_line(colour = 'black', linewidth = 0.5), 
      axis.ticks = element_line(colour = "black", linewidth = 0.5),  
      axis.ticks.length = unit(0.1, "cm"),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.caption = element_text(color = "Black", size = 30))
  
  png(file = paste("../../output/6_phewas figures/pca/", names(ldf)[[i]], ".png", sep = ""), 
      width = 1200, height = 900)
  
  print(phewas_plot)
  
  dev.off()
  
}
  