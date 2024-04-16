#Script for plotting histograms of allele frequency data created using Tuomas Hamala's poly_freq.c code where the name of the file is written in a .txt file called variables.txt
#Tuomas Hämälä, April 2023

library(ggplot2)
library(gridExtra)

# Read the names in the variables file 
filenames <-readLines("variables.txt")

# Select only filenames ending in.txt
txt_files <- filenames[grep("\\.txt$", filenames)]

#Create a folder to store plots
dir.create("AF_plots")

# Read population names from file
pop_names <- readLines("pops.txt")

# Define the number of populations
num_populations <- length(pop_names)

# Read in the AF data
for (i in seq_along(pop_names)) {
  population <- pop_names[i]
  file_path <- paste0(population, "_AF.txt")
  population_df <- read.delim(file_path, sep = "\t")
  
# Generate a color for each population using hsv
  color <- hsv((i - 1) / num_populations, 1, 1)
    
# Create the ggplot object using the population data
  Allele_freq_plot <- ggplot(population_df, aes(x = .data[[population]], fill = factor(1))) +
    geom_histogram(bins = 100, fill = color) +
    xlab("Allele Frequencies") +
    ggtitle(population) +
    theme_classic() 
  
  # Save each plot as a separate PDF file
   ggsave(filename = paste0("AF_plots/", population, "_allele_freq_plot.png"), plot = Allele_freq_plot, width = 3, height = 2)
}


