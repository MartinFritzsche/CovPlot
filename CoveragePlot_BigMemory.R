# Title: CoveragePlot_Bigmemory
# Author: Martin Fritzsche
# Date: 16/04/2018

library(bigmemory)
library(bigtabulate)
library(ggplot2)

setwd("Z:/Bioinformatic Misc/171-depth")

# Define input parameters
windowsize <- 10000
chrName <- "chrY"
project <- "171"
inputFile <- paste(chrName, ".split.depth", sep = "")

# Define binning function
slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats <- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)] 
    chunkmean <- mean(chunk)
    chunkstdv <- sd(chunk)
    chunkbps[i] <- chunkmean
    chunkstats[i] <- chunkstdv
  }
  return (list(starts,chunkbps,chunkstats))
}

# Load coverage file into big.matrix 
big_matrix <- read.big.matrix(inputFile, header = FALSE, type = "short", sep = "\t", backingfile = paste(chrName, ".bin", sep = ""), descriptorfile = paste(chrName, ".dest", sep = ""))

# Apply binning function to second column
window <- slidingwindowplot(windowsize, big_matrix[ , 2])

# Collapse into single dataframe and add chromosome name column
df <- as.data.frame(window, col.names = c("x", "mean", "sd"))
df$chr <- rep(chrName, nrow(df))
df$project <- rep(project, nrow(df))

# Assign chromosome and project-specific name to dataframe
assign(paste0(project, chrName, "_df"), df)
## comb_df <- rbind(`162chr14_df`,`162chr15_df`, `162chr16_df`, `162chr17_df`, `162chr18_df`, `162chr19_df`, `162chr20_df`, `162chr21_df`, `162chr22_df`, `162chrY_df`, `171chr14_df`, `171chr15_df`, `171chr16_df`, `171chr17_df`, `171chr18_df`, `171chr19_df`, `171chr20_df`, `171chr21_df`, `171chr22_df`, `171chrY_df`)
## comb_df <- rbind(`162chr01_df`, `162chr02_df`, `162chr03_df`, `162chr04_df`, `162chr05_df`, `162chr06_df`, `162chr07_df`, `162chr08_df`, `162chrX_df`, `171chr01_df`, `171chr02_df`, `171chr03_df`, `171chr04_df`, `171chr05_df`, `171chr06_df`, `171chr07_df`, `171chr08_df`, `171chrX_df`)
## comb_df <- rbind(`162chr09_df`, `162chr10_df`, `162chr11_df`, `162chr12_df`, `162chr13_df`, `171chr09_df`, `171chr10_df`, `171chr11_df`, `171chr12_df`, `171chr13_df`) 
## comb_df$chr <- as.factor(comb_df$chr)
## comb_df$project <- as.factor(comb_df$project)

# Generate plot
ggplot(data = df, aes(x = x, y = mean)) + 
  geom_line(colour = "#0066CC", size = 0.1) + 
  geom_ribbon(aes(ymax = mean + sd, ymin = mean - sd), alpha = 0.5, fill = "#0066CC") +
  ## facet_grid(chr ~ ., scale = "free_y", space = "free_y") +
  theme_bw() + 
  xlab("Reference Start Position") + 
  scale_x_continuous(expand = c(0,0)) +
  ylim(c(0, 50)) +
  ylab("Coverage") +
  ggtitle(paste("Coverage Across Reference", chrName))


  
# Generate comparison plot
ggplot(data = comb_df, aes(x = x, y = mean, colour = project)) + 
  geom_line(size = 0.1) + 
  facet_grid(chr ~ ., scale = "free_y", space = "free_y") +
  theme_bw() + 
  xlab("Reference Start Position") + 
  scale_x_continuous(expand = c(0,0)) +
  ylim(c(0, 50)) +
  ylab("Coverage") +
  ggtitle(paste("Coverage Across Reference", chrName))