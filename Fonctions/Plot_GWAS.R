library(qqman)
### Code pour tracer les résultats après eQTL avec SNPTEST ###
## Répértoire où se trouve les fichiers .out de l'eQTL
setwd("~/Etudes/M2/GENOM/Projet")

generatePlots <- function(filename,chrom) {
  # Read data from the file
  results <- read.table(filename, header=TRUE)
  
  # Add a 'chromosome' column
  results$chromosome <- chrom
  
  # Create a new dataframe for plotting
  DATA <- data.frame(SNP = results$rsid,
                     CHR = results$chromosome,
                     BP = results$position,
                     P = results$frequentist_add_pvalue)
  
  
  # Save QQ plot to a PDF file
  pdf(paste0('QQ_plot_', filename, '.pdf'))
  print(qq(DATA$P))
  dev.off()
  
  # Get SNPs having a p-value inferior to 5E-8
  threshold <-  7E-8
  DATA_significant <- subset(DATA, P < threshold)
  snpsOfInterest<-DATA_significant$SNP
  snpsOfInterest
  
  # Save Manhattan plot to a PDF file
  #pdf(paste0('Manhattan_plot_', filename, '.pdf'))
  #manhattan(subset(DATA, CHR == chrom), ylim=c(0, 20),suggestiveline=FALSE, highlight = snpsOfInterest, annotatePval = threshold, annotateTop = FALSE)
  #dev.off()
}

# Example usage
generatePlots("chr16_QTL.non_filtered.out",16)
generatePlots("chr16_QTL.filtered.out",16)
generatePlots("chr16_chol_QTL.out",16)
