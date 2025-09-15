setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/AllCorrectionAnalysis")

read_corrections <- function(file) {
  correctiondf <- read.table(file, header=FALSE, sep=" ");
  names(correctiondf) <- c("refallele", "altallele");
  correctiondf$corrlength <- nchar(correctiondf$altallele)-nchar(correctiondf$refallele);
  
  return(correctiondf)
}

corrlength_count <- function(correctiondf, corrlength) {
  numcorrs <- length(correctiondf[correctiondf$corrlength==corrlength, "corrlength"])
  return(numcorrs)
}

phase1corrections <- read_corrections("hg002v0.7_to_hg002v0.9.allelecorrections.txt")
phase2corrections <- read_corrections("hg002v0.9_to_hg002v1.0.allelecorrections.txt")
phase3corrections <- read_corrections("hg002v1.0.1_to_hg002v1.1.allelecorrections.txt")

barplotlengths <- seq(-6, 6)
phase1counts <- sapply(barplotlengths, function(count) {return(corrlength_count(phase1corrections, count))})
phase2counts <- sapply(barplotlengths, function(count) {return(corrlength_count(phase2corrections, count))})
phase3counts <- sapply(barplotlengths, function(count) {return(corrlength_count(phase3corrections, count))})

pdf("CorrectionLengthCountsThreeRounds.pdf")
barplot(rbind(phase1counts, phase2counts, phase3counts), col=c("brown", "blue", "darkgreen"), names.arg=barplotlengths, main="Corrections lengths in three phases of polishing", xlab="Correction length", ylab="Number of corrections")
legend("topright", c("Phase 1", "Phase 2", "Phase 3"), col=c("brown", "blue", "darkgreen"), bty='n', pch=15)
dev.off()

# compare with net sum of all corrections (i.e., difference lengths between v0.7 and v1.1):
allphasecorrections <- read_corrections("netcorrectionv0.7tov1.1.allelecorrections.txt")
allphasecounts <- sapply(barplotlengths, function(count) {return(corrlength_count(allphasecorrections, count))})

barplot(rbind(phase1counts+phase2counts+phase3counts, allphasecounts), col=c("black", "darkgray"), beside=TRUE, names.arg=barplotlengths, main="Total corrections lengths all phases vs. net corrections made start to end", xlab="Correction length", ylab="Number of corrections")
legend("topright", c("Total polishing changes", "Net changes"), col=c("black", "darkgray"), pch=15)

# Compare the numbers:
df <- data.frame("changesize"=barplotlengths, "totalchanges"=phase1counts+phase2counts+phase3counts, "netchanges"=allphasecounts)
df$fractionpreserved <- df$netchanges/df$totalchanges

