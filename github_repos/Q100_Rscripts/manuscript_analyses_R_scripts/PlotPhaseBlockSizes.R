setwd("/Users/nhansen/HG002_diploid_benchmark/plots/benchhapmerdistances/blocksizes")

library("stringr")
readblockbedfile <- function(bedfile) {
  intervals <- read.table(bedfile, sep="\t")
  names(intervals) <- c("chrom", "start", "end", "haplotype", "score", "strand", "startthick", "endthick", "color")

  return(intervals)
}

intervals0.05_0.05 <- readblockbedfile("vary_alpha_beta/hprc_hg002_curated_0.05_0.05.v1.1.phasedscaffolds.bed")
blocklengths0.05_0.05 <- intervals0.05_0.05$end - intervals0.05_0.05$start

intervals0.01_0.05 <- readblockbedfile("vary_alpha_beta/hprc_hg002_curated_0.01_0.05.v1.1.phasedscaffolds.bed")
blocklengths0.01_0.05 <- intervals0.01_0.05$end - intervals0.01_0.05$start

histdata0.05_0.05 <- hist(blocklengths0.05_0.05, breaks=100)
histdata0.01_0.05 <- hist(blocklengths0.01_0.05, breaks=100)
plot(histdata0.05_0.05$mids, histdata0.05_0.05$counts, type='l')
points(histdata0.01_0.05$mids, histdata0.01_0.05$counts, type='l')

calc_n50 <- function(filename) {
  fileintervals <- readblockbedfile(filename)
  fileblocklengths <- fileintervals$end - fileintervals$start
  sortedintervallengths <- sort(fileblocklengths, decreasing=TRUE)
  cumsum = 0
  totalsum <- sum(sortedintervallengths)
  
  for(i in sortedintervallengths) {
    cumsum <- cumsum + i
    if (cumsum > 1500000000) {
      return(i)
    }  
  }
  totalsum <- sum(sortedintervallengths)
  
  return(sortedintervallengths[length(sortedintervallengths)])
}  

files <- list.files(path="/Users/nhansen/HG002_diploid_benchmark/plots/benchhapmerdistances/blocksizes/vary_alpha_beta", pattern="*.phasedscaffolds.bed", full.names=TRUE, recursive=FALSE)
fileparams <- sub('/Users/nhansen/HG002_diploid_benchmark/plots/benchhapmerdistances/blocksizes/vary_alpha_beta/hprc_hg002_curated_', '', files)
fileparams <- sub('.v1.1.phasedscaffolds.bed', '', fileparams)
alphavals <- sub('_.*', '', fileparams, perl=TRUE)
alphavals <- as.numeric(alphavals)
betavals <- sub('.*_', '', fileparams, perl=TRUE)
betavals <- as.numeric(betavals)

test_n50_vector <- sapply(files, function(x) {
  n50 <- calc_n50(x)
})

pointcolors <- c('red', 'blue', 'green', 'orange')
pointtype <- c(15, 16, 17, 18)
plot(alphavals, test_n50_vector/1000000, col=pointcolors[as.factor(betavals)], pch=pointtype[as.factor(betavals)], xlim=c(0, 0.3), xlab="Alternate haplotype emission probabilities", ylab="Phase block N50 (Megabases)", main="Phase block size for different HMM parameters", cex=1.1)
legend("topright", c("0.01", "0.05", "0.1", "0.2"), title="Transition probs", col=pointcolors, pch=pointtype)

