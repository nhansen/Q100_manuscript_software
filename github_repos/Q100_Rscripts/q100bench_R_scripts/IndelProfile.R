#setwd("/Users/nhansen/HG002_diploid_benchmark/plots/coverage_plots")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#if (!requireNamespace("karyoploteR", quietly = TRUE))
  #BiocManager::install("karyoploteR")

library(stringr)
args = commandArgs(trailingOnly=TRUE)

genomename <- ifelse(!( is.na(args[1])), args[1], "verkko_thic_assigned_only")
benchname <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], "/Users/nhansen/tmp/svs")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c(genomename, " Indel Size Profile", genomename), sep="", collapse=""))

svboundaries <- function(svfile) {
  svboundaries <- read.table(svfile, header=FALSE, sep="\t")
  names(svboundaries) <- c("chrom", "start", "end", "svtype", "contigs", "querystart", "queryend", "refstart", "refend", "strand")
  
  return(svboundaries)
}

plotindelprofile <- function(svs) {
  samecontigsvs <- svs[svs$svtype=="SameContigInsertion" | svs$svtype=="SameContigDeletion", ]

  refdiff <- samecontigsvs$refend - samecontigsvs$refstart + 1
  querydiff <- as.integer(samecontigsvs$queryend) - as.integer(samecontigsvs$querystart) + 1
  
  plot(querydiff, refdiff, xlim=c(-500, 500), ylim=c(-500, 500))
  
  return(samecontigsvs)
}

svdefaultfile <- paste(c(outputdir, "/", genomename, ".svs.bed"), sep="", collapse="")
svs <- svboundaries(svdefaultfile)

ash1svs <- svboundaries("/Users/nhansen/tmp/svs/ash1v2.svs.bed")

#profileplotname <- paste(c(outputdir, "/", genomename, ".largeindelprofile.", benchname, ".pdf"), sep="", collapse="")
#pdf(profileplotname, 8.5, 11.0)
#plot_coverage_vs_benchmark()
#dev.off()

#errorplotname <- paste(c(outputdir, "/", genomename, ".benchcoveredwitherrors.", benchname, ".pdf"), sep="", collapse="")
#pdf(errorplotname, 8.5, 11.0)
#plot_coverage_vs_benchmark(phasingerrorgranges, snperrorgranges, indelerrorgranges)
#dev.off()

