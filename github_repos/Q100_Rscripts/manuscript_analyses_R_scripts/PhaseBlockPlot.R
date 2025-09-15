setwd("/Users/nhansen/HG002_diploid_benchmark/plots/benchhapmerdistances/blocksizes")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#if (!requireNamespace("karyoploteR", quietly = TRUE))
  #BiocManager::install("karyoploteR")

library(stringr)
library(karyoploteR)
args = commandArgs(trailingOnly=TRUE)

genomename <- ifelse(!( is.na(args[1])), args[1], "hprc_hg002_curated_0.05_0.01")
outputdir <- ifelse(!( is.na(args[2])), args[2], ".")
resourcedir <- ifelse(! ( is.na(args[3])), args[3], ".")
benchname <- ifelse(!( is.na(args[4])), args[4], "v1.1")
plottitle <- ifelse(!( is.na(args[5])), args[5], paste(c(genomename, " aligned coverage vs. v1.1"), sep="", collapse=""))
genomefile <- paste(c(outputdir, "/genome.", genomename, ".bed"), sep="", collapse="")
scaffolds <- toGRanges(genomefile)
scaffolds <- sort(scaffolds)
scaffolddf <- read.table(genomefile, sep="\t", comment.char="", header=FALSE)
seqlengths(scaffolds) <- scaffolddf[,3]
bigscaffolds <- sort(seqlengths(scaffolds), decreasing=TRUE)
chroms <- names(bigscaffolds)
if (length(chroms) > 40){
  chroms <- chroms[1:40]
}

matphaseblockfile <- paste(c(outputdir, "/vary_alpha_beta/", genomename, "." , benchname, ".phasedscaffolds.mat.bed"), sep="", collapse="")
matphaseblockranges <- toGRanges(matphaseblockfile)

patphaseblockfile <- paste(c(outputdir, "/vary_alpha_beta/", genomename, "." , benchname, ".phasedscaffolds.pat.bed"), sep="", collapse="")
patphaseblockranges <- toGRanges(patphaseblockfile)

nlocfile <- paste(c(outputdir, "/nlocs.", genomename, ".bed"), sep="", collapse="")
nlocranges <- toGRanges(nlocfile)

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 100
pp$topmargin <- 300

plotname <- paste(c(outputdir, "/", genomename, ".phaseblockcovered.", benchname, ".pdf"), sep="", collapse="")
pdf(plotname, 8.5, 11.0)
#kp <- plotKaryotype(genome=scaffolds,chromosomes=paste0("chr", c(1:22, "X", "Y")), plot.type=1, main=plottitle)
kp <- plotKaryotype(genome=scaffolds, plot.type=1, chromosomes=chroms, main=plottitle, cex=0.5, plot.params=pp)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red")
if (length(nlocranges) > 0)
    kpRect(kp, data=nlocranges, col="red", data.panel="ideogram", y0=rep(0, length(nlocranges)), y1=rep(1, length(nlocranges)))
if (length(patphaseblockranges) > 0)
    kpRect(kp, data=patphaseblockranges, col="blue", data.panel="ideogram", y0=rep(0, length(patphaseblockranges)), y1=rep(1, length(patphaseblockranges)))
if (length(matphaseblockranges) > 0)
    kpRect(kp, data=matphaseblockranges, col="green", data.panel="ideogram", y0=rep(0, length(matphaseblockranges)), y1=rep(1, length(matphaseblockranges)))

legend("bottomright", c("Paternal", "Maternal", "Ns"), fill=c("blue", "green", "red"))

dev.off()

