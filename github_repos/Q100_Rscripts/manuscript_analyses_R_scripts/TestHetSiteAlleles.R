setwd("/Users/nhansen/HG002_diploid_benchmark/plots/hetsite_plots")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#if (!requireNamespace("karyoploteR", quietly = TRUE))
  #BiocManager::install("karyoploteR")

library(stringr)
library(karyoploteR)
args = commandArgs(trailingOnly=TRUE)

genomename <- ifelse(!( is.na(args[1])), args[1], "v2_trio.hap1")
outputdir <- ifelse(!( is.na(args[2])), args[2], ".")
resourcedir <- ifelse(! ( is.na(args[3])), args[3], ".")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c(genomename, " alleles at v1.0.1 het sites"), sep="", collapse=""))
genomefile <- paste(c(outputdir, "/genome.", genomename, ".bed"), sep="", collapse="")
scaffolds <- toGRanges(genomefile)
scaffolds <- sort(scaffolds)
scaffolddf <- read.table(genomefile, sep="\t", comment.char="", header=FALSE)
seqlengths(scaffolds) <- scaffolddf[,3]
bigscaffolds <- sort(seqlengths(scaffolds), decreasing=TRUE)
chroms <- names(bigscaffolds)
if (length(chroms) > 25){
  chroms <- chroms[1:25]
}

hetallelefile <- paste(c(outputdir, "/", "v1.0.1", ".coveredhetalleles.", genomename, ".bed"), sep="", collapse="")
matcoveredfile <- paste(c(outputdir, "/testmatcovered.", genomename, ".merged.bed"), sep="", collapse="")
patcoveredfile <- paste(c(outputdir, "/testpatcovered.", genomename, ".merged.bed"), sep="", collapse="")
matcoveredranges <- toGRanges(matcoveredfile)
patcoveredranges <- toGRanges(patcoveredfile)

nlocfile <- paste(c(outputdir, "/nlocs.", genomename, ".bed"), sep="", collapse="")
nlocranges <- toGRanges(nlocfile)

hetalleles <- read.table(hetallelefile, sep="\t")
names(hetalleles) <- c("chrom", "start", "end", "name", "allele", "alignedchrom", "alignedstart", "alignedend", "contig", "contigstart", "contigend", "samediff")

hetmatallelessamehap <- hetalleles[str_detect(hetalleles$alignedchrom, "MATERNAL") & hetalleles$samediff=="SAMEHAP",]
hetmatallelesalthap <- hetalleles[str_detect(hetalleles$alignedchrom, "MATERNAL") & hetalleles$samediff=="ALTHAP",]

hetpatallelessamehap <- hetalleles[str_detect(hetalleles$alignedchrom, "PATERNAL") & hetalleles$samediff=="SAMEHAP",]
hetpatallelesalthap <- hetalleles[str_detect(hetalleles$alignedchrom, "PATERNAL") & hetalleles$samediff=="ALTHAP",]

hetmatallelessamehap$end <- hetmatallelessamehap$start + 1
hetmatallelesalthap$end <- hetmatallelesalthap$start + 1

hetpatallelessamehap$end <- hetpatallelessamehap$start + 1
hetpatallelesalthap$end <- hetpatallelesalthap$start + 1

hetmatranges <- as(rbind(hetmatallelessamehap, hetpatallelesalthap),  "GRanges")
hetpatranges <- as(rbind(hetpatallelessamehap, hetmatallelesalthap),  "GRanges")
#hetmatranges <- as(rbind(hetmatallelessamehap, hetmatallelesalthap),  "GRanges")
#hetpatranges <- as(rbind(hetpatallelessamehap, hetpatallelesalthap),  "GRanges")

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 100
pp$data1height <- 1000
#pp$data2height <- 1000
pp$topmargin <- 300

plotname <- paste(c(outputdir, "/", genomename, ".testcovered.v1.0.1.pdf"), sep="", collapse="")
#pdf(plotname, 8.5, 11.0)

kp <- plotKaryotype(genome=scaffolds, plot.type=2, chromosomes=chroms, main=plottitle, cex=0.5, plot.params=pp)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=0.4)

hetpatplusmat <- as(rbind(hetmatallelessamehap, hetpatallelesalthap, hetpatallelessamehap, hetmatallelesalthap), "GRanges")
kpPlotDensity(kp, data=hetpatplusmat, data.panel=1, r0=0.05, r1=0.9, col="blue")
kpPlotDensity(kp, data=hetmatranges, data.panel=1, r0=0.05, r1=0.9, col="green")

if (length(nlocranges) > 0)
    kpRect(kp, data=nlocranges, col="red", data.panel="ideogram", y0=rep(0, length(nlocranges)), y1=rep(1, length(nlocranges)))

legend("bottomright", c("Paternal", "Maternal", "Ns"), fill=c("blue", "green", "red"))

dev.off()
