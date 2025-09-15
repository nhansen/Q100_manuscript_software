setwd("/Users/nhansen/HG002_diploid_benchmark/plots")
library("stringr")

alignments <- read.table("v1.0.1.mataut_vs_pataut.1to1alignable.bed", sep="\t")
names(alignments) <-c("chrom1", "start1", "end1", "coords2", "score", "strand")

chromlengths <- read.table("v1.0.1.genome.bed", sep="\t", header=FALSE)
names(chromlengths) <- c("chrom", "start", "length")

alignments$chrom2 = str_extract(alignments$coords2, regex("[^:]+(?=:)"))
alignments$start2 = ifelse(alignments$strand=="F", as.integer(str_extract(alignments$coords2, regex("\\d+(?=-\\d+)"))), as.integer(str_extract(alignments$coords2, regex("\\d+(?=/[FR]$)"))))
alignments$end2 = ifelse(alignments$strand=="F", as.integer(str_extract(alignments$coords2, regex("\\d+(?=/[FR]$)"))), as.integer(str_extract(alignments$coords2, regex("\\d+(?=-\\d+)"))))


plot_alignsegments <- function(aligns, chrlengths, chrom1, chrom2) {
  chrom1length <- as.integer(chrlengths[chrlengths$chrom == chrom1, "length"]/1000000)
  chrom2length <- as.integer(chrlengths[chrlengths$chrom == chrom2, "length"]/1000000)
  plot(1, type="l", xlab=chrom1, ylab=chrom2, xlim=c(0,chrom1length), ylim=c(0, chrom2length))
  chromaligns <- aligns[aligns$chrom1==chrom1 & aligns$chrom2==chrom2, ]
  segments(as.integer(chromaligns$start1/1000000), as.integer(chromaligns$start2/1000000), as.integer(chromaligns$end1/1000000), as.integer(chromaligns$end2/1000000), col=ifelse(chromaligns$strand=="F", "blue", "red"))
}
plot(1, type="l", xlab="chrom1", ylab="chrom2", xlim=c(0, 10), ylim=c(0, 10))
