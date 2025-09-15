setwd("/Users/nhansen/OneDrive/HPRC_assembly_comparison/data_for_Rplots/continuity")

numclustersbychrom <- read.table("samples116.numclustersbychrom.txt", sep="\t")
names(numclustersbychrom) <- c("chrom", "numclusters", "numsamples")

collapsedbasecounts <- read.table("samples116.multicovered.txt", sep=" ")
names(collapsedbasecounts) <- c("sample", "collapsedbases")

uncoveredbasecounts <- read.table("samples116.uncovered.txt", sep=" ")
names(uncoveredbasecounts) <- c("sample", "uncoveredbases")

inversionchromosomes <- read.table("samples116.inversionchroms.txt")
inversionchromhist <- read.table("inversionchromhist.txt", sep="\t")
candidateinversions <- read.table("samples116.inversioncandidates.txt", sep="\t", comment=":")
names(candidateinversions) <- 

chimericchromhist <- read.table("chimeric_chroms.hist.txt", sep="\t")

svnumber <- read.table("samples116.numsvs.txt", sep="\t")
names(svnumber) <- c("SVNumber", "sample")

discrepancynumber <- read.table("samples116.discrepancycounts.txt", sep="\t")
names(discrepancynumber) <- c("sample", "discrepancynumber")

safe_colorblind_palette <- c("#332288", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
pardefault <- par()

plotnumclusterbars <- function(numclustersbychrom) {
  allchroms <- unique(numclustersbychrom$chrom)

  barcolors <- safe_colorblind_palette[0:6]
    
  onepiece <- sapply(allchroms, function(x) { if (length(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters==1, "numsamples"]) > 0) {return(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters==1, "numsamples"])} else {return(0)} })
  twopiece <- sapply(allchroms, function(x) { if (length(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters==2, "numsamples"]) > 0) {return(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters==2, "numsamples"])} else {return(0)} })
  threepiece <- sapply(allchroms, function(x) { if (length(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters==3, "numsamples"]) > 0) {return(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters==3, "numsamples"])} else {return(0)} })
  fourpiece <- sapply(allchroms, function(x) { if (length(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters==4, "numsamples"]) > 0) {return(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters==4, "numsamples"])} else {return(0)} })
  fiveormorepiece <- sapply(allchroms, function(x) { if (length(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters>=5, "numsamples"]) > 0) {return(sum(numclustersbychrom[numclustersbychrom$chrom==x & numclustersbychrom$numclusters>=5, "numsamples"]))} else {return(0)} })
  
  barplot(rbind(onepiece, twopiece, threepiece, fourpiece, fiveormorepiece), names=allchroms, beside=FALSE, col=barcolors)

  return(cbind(onepiece, twopiece, threepiece, fourpiece, fiveormorepiece))
}

pdf("EstimatedCollapsedBasesSamples116.pdf")
hist(collapsedbasecounts$collapsedbases/1000000, breaks=seq(0, 50, 2), xlab="Estimated collapsed bases (mb)", ylab ="Number of HPRC release 2 samples", main="Megabases of collapsed sequence", col="#CC6677")
dev.off()

pdf("EstimatedUncoveredBasesSamples116.pdf")
hist(uncoveredbasecounts$uncoveredbases/1000000, breaks=seq(0, 110, 5), xlab="Estimated uncovered bases (mb)", ylab ="Number of HPRC release 2 samples", main="Megabases of uncovered sequence", col="#6699CC")
dev.off()

pdf("InversionDiscrepanciesSamples116.pdf", height=8.0, width=18.0)
barplot(inversionchromhist$V2/2, names=inversionchromhist$V1, main="113 inversion discrepancies in 116 samples", ylab="Number of samples", col="#44AA99")
dev.off()

pdf("InversionSizesSamples116.pdf", height=8.0, width=8.0)
hist(candidateinversions[seq(1, length(candidateinversions$V1), by=2), "V1"], breaks=seq(0, 1300000, by=50000), xlab="Inversion size (with repeat)", ylab="Number of inversions", col="#999933", main="Sizes of 113 inversion discrepancies")
dev.off()

par(mar = c(8.1, 4.1, 4.1, 2.1))
pdf("ChimericChromosomeAlignments.pdf", height=8.0, width=18.0)
barplot(chimericchromhist$V1, names=chimericchromhist$V2, las=2, ylab="Number of scaffolds", cex.names=0.8, col="#332288", main="Chimeric chromosome alignments by chromosome combination")
dev.off()

par(pardefault)

pdf("SVNumberSamples116.pdf", height=8.0, width=18.0)
orderedsvnumber <- svnumber[order(svnumber$SVNumber, decreasing=TRUE), ]
barplot(orderedsvnumber$SVNumber, names.arg=orderedsvnumber$sample, las=2, cex.names=0.4, col="#88CCEE", main="Number of alignment breaks in verkko2/hifiasm comparisons")
dev.off()

pdf("DiscrepancyNumberSamples116.pdf", height=8.0, width=18.0)
ordereddiscnumber <- discrepancynumber[order(discrepancynumber$discrepancynumber, decreasing=TRUE), ]
barplot(ordereddiscnumber$discrepancynumber, names.arg=ordereddiscnumber$sample, las=2, cex.names=0.4, col="#661100", main="Number of discrepancies within verkko2/hifiasm alignments")
dev.off()
