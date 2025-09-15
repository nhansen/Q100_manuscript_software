setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/PaperFigures/Figures/VCFgenome")

library(colorspace)
library(Hmisc)

# no scientific notation on axes
options(scipen=5)

qualcolors <- c("#44AA99", "#332288", "#882255", "#888888")
barcolors <- sapply(qualcolors, function(x) {c(x, x)})

errorstats <- read.table("totalhcbasesandvariantcounts.v3issues.txt", sep="\t")
names(errorstats) <- c("ref", "mingq", "covered", "mergedcovered", "totalerrors", "snperrors", "indelerrors")

errorstats$percentcovered <- errorstats$covered/5999408148
errorstats$totalqv <- -10.0*log10(errorstats$totalerrors/errorstats$covered)
errorstats$snpqv <- -10.0*log10(errorstats$snperrors/errorstats$covered)
errorstats$indelqv <- -10.0*log10(errorstats$indelerrors/errorstats$covered)

grch38errorstats <- errorstats[errorstats$ref=="GRCh38",]
chm13errorstats <- errorstats[errorstats$ref=="CHM13",]

pdf("PBRevioDVBasesCovered.pdf")
q10covered <- c(errorstats[errorstats$mingq==10 & errorstats$ref=="CHM13", "covered"]/1000000, errorstats[errorstats$mingq==10 & errorstats$ref=="GRCh38", "covered"]/1000000) 
q20covered <- c(errorstats[errorstats$mingq==20 & errorstats$ref=="CHM13", "covered"]/1000000, errorstats[errorstats$mingq==20 & errorstats$ref=="GRCh38", "covered"]/1000000) 
q30covered <- c(errorstats[errorstats$mingq==30 & errorstats$ref=="CHM13", "covered"]/1000000, errorstats[errorstats$mingq==30 & errorstats$ref=="GRCh38", "covered"]/1000000) 
q40covered <- c(errorstats[errorstats$mingq==40 & errorstats$ref=="CHM13", "covered"]/1000000, errorstats[errorstats$mingq==40 & errorstats$ref=="GRCh38", "covered"]/1000000) 
barplot(rbind(q10covered, q20covered, q30covered, q40covered), ylab="Mb covered", ylim=c(0,10000), names.arg=c("CHM13", "GRCh38"), beside=TRUE, col=t(barcolors), main="Total bases covered by variant-constructed genomes")
legend("topright", c("10", "20", "30", "40"), col=qualcolors, pch=15, title="Min GQ", cex=1.4)
dev.off()

pdf("PBRevioDVCGTotalErrors.pdf")
q10errors <- c(errorstats[errorstats$mingq==10 & errorstats$ref=="CHM13", "totalerrors"], errorstats[errorstats$mingq==10 & errorstats$ref=="GRCh38", "totalerrors"]) 
q20errors <- c(errorstats[errorstats$mingq==20 & errorstats$ref=="CHM13", "totalerrors"], errorstats[errorstats$mingq==20 & errorstats$ref=="GRCh38", "totalerrors"]) 
q30errors <- c(errorstats[errorstats$mingq==30 & errorstats$ref=="CHM13", "totalerrors"], errorstats[errorstats$mingq==30 & errorstats$ref=="GRCh38", "totalerrors"]) 
q40errors <- c(errorstats[errorstats$mingq==40 & errorstats$ref=="CHM13", "totalerrors"], errorstats[errorstats$mingq==40 & errorstats$ref=="GRCh38", "totalerrors"]) 
barplot(rbind(q10errors, q20errors, q30errors, q40errors), ylab="Errors", ylim=c(0,1800000), names.arg=c("CHM13", "GRCh38"), beside=TRUE, col=t(barcolors), main="Total errors in variant-constructed genomes")
legend("topright", c("10", "20", "30", "40"), col=qualcolors, pch=15, title="Min GQ", cex=1.4)
dev.off()

pdf("PBRevioDVCGSubsRates.pdf")
q10subsrates <- c(errorstats[errorstats$mingq==10 & errorstats$ref=="CHM13", "snperrors"]/(errorstats[errorstats$mingq==10 & errorstats$ref=="CHM13", "covered"]/1000000), errorstats[errorstats$mingq==10 & errorstats$ref=="GRCh38", "snperrors"]/(errorstats[errorstats$mingq==10 & errorstats$ref=="GRCh38", "covered"]/1000000)) 
q20subsrates <- c(errorstats[errorstats$mingq==20 & errorstats$ref=="CHM13", "snperrors"]/(errorstats[errorstats$mingq==20 & errorstats$ref=="CHM13", "covered"]/1000000), errorstats[errorstats$mingq==20 & errorstats$ref=="GRCh38", "snperrors"]/(errorstats[errorstats$mingq==20 & errorstats$ref=="GRCh38", "covered"]/1000000)) 
q30subsrates <- c(errorstats[errorstats$mingq==30 & errorstats$ref=="CHM13", "snperrors"]/(errorstats[errorstats$mingq==30 & errorstats$ref=="CHM13", "covered"]/1000000), errorstats[errorstats$mingq==30 & errorstats$ref=="GRCh38", "snperrors"]/(errorstats[errorstats$mingq==30 & errorstats$ref=="GRCh38", "covered"]/1000000)) 
q40subsrates <- c(errorstats[errorstats$mingq==40 & errorstats$ref=="CHM13", "snperrors"]/(errorstats[errorstats$mingq==40 & errorstats$ref=="CHM13", "covered"]/1000000), errorstats[errorstats$mingq==40 & errorstats$ref=="GRCh38", "snperrors"]/(errorstats[errorstats$mingq==40 & errorstats$ref=="GRCh38", "covered"]/1000000)) 
barplot(rbind(q10subsrates, q20subsrates, q30subsrates, q40subsrates), ylab="Substitution errors per Mb", ylim=c(0,350), names.arg=c("CHM13", "GRCh38"), beside=TRUE, col=t(barcolors), main="Substitution error rates in variant-constructed genomes")
legend("topright", c("10", "20", "30", "40"), col=qualcolors, pch=15, title="Min GQ", cex=1.4)
dev.off()

pdf("PBRevioDVCGIndelRates.pdf")
q10indelrates <- c(errorstats[errorstats$mingq==10 & errorstats$ref=="CHM13", "indelerrors"]/(errorstats[errorstats$mingq==10 & errorstats$ref=="CHM13", "covered"]/1000000), errorstats[errorstats$mingq==10 & errorstats$ref=="GRCh38", "indelerrors"]/(errorstats[errorstats$mingq==10 & errorstats$ref=="GRCh38", "covered"]/1000000)) 
q20indelrates <- c(errorstats[errorstats$mingq==20 & errorstats$ref=="CHM13", "indelerrors"]/(errorstats[errorstats$mingq==20 & errorstats$ref=="CHM13", "covered"]/1000000), errorstats[errorstats$mingq==20 & errorstats$ref=="GRCh38", "indelerrors"]/(errorstats[errorstats$mingq==20 & errorstats$ref=="GRCh38", "covered"]/1000000)) 
q30indelrates <- c(errorstats[errorstats$mingq==30 & errorstats$ref=="CHM13", "indelerrors"]/(errorstats[errorstats$mingq==30 & errorstats$ref=="CHM13", "covered"]/1000000), errorstats[errorstats$mingq==30 & errorstats$ref=="GRCh38", "indelerrors"]/(errorstats[errorstats$mingq==30 & errorstats$ref=="GRCh38", "covered"]/1000000)) 
q40indelrates <- c(errorstats[errorstats$mingq==40 & errorstats$ref=="CHM13", "indelerrors"]/(errorstats[errorstats$mingq==40 & errorstats$ref=="CHM13", "covered"]/1000000), errorstats[errorstats$mingq==40 & errorstats$ref=="GRCh38", "indelerrors"]/(errorstats[errorstats$mingq==40 & errorstats$ref=="GRCh38", "covered"]/1000000)) 
barplot(rbind(q10indelrates, q20indelrates, q30indelrates, q40indelrates), ylab="Insertion/deletion errors per Mb", ylim=c(0,150), names.arg=c("CHM13", "GRCh38"), beside=TRUE, col=t(barcolors), main="Indel error rates in variant-constructed genomes")
legend("topright", c("10", "20", "30", "40"), col=qualcolors, pch=15, title="Min GQ", cex=1.4)
dev.off()

pdf("PBRevioDVCGQVscores.pdf")
q10qvscores <- c(-10.0*log10(errorstats[errorstats$mingq==10 & errorstats$ref=="CHM13", "totalerrors"]/(errorstats[errorstats$mingq==10 & errorstats$ref=="CHM13", "covered"])), -10.0*log10(errorstats[errorstats$mingq==10 & errorstats$ref=="GRCh38", "totalerrors"]/(errorstats[errorstats$mingq==10 & errorstats$ref=="GRCh38", "covered"]))) 
q20qvscores <- c(-10.0*log10(errorstats[errorstats$mingq==20 & errorstats$ref=="CHM13", "totalerrors"]/(errorstats[errorstats$mingq==20 & errorstats$ref=="CHM13", "covered"])), -10.0*log10(errorstats[errorstats$mingq==20 & errorstats$ref=="GRCh38", "totalerrors"]/(errorstats[errorstats$mingq==20 & errorstats$ref=="GRCh38", "covered"]))) 
q30qvscores <- c(-10.0*log10(errorstats[errorstats$mingq==30 & errorstats$ref=="CHM13", "totalerrors"]/(errorstats[errorstats$mingq==30 & errorstats$ref=="CHM13", "covered"])), -10.0*log10(errorstats[errorstats$mingq==30 & errorstats$ref=="GRCh38", "totalerrors"]/(errorstats[errorstats$mingq==30 & errorstats$ref=="GRCh38", "covered"]))) 
q40qvscores <- c(-10.0*log10(errorstats[errorstats$mingq==40 & errorstats$ref=="CHM13", "totalerrors"]/(errorstats[errorstats$mingq==40 & errorstats$ref=="CHM13", "covered"])), -10.0*log10(errorstats[errorstats$mingq==40 & errorstats$ref=="GRCh38", "totalerrors"]/(errorstats[errorstats$mingq==40 & errorstats$ref=="GRCh38", "covered"]))) 
barplot(rbind(q10qvscores, q20qvscores, q30qvscores, q40qvscores), ylab="GQC QV score", ylim=c(0,45), names.arg=c("CHM13", "GRCh38"), beside=TRUE, col=t(barcolors), main="GQC QV scores of variant-constructed genomes")
legend("bottomright", c("10", "20", "30", "40"), col=qualcolors, pch=15, title="Min GQ", cex=1.4, bg="white")
dev.off()

## file of counts considering only covered bases/errors within the intersected covered regions for GRCh38-constructed and CHM13-constructed genomes:
errorstats <- read.table("totalintersectedcoverageanderrors.txt", sep="\t")
names(errorstats) <- c("mingq", "intersectcovered", "totalgrch38errors", "totalchm13errors")
errorstats$percentcovered <- errorstats$intersectcovered/5999408148

#> q10qvscores
#[1] 37.21716 38.26332
#> q20qvscores
#[1] 38.47365 39.28307
#> q30qvscores
#[1] 39.40000 40.07201
#> q40qvscores
#[1] 40.52527 41.01121

q10qvscores <- c(-10.0*log10(errorstats[errorstats$mingq==10, "totalgrch38errors"]/(errorstats[errorstats$mingq==10, "intersectcovered"])), -10.0*log10(errorstats[errorstats$mingq==10, "totalchm13errors"]/(errorstats[errorstats$mingq==10, "intersectcovered"]))) 
q20qvscores <- c(-10.0*log10(errorstats[errorstats$mingq==20, "totalgrch38errors"]/(errorstats[errorstats$mingq==20, "intersectcovered"])), -10.0*log10(errorstats[errorstats$mingq==20, "totalchm13errors"]/(errorstats[errorstats$mingq==20, "intersectcovered"])))
q30qvscores <- c(-10.0*log10(errorstats[errorstats$mingq==30, "totalgrch38errors"]/(errorstats[errorstats$mingq==30, "intersectcovered"])), -10.0*log10(errorstats[errorstats$mingq==30, "totalchm13errors"]/(errorstats[errorstats$mingq==30, "intersectcovered"])))
q40qvscores <- c(-10.0*log10(errorstats[errorstats$mingq==40, "totalgrch38errors"]/(errorstats[errorstats$mingq==40, "intersectcovered"])), -10.0*log10(errorstats[errorstats$mingq==40, "totalchm13errors"]/(errorstats[errorstats$mingq==40, "intersectcovered"])))

