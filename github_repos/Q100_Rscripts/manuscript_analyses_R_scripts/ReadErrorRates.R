setwd("/Users/nhansen/HG002_diploid_benchmark/plots/readaligns")

singlebaseerrorfiles <- c("elementins1000mat.singlebaseerrorrates.txt", "hifirevio3cell.singlebaseerrorrates.txt", "ont_q28_corrected.singlebaseerrorrates.txt", "hifidcv1.1.singlebaseerrorrates.txt", "ont_epi2me_q28.singlebaseerrorrates.txt")

readerrorcounts <- function(filename) {
  ecdf <- read.table(filename, header=FALSE, sep="\t")
  names(ecdf) <- c("rate", "refbase", "readbase")
  ecdf$errorspermillion <- 1000000*ecdf$rate
   
  return(ecdf)
}

readmononuchist <- function(filename) {
  mndf <- read.table(filename, header=TRUE, sep="")
  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  mndf <- mndf[mndf$readlength != -1, ]
  mndf$reflength <- as.integer(mndf$reflength)
  mndf$readlength <- as.integer(mndf$readlength)
  mndf$numcorrect <- as.integer(mndf$numcorrect)
  mndf$numerror <- as.integer(mndf$numerror)
  
  return(mndf)
}

readmononuccounts <- function(filename) {
  mndf <- read.table(filename, header=FALSE, sep="")
  names(mndf) <- c("count", "base", "reflength", "readlength", "type")
  mndf <- mndf[mndf$readlength != "COMPLEX", ]
  
  return(mndf)
}

platformnames <- c("Element", "Illumina", "Sequel", "R10 Q28", "R10 Q28 Corrected", "R10 Duplex")
plotcolors <- c("green", "violet", "red", "lightblue", "blue", "orange")
pointtype <- c(15, 20, 16, 17, 23, 19)
plot(seq(10,25), accuracylist(countfiles[1]), ylim=c(0,1.0), col=plotcolors[1], pch=pointtype[1], main="Accuracy of Mononucleotide Runs >= 10 bp", xlab="Length of Run", ylab="Fraction Correct")
for (i in seq(2,6)) {
  points(seq(10,25), accuracylist(countfiles[i]), col=plotcolors[i], pch=pointtype[i])
}
legend("bottomleft", platformnames, col=plotcolors, pch=pointtype)

mndiffs <- function(filename, i) {
  mndiff <- readmononuchist(filename);

  if (i==0) {
    diffcount <- sum(mndiff[mndiff$readlength==mndiff$reflength, "numcorrect"])
  }
  else {
    diffcount <- sum(mndiff[mndiff$readlength-mndiff$reflength==i, "numerror"])
  }
  #diffcount <- sum(mndf[mndf$type!="HET" & mndf$type !="COMPLEX" & as.integer(mndf$readlength)-as.integer(mndf$reflength)==i, "count"])

  return(diffcount)
}

diffarray <- function(filename) {
  platformdiffs <- sapply(seq(-4, 4), function(i) {mndiffs(filename, i)})
  platformdiffs <- platformdiffs/sum(platformdiffs)
  
  return(platformdiffs)
}
plotcolors <- c("blue", "lightblue", "lightgreen", "darkgreen", "violet", "darkviolet")
countfiles <- c("hifirevio3cell.mononucsummary.txt", "elementins1000mat.mononucsummary.txt", "illumina2x250.mononucsummary.txt" )
par(oma=c(5,5,5,5))
par(mfrow = c(2, 3))
barplot(diffarray("sequelII_15k.mononuchist.txt"), names=seq(-4,4), col=plotcolors[1], cex.names=1.0, main="SequelII_15k", xlab="Difference from benchmark", ylab="Fraction of reads")
barplot(diffarray("sequelII_20kb_20210812.mononuchist.txt"), names=seq(-4,4), col=plotcolors[2], cex.names=1.0, main="SequelII_20kb_20210812", xlab="Difference from benchmark")
barplot(diffarray("hifidcv1.1.mononuchist.txt"), names=seq(-4,4), col=plotcolors[3], cex.names=1.0, main="HiFi/DCv1.1", xlab="Difference from benchmark")
barplot(diffarray("illumina2x250mat.mononuchist.txt"), names=seq(-4,4), col=plotcolors[5], cex.names=0.7, main="Illumina 2x250 (maternal)", xlab="Difference from benchmark")
barplot(diffarray("element_avitilongmat.mononuchist.txt"), names=seq(-4,4), col=plotcolors[6], cex.names=0.7, main="Element Aviti (maternal)", xlab="Difference from benchmark")
mtext("Homopolymer length concordance", side = 3, line = 1, outer = TRUE)

par(mfrow=c(1,1))
barplot(diffarray("hifirevio3cell.mononucsummary.txt"), names=seq(-10,10), col=plotcolors[1], cex.names=0.7, main=paste(c("Revio (", as.integer(diffarray("hifirevio3cell.mononucsummary.txt")[11]*100 + 0.5), "% Accurate)"), sep="", collapse=""))
barplot(diffarray("elementins1000mat.mononucsummary.txt"), names=seq(-10,10), col=plotcolors[1], cex.names=0.7, main=paste(c("Element (", as.integer(diffarray("elementins1000mat.mononucsummary.txt")[11]*100 + 0.5), "% Accurate)"), sep="", collapse=""), ylim=c(0,0.05))
barplot(diffarray("illumina2x250.mononucsummary.txt"), names=seq(-10,10), col="violet", cex.names=0.7, main=paste(c("Illumina (", as.integer(diffarray("illumina2x250.mononucsummary.txt")[11]*100 + 0.5), "% Accurate)"), sep="", collapse=""), ylim=c(0,0.05))
