setwd("/Users/nhansen/HG002_diploid_benchmark/plots/readtandemreps")

args = commandArgs(trailingOnly=TRUE)

readsetname <- ifelse(!( is.na(args[1])), args[1], "ont_epi2me_q28")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c("length concordance for ", readsetname, " vs ", genomename), sep="", collapse=""))

readsetlabel=readsetname

strtypes <- c("mononuc", "dinuc", "trinuc", "tetranuc")
strlabels <- c("Homopolymer", "Dinucleotide", "Trinucleotide", "Tetranucleotide")
strnumbers <- c(1, 2, 3, 4)
strminlength <- c(10, 10, 12, 12)

accbylengthdifffiles <- sapply(strtypes, function(x) {paste(c(outputdir, "/", readsetname, ".strlengthdiffaccuracy.", x, ".txt"), sep="", collapse="")})
accbylengthfiles <- sapply(strtypes, function(x) {paste(c(outputdir, "/", readsetname, ".strlengthaccuracy.", x, ".txt"), sep="", collapse="")})

lengthdifftitles <- sapply(strtypes, function(x) {paste(c(strlabels[which(strtypes==x)], " length concordance for ", readsetlabel), sep="", collapse="")})
readlengthdiffcounts <- function(file) {
  countdf <- read.table(file, header=TRUE, sep="\t")
  
  return(countdf)
}

reflengthaccuracycounts <- function(file) {
  countdf <- read.table(file, header=TRUE, sep="\t")
  
  return(countdf)
}

# lists of length 4 with df for each str type:
lengthdiffacc <- lapply(seq(length(strtypes)), function(i) {readlengthdiffcounts(accbylengthdifffiles[i])})
lengthacc <- lapply(seq(length(strtypes)), function(i) {reflengthaccuracycounts(accbylengthfiles[i])})
illuminaaccbylengthdifffiles <- sapply(strtypes, function(x) {paste(c(outputdir, "/illumina_googlepcrfree.strlengthdiffaccuracy.", x, ".txt"), sep="", collapse="")})
illuminalengthdiffacc <- lapply(seq(length(strtypes)), function(i) {readlengthdiffcounts(accbylengthdifffiles[i])})
sprq30accbylengthdifffiles <- sapply(strtypes, function(x) {paste(c(outputdir, "/hg002v1.1_revio_sprq_30hr.strlengthdiffaccuracy.", x, ".txt"), sep="", collapse="")})
sprq30lengthdiffacc <- lapply(seq(length(strtypes)), function(i) {readlengthdiffcounts(accbylengthdifffiles[i])})

plotlengthdiffhistogram <- function(df, maxdiff=10, plottitle="Length concordance accuracy rates") {
  differences <- df[df$Type=="LENGTHERROR" | df$Type=="CORRECT", "LengthDiff"]
  totalreads <- sum(df[df$Type != "HET", "NumReads"])
  accrates <- df[df$Type=="LENGTHERROR" | df$Type=="CORRECT", "NumReads"]/totalreads
  overallaccuracy <- sum(df[df$Type =="CORRECT", "NumReads"])/totalreads
  accrate <- round(100*overallaccuracy, 2)
  overallaccuracy <- sum(df[df$Type =="CORRECT", "NumReads"])/totalreads
  lengtherrorrate <- round(100*sum(df[df$Type =="LENGTHERROR", "NumReads"])/totalreads, 2)
  complexrate <- round(100*sum(df[df$Type =="COMPLEX", "NumReads"])/totalreads, 2)
  flankerrorrate <- round(100*sum(df[df$Type =="FLANKERROR", "NumReads"])/totalreads, 2)
  
  plotdiffs <- seq(-1*maxdiff, maxdiff)
  plotaccrates <- sapply(plotdiffs, function(x) { if (length(which(differences==x)) > 0) { return(accrates[which(differences==x)]) } else { return(0) }})

  barplot(plotaccrates, names=seq(-1*maxdiff,maxdiff), col="blue", cex.names=1.0, main=plottitle, xlab="Difference from benchmark", ylab="Fraction of reads", ylim=c(0,1.0))
  mtext(paste(c("Accuracy Rate: ", as.character(accrate), "%"), sep="", collapse=""), side=3, adj=0.85, line=-4)
  mtext(paste(c("Length Error Rate: ", as.character(lengtherrorrate), "%"), sep="", collapse=""), side=3, adj=0.85, line=-5)
  mtext(paste(c("Complex Error Rate: ", as.character(complexrate), "%"), sep="", collapse=""), side=3, adj=0.85, line=-6)
  mtext(paste(c("Flank Error Rate: ", as.character(flankerrorrate), "%"), sep="", collapse=""), side=3, adj=0.85, line=-7)
  
  return(plotaccrates)
}

pdf("ONTQ28DinucleotideAccuracy.pdf", 9.0, 6.0)
plotlengthdiffhistogram(lengthdiffacc[[2]], plottitle="Dinucleotide run length accuracy in ONT Q28 reads")
dev.off()
plotlengthdiffhistogram(lengthdiffacc[[4]], plottitle="Tetranucleotide run length accuracy in ONT Q28 reads")
pdf("IlluminaTetranucleotideAccuracy.pdf", 9.0, 6.0)
plotlengthdiffhistogram(illuminalengthdiffacc[[4]], plottitle="Tetranucleotide run length accuracy in Illumina PCR free reads")
dev.off()
pdf("RevioSPRQ30TrinucleotideAccuracy.pdf", 9.0, 6.0)
plotlengthdiffhistogram(sprq30lengthdiffacc[[3]], plottitle="Trinucleotide run length accuracy in Revio SPRQ 30hr reads")
dev.off()

plotstraccuracybylength <- function(df, plottitle=NA, minlength=NA, maxlength=NA, strlabel="Tandem repeat", color="black"){
  if(is.na(minlength)) {
    minlength <- min(df$RefLength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(df$RefLength)
  }
  if(is.na(plottitle)) {
    plottitle <- paste(c("Accuracy of tandem repeat runs in ", readsetname), sep="", collapse="")
  }
  accarray <- sapply(seq(minlength, maxlength), function(x) {df[df$RefLength==x, "Correct"]/(df[df$RefLength==x, "Correct"] + df[df$RefLength==x, "LengthError"] + df[df$RefLength==x, "ComplexError"] + df[df$RefLength==x, "FlankError"])})
  #totalaccrate <- as.integer(sum(mncounts$numcorrect)*1000/sum(mncounts$numerror+mncounts$numcorrect))/10
  plot(seq(minlength, maxlength), accarray, pch=16, xlab=c(paste(c(strlabel, " run length"), sep="", collapse="")), ylab=c("Accuracy"), main=plottitle, ylim=c(0,1), col=color)
  #if (totalaccrate < 75) {
    #textheight <- 0.9
    #textpos <- minlength + 3*(maxlength-minlength)/4    
  #}
  #else {
    ##textheight <- accarray[as.integer((maxlength-minlength)/3)]/2
    #textheight <- 0.2
    #textpos <- minlength + (maxlength-minlength)/4
  #}
  #text(textpos, textheight, labels= paste(c("Overall: ", totalaccrate, "%"), sep="", collapse=""))
}

multiplotstraccuracybylength <- function(accdfs, plottitle=NA, minlength=NA, maxlength=NA, strlabel="Tandem repeat", pchvals=NA, labels=NA, colors=c("black")){
  df <- accdfs[[1]]
  if(is.na(minlength)) {
    minlength <- min(df$RefLength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(df$RefLength)
  }
  if(is.na(plottitle)) {
    plottitle <- paste(c("Accuracy of tandem repeat runs in ", readsetname), sep="", collapse="")
  }
  accarray <- sapply(seq(minlength, maxlength), function(x) {df[df$RefLength==x, "Correct"]/(df[df$RefLength==x, "Correct"] + df[df$RefLength==x, "LengthError"] + df[df$RefLength==x, "ComplexError"] + df[df$RefLength==x, "FlankError"])})
  plot(seq(minlength, maxlength), accarray, pch=pchvals[1], xlab=c(paste(c(strlabel, " run length"), sep="", collapse="")), ylab=c("Accuracy"), main=plottitle, ylim=c(0,1), col=colors[1])

  if (length(accdfs)>=2) {
    for (i in seq(2, length(accdfs))) {
      df <- accdfs[[i]]
      accarray <- sapply(seq(minlength, maxlength), function(x) {df[df$RefLength==x, "Correct"]/(df[df$RefLength==x, "Correct"] + df[df$RefLength==x, "LengthError"] + df[df$RefLength==x, "ComplexError"] + df[df$RefLength==x, "FlankError"])})
      points(seq(minlength, maxlength), accarray, pch=pchvals[i], col=colors[i])
    }
  }
  legend("bottomleft", labels, pch=pchvals, col=colors)
}

platformcolors <- c( "#117733", "#CC6677", "#DDCC77", "#332288", "#88CCEE")
platformnames <- c("ONT Q28", "HiFi SPRQ/30 hrs", "Illumina GooglePCRFree", "NIST Onso", "Element Ultraq")
platformpchvals <- c(0, 1, 2, 5, 6, 8)
filledplatformpchvals <- c(15, 16, 17, 18, 20, 8)

multiplotdatasets <- c("ont_epi2me_q28", "hg002v1.1_revio_sprq_30hr", "illumina_googlepcrfree", "NIST_Onso_2024Q1", "element_ultraq_jun2024")
multiplotmononucfiles <- lapply(multiplotdatasets, function(x) {paste(c(x, ".strlengthaccuracy.mononuc.txt"), sep="", collapse="")})
mononucaccdfs <- lapply(multiplotmononucfiles, function(x) {reflengthaccuracycounts(x)})
multiplotdinucfiles <- lapply(multiplotdatasets, function(x) {paste(c(x, ".strlengthaccuracy.dinuc.txt"), sep="", collapse="")})
dinucaccdfs <- lapply(multiplotdinucfiles, function(x) {reflengthaccuracycounts(x)})
multiplottrinucfiles <- lapply(multiplotdatasets, function(x) {paste(c(x, ".strlengthaccuracy.trinuc.txt"), sep="", collapse="")})
trinucaccdfs <- lapply(multiplottrinucfiles, function(x) {reflengthaccuracycounts(x)})
multiplottetranucfiles <- lapply(multiplotdatasets, function(x) {paste(c(x, ".strlengthaccuracy.tetranuc.txt"), sep="", collapse="")})
tetranucaccdfs <- lapply(multiplottetranucfiles, function(x) {reflengthaccuracycounts(x)})

pdf("HomopolymerAccMultiplatform.pdf")
multiplotstraccuracybylength(mononucaccdfs, colors=platformcolors, labels=platformnames, pchvals=filledplatformpchvals, maxlength=35, plottitle="Homopolymer accuracy by reference length", strlabel=strlabels[1])
dev.off()
pdf("DinucleotideAccMultiplatform.pdf")
multiplotstraccuracybylength(dinucaccdfs, colors=platformcolors, labels=platformnames, pchvals=filledplatformpchvals, maxlength=50, plottitle="Dinucleotide accuracy by reference length", strlabel=strlabels[2])
dev.off()
pdf("TrinucleotideAccMultiplatform.pdf")
multiplotstraccuracybylength(trinucaccdfs, colors=platformcolors, labels=platformnames, pchvals=filledplatformpchvals, maxlength=50, plottitle="Trinucleotide accuracy by reference length", strlabel=strlabels[3])
dev.off()
pdf("TetranucleotideAccMultiplatform.pdf")
multiplotstraccuracybylength(tetranucaccdfs, colors=platformcolors, labels=platformnames, pchvals=filledplatformpchvals, maxlength=50, plottitle="Tetranucleotide accuracy by reference length", strlabel=strlabels[4])
dev.off()

onebaseindelratio <- function(filename) {
  onebaseratio <- mndiffs(filename, 1)*1.0/mndiffs(filename, -1)
  
  return(onebaseratio)
}

makemultiplot <- function(maxdiff=4) {
  plotcolors <- c("blue", "lightblue", "lightgreen", "darkgreen", "violet", "darkviolet")
  par(oma=c(5,5,5,5))
  par(mfrow = c(2, 3))
  barplot(diffarray("sequelII_15k.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[1], cex.names=1.0, main="SequelII_15k", xlab="Difference from benchmark", ylab="Fraction of reads")
  barplot(diffarray("sequelII_20kb_20210812.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[2], cex.names=1.0, main="SequelII_20kb_20210812", xlab="Difference from benchmark")
  barplot(diffarray("hifidcv1.1.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[3], cex.names=1.0, main="HiFi/DCv1.1", xlab="Difference from benchmark")
  barplot(diffarray("hifirevio3cells.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[4], cex.names=1.0, main="HiFi Revio", xlab="Difference from benchmark")
  barplot(diffarray("illumina2x250mat.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[5], cex.names=0.7, main="Illumina 2x250 (maternal)", xlab="Difference from benchmark")
  barplot(diffarray("element_avitilongmat.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[6], cex.names=0.7, main="Element Aviti (maternal)", xlab="Difference from benchmark")
  mtext("Homopolymer length concordance", side = 3, line = 1, outer = TRUE)
}

plotmononuchist <- function(file, plottitle="", maxdiff=4){
  differences <- diffarray(file, maxdiff=maxdiff)
  accrate <- as.integer(differences[maxdiff+1]*1000)/10
  barplot(differences, names=seq(-1*maxdiff,maxdiff), col="blue", cex.names=1.0, main=plottitle, xlab="Difference from benchmark", ylab="Fraction of reads", ylim=c(0,1.0))
  mtext(paste(c("Accuracy: ", as.character(accrate), "%"), sep="", collapse=""), side=3, adj=0.85, line=-4)
  indelratio <- as.integer(onebaseindelratio(file)*1000)/1000
  mtext(paste(c("1bp Ins/Dels: ", as.character(indelratio)), sep="", collapse=""), side=3, adj=0.85, line=-5 )
}

plotmononucaccuracy <- function(file, plottitle=NA, minlength=NA, maxlength=NA, color="black"){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- paste(c("Accuracy of mononucleotide runs in ", genomename), sep="", collapse="")
  }
  accarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numcorrect"])/sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})
  #  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  totalaccrate <- as.integer(sum(mncounts$numcorrect)*1000/sum(mncounts$numerror+mncounts$numcorrect))/10
  plot(seq(minlength, maxlength), accarray, pch=16, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), main=plottitle, ylim=c(0,1), col=color)
  if (totalaccrate < 75) {
    textheight <- 0.9
    textpos <- minlength + 3*(maxlength-minlength)/4    
  }
  else {
    #textheight <- accarray[as.integer((maxlength-minlength)/3)]/2
    textheight <- 0.2
    textpos <- minlength + (maxlength-minlength)/4
  }
  text(textpos, textheight, labels= paste(c("Overall: ", totalaccrate, "%"), sep="", collapse=""))
}

par(mfrow=c(1,1))
mononucfile <- paste(c(outputdir, "/", readsetname, ".mononuchist.txt"), sep="", collapse="")

plotname <- paste(c(outputdir, "/", readsetname, ".mononucsizehist.pdf"), sep="", collapse="")
pdf(plotname, 8.5, 5.5)
plotmononuchist(mononucfile, plottitle, 4)
dev.off()

plotname <- paste(c(outputdir, "/", readsetname, ".mononucaccuracy.pdf"), sep="", collapse="")
pdf(plotname, 8.5, 5.5)
plotmononucaccuracy(mononucfile, plottitle, maxlength=40)
dev.off()

comparison_plot_may_2024 <- function(maxdiff=4) {
  plotcolors <- palette.colors(n=6)
  par(oma=c(5,5,5,5))
  par(mfrow = c(2, 3))
  barplot(diffarray("element_avitilongmat.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[1], cex.names=0.7, main="Element Aviti (maternal)", xlab="Difference from benchmark", ylim=c(0,1))
  barplot(diffarray("illumina2x250mat_nonexc.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[2], cex.names=0.7, main="Illumina 2x250 (maternal)", xlab="Difference from benchmark", ylim=c(0,1))
  barplot(diffarray("hifi_revio_pbmay24.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[3], cex.names=1.0, main="HiFi Revio/Late 2023", xlab="Difference from benchmark", ylab="Fraction of reads", ylim=c(0,1))
  barplot(diffarray("ONT_Q28.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[4], cex.names=1.0, main="ONT Q28/Nov 2023", xlab="Difference from benchmark", ylim=c(0,1))
  barplot(diffarray("ONT_R10_duplex_nonexc.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[5], cex.names=1.0, main="ONT R10 Duplex", xlab="Difference from benchmark", ylim=c(0,1))
  barplot(diffarray("ONT_Q28_herrocorrected.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[6], cex.names=1.0, main="Herro-corrected ONT", xlab="Difference from benchmark", ylim=c(0,1))
  mtext("Homopolymer length concordance", side = 3, line = 1, outer = TRUE)
}
accuracy_comparison_plot_may_2024 <- function() {
  plotcolors <- palette.colors(n=6)
  par(oma=c(5,5,5,5))
  par(mfrow = c(2, 3))
  plotmononucaccuracy("element_avitilongmat.mononuchist.txt", maxlength=40, color=plotcolors[1], plottitle="Element Aviti (maternal)")
  plotmononucaccuracy("illumina2x250mat_nonexc.mononuchist.txt", maxlength=40, color=plotcolors[2], plottitle="Illumina 2x250 (maternal)")
  plotmononucaccuracy("hifi_revio_pbmay24.mononuchist.txt", maxlength=40, color=plotcolors[3], plottitle="HiFi Revio/Late 2023")
  plotmononucaccuracy("ONT_Q28.mononuchist.txt", maxlength=40, color=plotcolors[4], plottitle="ONT Q28/Nov 2023")
  plotmononucaccuracy("ONT_R10_duplex_nonexc.mononuchist.txt", maxlength=40, color=plotcolors[5], plottitle="ONT R10 Duplex")
  plotmononucaccuracy("ONT_Q28_herrocorrected.mononuchist.txt", maxlength=40, color=plotcolors[6], plottitle="Herro-corrected ONT")
  mtext("Homopolymer length accuracy", side = 3, line = 1, outer = TRUE)
}

