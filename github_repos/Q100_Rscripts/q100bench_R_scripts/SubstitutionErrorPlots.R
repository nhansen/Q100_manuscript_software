setwd("/Users/nhansen/HG002_diploid_benchmark/plots/readerrors")

library(stringr)
args = commandArgs(trailingOnly=TRUE)

readsetname <- ifelse(!( is.na(args[1])), args[1], "elementins1000mat")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.0.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste("Substitution Error Rates for ", readsetname, sep="", collapse=""))

# Plot substitution rate-by-type histogram

typeorder <- c("A_C", "A_G", "A_T", "T_C", "T_G", "T_A", "G_A", "G_T", "G_C", "C_A", "C_T", "C_G")
titv <- c("tv", "ti", "tv", "ti", "tv", "tv", "ti", "tv", "tv", "tv", "ti", "tv" )
plottypesubsrates <- function(readsetname, outputdir, xlabval="Substitution Type", ylabval="Errors per mb", titleval="", ymax=NA, legend=FALSE) {
  subsfile <- paste(c(outputdir, "/", readsetname, ".singlenucerrorstats.txt"), sep="", collapse="")
  
  if (titleval != "") {
    barplottitle <- titleval
  }
  else {
    barplottitle <- plottitle
  }
  subshist <- read.table(subsfile, sep="\t")
  names(subshist) <- c("errortype", "errorcount", "errorspermbaligned")

  #typecolors <- c(1,1,1) * palette.colors(palette = "Alphabet", n = 4)
  colorpalette <- palette.colors(palette="Alphabet", n=12)

  twocolors <- c("#0000FF", "#FF0000")
  
  typecolors <- as.character(sapply(seq(1,12), function(i) {ifelse(titv[i]=="tv", twocolors[1], twocolors[2])}))
  #typecolors <- as.character(colorpalette)
  indexorder <- sapply(typeorder, function(x) {which(subshist$errortype==x)})
  valueorder <- subshist[indexorder,]$errorspermbaligned
  if (is.na(ymax)) {
    barplot(valueorder, names.arg=typeorder, col=typecolors, main=barplottitle, xlab=xlabval, ylab=ylabval)
  }
  else {
    barplot(valueorder, names.arg=typeorder, col=typecolors, main=barplottitle, xlab=xlabval, ylab=ylabval, ylim=c(0,ymax))
  }
  if (legend) {
    legend("topright", c("Transversion","Transition"), col=as.character(c(twocolors[1], twocolors[2])), pch=15)
  }
}

setwd("/Users/nhansen/HG002_diploid_benchmark/plots/elementcomparison")
plottypesubsrates("element_ins500_600mat/element_ins500_600mat", outputdir, titleval="Element 500-600bp ins substitution rates", legend=TRUE)

