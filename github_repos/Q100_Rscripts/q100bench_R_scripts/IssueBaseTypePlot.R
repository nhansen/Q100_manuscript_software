setwd("/Users/nhansen/HG002_diploid_benchmark/plots/readaligns")

#benchcovereddf <- read.table(benchcoveredfile, header=FALSE, sep="\t")
countfiles <- list.files(path=".", pattern=".issuebasecounts.txt")

#namestoplot <- c("Sequel_DCv1.1", "Revio_DCv1.2", "R10_Duplex", "R10_Q28", "Q28_corrected")
namestoplot <- c("Sequel_DCv1.1", "R10_Duplex", "R10_Q28", "Q28_corrected")
filenames <- paste(namestoplot, ".issuebasecounts.txt", sep="")

readissuebases <- function(filename) {
  data <- read.table(filename, header=FALSE, sep="\t")
  data$V1 <- data$V1/1000000
  return(data)
}

allissuetypes <- unique(unlist(lapply(filenames, function(x) {readissuebases(x)$V2})))
plotissuetypes <- sort(allissuetypes[allissuetypes!="Clipped"])

## NOTE: THESE LABELS SHOULD BE CHANGED/REORDERED ETC IF THE ISSUE TYPES CHANGE!!!! ##
plotissuenames <- c("High", "High AT Rich", "High GA/TC Rich", "High Poor Qual", "Low", "Low AT Rich", "Low GA/TC Rich", "Low GC Rich", "Low Poor Qual")
issuetypevalues <- function(filename, issuetypes) {
  filedf <- readissuebases(filename)
  typevalues <- sapply(issuetypes, function(x) {ifelse(length(filedf[filedf$V2==x, "V1"])>0, filedf[filedf$V2==x, "V1"],  0)})
  return(typevalues)
}

#platformcolors = c("red", "orange", "lavender", "lightblue", "darkblue")
platformcolors = c("red", "orange", "lightblue", "darkblue")
barplot(t(as.matrix(sapply(filenames, function(x){issuetypevalues(x, plotissuetypes)}))), names=plotissuenames, beside=TRUE, cex.names=0.5, cex.axis=0.7, col=platformcolors, xlab="Issue Type", ylab="Number of Megabases")
legend("topright", namestoplot, col=platformcolors, pch=15)

