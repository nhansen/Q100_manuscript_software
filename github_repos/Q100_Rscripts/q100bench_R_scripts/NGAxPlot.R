setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/plots/ngaplots")

args = commandArgs(trailingOnly=TRUE)
assemblycolors <- c("#44AA99", "#332288", "#882255", "#888888", "#DDCC77", "#661100", "#6699CC")
methodcolors <- c("#DDCC77", "#661100", "#6699CC")

assemblyname <- ifelse(!( is.na(args[1])), args[1], "lc24_ul_herro_ulk")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
idealfile <- ifelse(!( is.na(args[4])), args[4], "v1.1.non_n_seq.bed")
plottitle <- ifelse(!( is.na(args[5])), args[5], paste(c("NGAx ", assemblyname, " vs ", genomename), sep="", collapse=""))
lengthtypes <- c("alignclusterlengths", "contiglengths", "scaffoldlengths")
assemblysizefiles <- sapply(lengthtypes, function(x) {file=paste(c(outputdir, "/", assemblyname, ".", x, ".txt"), sep="", collapse=""); return(file)})

readlengths <- function(sizefile) {
  clusterlengths <- read.table(sizefile, sep="\t", header=FALSE)
  names(clusterlengths) <- c("perc", "clusterlength", "totallength", "chromosome")
  clusterlengths$clusterlength <- clusterlengths$clusterlength/1000000
  
  return(clusterlengths)
}

readideallengths <- function(bedfile) {
  idealdf <- read.table(bedfile)
  idealdf$lengths <- (as.integer(idealdf$V3)-as.integer(idealdf$V2))/1000000
  
  totallengths <- sum(idealdf$lengths)
  sortedlengths <- sort(idealdf$lengths, decreasing=TRUE)
  perccovered <- sapply(seq(1, length(sortedlengths)), function(i) {return(100*sum(sortedlengths[seq(1,i)])/totallengths)})
  
  return(data.frame(clusterlength=sortedlengths, perc=perccovered))
}

plotlengths <- function(clusterlengths, color="red", title="NGAx", dashed=FALSE, ltyval=NA, cexval=1.0) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  if (is.na(ltyval)) {
    ltyval <- ifelse(dashed, 2, 1)
  }
  plot(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,100), ylim=c(0,250), xlab="Percent of Diploid Genome", ylab="Aligned length", main=title, cex=cexval)
}

addclusterlengths <- function(clusterlengths, color="blue", dashed=FALSE, ltyval=NA) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  if (is.na(ltyval)) {
    ltyval <- ifelse(dashed, 2, 1)
  }
  lines(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,100))
}


assembly_ngax_plot <- function(clusterfiles, contigfiles=c(), scaffoldfiles=c(), assemblylabels=c(), ideal=FALSE, haplotype="MAT", plottitle="", cexval=1.0) {
  
  firstclusters <- readlengths(clusterfiles[1]) 
  plotlengths(firstclusters, col=assemblycolors[1], title=plottitle, lty=1, cexval=cexval)
  if (length(contigfiles) > 1) {
    firstcontigs <- readlengths(contigfiles[1]) 
    addclusterlengths(firstcontigs, col=assemblycolors[1])    
  }
  if (length(scaffoldfiles) > 1) {
    firstscaffolds <- readlengths(scaffoldfiles[1]) 
    addclusterlengths(firstscaffolds, col=assemblycolors[1])    
  }
  
  if (length(clusterfiles) > 1) {
    for (i in seq(2, length(clusterfiles))) {
      #addclusterlengths(readlengths(clusterfiles[i]), col=assemblycolors[i], lty=i)
      addclusterlengths(readlengths(clusterfiles[i]), col=assemblycolors[i], lty=1)
    }
  }
  if (ideal) {
    ideallengths <- readideallengths(idealfile)
    addclusterlengths(ideallengths, col="black")
    assemblycolors <- c(assemblycolors[1:length(assemblylabels)], "black")
    assemblylabels <- c(assemblylabels, "Ideal (HG002v1.1)")
  }
  legend("topright", assemblylabels, col=assemblycolors, lty=rep(1, length(clusterfiles)))
  
}

# Make NGAx plot:
plotname <- paste(c(outputdir, "/", assemblyname, ".continuitystats.pdf"), sep="", collapse="")
pdf(plotname)
assembly_ngax_plot(assemblysizefiles, assemblylabels=c("Aligned NGAx", "Contig NGx", "Scaffold NGx"), ideal=TRUE, haplotype=NA, plottitle=plottitle)
dev.off()

