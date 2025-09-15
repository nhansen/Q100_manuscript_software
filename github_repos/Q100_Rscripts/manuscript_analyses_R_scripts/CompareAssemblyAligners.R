setwd("/Users/nhansen/HG002_diploid_benchmark/PaperFigures")

library(colorspace)
library(Hmisc)

################
### Supplemental Figure ?? ###
### Adam's description: ###
### Historical assemblies evaluated against the HG002 genome benchmark. ###
### [Comparison of historical assemblies against the benchmark to include ###
### things like Ash1, HPRCv1, new Verkko, new Hifiasm, and whatever older ###
### HG002 assemblies we can find in GenBank to show progress towards ###
### personalized genomes.] ###
################

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
assemblycolors <- c("#44AA99", "#332288", "#882255", "#888888", "#AA4499", "#88CCEE", "#661100", "#999933")
methodcolors <- c("#DDCC77", "#661100", "#6699CC")
assemblynames <- c("ash1v2.mm2def.splits", "ash1v2.wm2def.splits", "hprc_year1_mat.mm2def.splits", "hprc_year1_mat.wm2def.splits", "hprc_year2_hap1.mm2def.splits", "hprc_year2_hap1.wm2def.splits", "hifi_q28_trio_hic_hap1.mm2def.splits", "hifi_q28_trio_hic_hap1.wm2def.splits")
assemblylabels <- c("Ash1v2 2020/MM", "Ash1v2 2020/WM", "Yr1 HPRC/MM", "Yr1 HPRC/WM", "Yr2 HPRC/MM", "Yr2 HPRC/WM", "HiFi Q28/MM", "HiFi Q28/WM")
#assemblylabels <- c("Shumate et al, 2020", "hifiasm, HPRC 2021", "verkko2.2, 80x Q28 ONT 2024", "verkko2.2, 60x Revio 2024")
# NGA plot
assemblysizefiles <- sapply(assemblynames, function(x) {file=paste(c("CompareAlignerOutput/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})

readlengths <- function(sizefile) {
  clusterlengths <- read.table(sizefile, sep="\t", header=FALSE)
  names(clusterlengths) <- c("perc", "clusterlength", "totallength", "chromosome")
  clusterlengths$clusterlength <- clusterlengths$clusterlength/1000000
  
  return(clusterlengths)
}

plotclusterlengths <- function(clusterlengths, color="red", title="NGAx", dashed=FALSE, ltyval=NA, cexval=1.0) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  if (is.na(ltyval)) {
    ltyval <- ifelse(dashed, 2, 1)
  }
  plot(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,100), ylim=c(0,250), xlab="Percent of Haploid Genome", ylab="Aligned length", main=title, cex=cexval)
}

addclusterlengths <- function(clusterlengths, color="blue", dashed=FALSE, ltyval=NA) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  if (is.na(ltyval)) {
    ltyval <- ifelse(dashed, 2, 1)
  }
  lines(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,100))
}

plotsingleassemblyresults <- function(clusterfile, contigfile, scaffoldfile, ideal=TRUE, haplotype="MAT", plottitle="") {
  cluster_lengths <- readlengths(clusterfile)
  contig_lengths <- readlengths(contigfile)
  scaffold_lengths <- readlengths(scaffoldfile)
  
  plotclusterlengths(cluster_lengths, col="red", title=plottitle)
  addclusterlengths(contig_lengths, col="blue")
  addclusterlengths(scaffold_lengths, col="darkgreen")
  if (ideal) {
    addclusterlengths(ideallengths, col="black")
  }
  
  if (ideal) {
    legend("topright", c("ideal", "NGx (scaffolds)", "NGx (contigs)", "NGAx"), col=c("black", "darkgreen", "blue", "red"), pch=15)
  }
  else {
    legend("topright", c("NGx (scaffolds)", "NGx (contigs)", "NGAx"), col=c("darkgreen", "blue", "red"), pch=15)
  }
}

assembly_ngax_plot <- function(clusterfiles, contigfiles=c(), scaffoldfiles=c(), assemblylabels=c(), ideal=FALSE, haplotype="MAT", plottitle="", cexval=1.0) {

  firstclusters <- readlengths(clusterfiles[1]) 
  plotclusterlengths(firstclusters, col=assemblycolors[1], title=plottitle, lty=1, cexval=cexval)
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
    ideallengths <- readlengths("CompareAlignerOutput/v1.1.mat.alignclusterlengths.txt")
    addclusterlengths(ideallengths, col="black")
  }
  #legend("topright", assemblylabels, col=assemblycolors, lty=seq(1, length(clusterfiles)))
  legend("topright", assemblylabels, col=assemblycolors, lty=rep(1, length(clusterfiles)))
  
}

# Make NGAx plot:

pdf("CompareAlignerOutput/NGAxPlot.pdf")
assembly_ngax_plot(assemblysizefiles, assemblylabels=assemblylabels, ideal=TRUE, plottitle="Assembly NGAx for minimap2, winnowmap2")
dev.off()

# Mononucleotide accuracy:
mnstatsfiles <- sapply(assemblynames, function(x) {file=paste(c("CompareAlignerOutput/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})

readmnstatsfile <- function(filename) {
  mnstats <- read.table(filename, header=FALSE, sep="\t")
  names(mnstats) <- c("name", "base", "reflength", "assemblylength", "type")
 
  return(mnstats) 
}

mononucaccuracystats <- function(mnstats) {
  consensuserrors <- mnstats[(mnstats$assemblylength != -1) & (mnstats$type == "CONSENSUS"), ]
  noncomplexcovered <- mnstats[mnstats$assemblylength != -1, ]
  consensuserrorcounts <- hist(consensuserrors$reflength, plot=FALSE, breaks=seq(10, 100, 1))
  noncomplexcovcounts <- hist(noncomplexcovered$reflength, plot=FALSE, breaks=seq(10, 100, 1))
  accrate <- 1.0 - consensuserrorcounts$counts/noncomplexcovcounts$counts
  correctcalls <- noncomplexcovcounts$counts - consensuserrorcounts$counts

  return(data.frame(lengths=consensuserrorcounts$mids, errors=consensuserrorcounts$counts, correctcalls=correctcalls, totals=noncomplexcovcounts$counts, accuracy=accrate) ) 
}

assembly_mononucacc_plot <- function(mnstatsfiles=c(), assemblylabels=c(), plottitle="", pointcex=1.0, errorbars=FALSE) {
  pchvalues <- c(15, 16, 17, 18, 19, 20, 21, 22)
  multiplevector <- c(1.4, 1.4, 1.4, 1.9, 1.4, 1.4, 1.4, 1.4)
  ptexp <- pointcex*multiplevector

  if (length(mnstatsfiles) > 0) {
    firststats <- readmnstatsfile(mnstatsfiles[1])
    genomename <- assemblylabels[1]

    plotdata <- mononucaccuracystats(firststats)    
    mnlengths <- plotdata[seq(1,30), "lengths"]
    accuracies <- plotdata[seq(1,30), "accuracy"]
    accconfints <- binconf(plotdata[seq(1,30), "correctcalls"], plotdata[seq(1,30), "totals"], return.df=TRUE)
    lowaccs <- accconfints$Lower
    highaccs <- accconfints$Upper
    plot(mnlengths, accuracies, xlim=c(10,40), ylim=c(0.0,1.0), pch=pchvalues[1], cex=ptexp[1], col=assemblycolors[1], xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), main=c("Accuracy of mononucleotide runs"))
    #text(20, 0.5, labels= paste(c("Overall error rate: ", mononucerrorperc, "%"), sep="", collapse=""))
    if (errorbars) {
      arrows(x0=mnlengths, y0=lowaccs, x1=mnlengths, y1=highaccs, code=3, angle=90, length=0.1, col=assemblycolors[1])
    }
  }
  if (length(mnstatsfiles) > 1) {
      for (i in seq(2, length(mnstatsfiles))) {
        nextstats <- readmnstatsfile(mnstatsfiles[i])
        plotdata <- mononucaccuracystats(nextstats)
        mnlengths <- plotdata[seq(1,30), "lengths"]
        accuracies <- plotdata[seq(1,30), "accuracy"]
        accconfints <- binconf(plotdata[seq(1,30), "correctcalls"], plotdata[seq(1,30), "totals"], return.df=TRUE)
        lowaccs <- accconfints$Lower
        highaccs <- accconfints$Upper
        points(mnlengths, accuracies, xlim=c(10,40), cex=ptexp[i], pch=pchvalues[i], col=assemblycolors[i])
        if (errorbars) {
          arrows(x0=mnlengths, y0=lowaccs, x1=mnlengths, y1=highaccs, code=3, angle=90, length=0.1, col=assemblycolors[i])
        }
      }    
  }

  legend(15, 0.6, assemblylabels, col=assemblycolors, pch=pchvalues, pt.cex=ptexp)
}

assembly_mononucqv_plot <- function(mnstatsfiles=c(), assemblylabels=c(), plottitle="Accuracy of mononucleotide runs", errorbars=FALSE, pointcex=1.0) {
  pchvalues <- c(15, 16, 17, 18, 19, 20, 21, 22)
  multiplevector <- c(1.4, 1.4, 1.4, 1.9, 1.4, 1.4, 1.4, 1.4)
  ptexp <- pointcex*multiplevector
  
  if (length(mnstatsfiles) > 0) {
    firststats <- readmnstatsfile(mnstatsfiles[1])
    genomename <- assemblylabels[1]
    plotdata <- mononucaccuracystats(firststats)    
    mnlengths <- plotdata[seq(1,30), "lengths"]
    errorrates <- plotdata[seq(1,30), "errors"]/plotdata[seq(1,30), "totals"]
    errorconfints <- binconf(plotdata[seq(1,30), "errors"], plotdata[seq(1,30), "totals"], return.df=TRUE)
    lowqvals <- -10.0*log10(errorconfints$Lower)
    highqvals <- -10.0*log10(errorconfints$Upper)
    
    qvals <- sapply(errorrates, function(x) { qv=-10.0*log10(x); return(qv) })
    plot(mnlengths, qvals, xlim=c(10,40), ylim=c(0, 40), pch=pchvalues[1], col=assemblycolors[1], xlab=c("Mononucleotide run length"), ylab=c("Phred QV Score"), main=plottitle, cex=ptexp[1])
    if (errorbars) {
      arrows(x0=mnlengths, y0=lowqvals, x1=mnlengths, y1=highqvals, code=3, angle=90, length=0.1, col=assemblycolors[1])
    }
  }
  if (length(mnstatsfiles) > 1) {
    for (i in seq(2, length(mnstatsfiles))) {
      nextstats <- readmnstatsfile(mnstatsfiles[i])
      plotdata <- mononucaccuracystats(nextstats)
      mnlengths <- plotdata[seq(1,30), "lengths"]
      errorrates <- plotdata[seq(1,30), "errors"]/plotdata[seq(1,30), "totals"]
      errorconfints <- binconf(plotdata[seq(1,30), "errors"], plotdata[seq(1,30), "totals"], return.df=TRUE)
      lowqvals <- -10.0*log10(errorconfints$Lower)
      highqvals <- -10.0*log10(errorconfints$Upper)
      
      qvals <- sapply(errorrates, function(x) { qv=-10.0*log10(x); return(qv) })
      points(mnlengths, qvals, xlim=c(10,40), cex=ptexp[i], pch=pchvalues[i], col=assemblycolors[i])
      if (errorbars) {
        arrows(x0=mnlengths, y0=lowqvals, x1=mnlengths, y1=highqvals, code=3, angle=90, length=0.1, col=assemblycolors[i])
      }
    }    
  }
  
  legend(28, 38, assemblylabels, col=assemblycolors, pch=pchvalues, pt.cex=pointcex*multiplevector)
  #legend(30, 30, assemblylabels, col=assemblycolors, pch=pchvalues)
}

# Make mononucleotide accuracy/QV plots

pdf("CompareAlignerOutput/AssemblyMononucAccuracy.pdf")  
assembly_mononucacc_plot(mnstatsfiles, assemblylabels)
dev.off()

pdf("CompareAlignerOutput/AssemblyMononucErrorQVs.pdf")  
assembly_mononucqv_plot(mnstatsfiles, assemblylabels, plottitle="Mononucleotide quality score by length for minimap2, winnowmap2")
dev.off()

# Indel errors:
indelstatsfiles <- sapply(assemblynames, function(x) {file=paste(c("CompareAlignerOutput/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})

totalindelrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")

  return(sum(indellengthhist$indelspermbaligned, na.rm=TRUE))  
}
totalinsertionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist[indellengthhist$indellength>0, "indelspermbaligned"], na.rm=TRUE))  
}
totaldeletionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist[indellengthhist$indellength<0, "indelspermbaligned"], na.rm=TRUE))  
}
assembly_indels_plot <- function(indellengthfiles, assemblylabels, xlabval="Assembly", ylabval="Indel Errors per mb", legendposxindex=NA, legendposy=NA, titleval="", ymax=NA) {
  totalinsertionrates <- sapply(indellengthfiles, totalinsertionrate)
  totaldeletionrates <- sapply(indellengthfiles, totaldeletionrate)

  assemblylabelswithinsdel <- sapply(assemblylabels, function(x) {label = paste(c(x, "\nIns/Del"), sep="", collapse=""); return(label)})
  
  barcolors <- sapply(assemblycolors, function(x) {c(darken(x), lighten(x, 0.3))})
  out <- barplot(rbind(totalinsertionrates, totaldeletionrates), beside=TRUE, names.arg=assemblylabelswithinsdel, col=barcolors, main=titleval, xlab=xlabval, ylab=ylabval)
  
  if (is.na(legendposxindex)) {
    legendposx = out[5]
  }
  else {
    legendposx = out[legendposxindex]
  }
  if (is.na(legendposy)) {
    legendposy = 13.9
  }
  legend(legendposx, legendposy, assemblylabels, col=assemblycolors, pch=15)
  formattedrates <- as.integer(rbind(totalinsertionrates, totaldeletionrates)*10)/10
  text(out, rbind(totalinsertionrates, totaldeletionrates), formattedrates, pos=3, xpd=NA)
  
}

# Make indels plot:

pdf("CompareAlignerOutput/IndelRates.pdf", width=16, height=8)  
assembly_indels_plot(indelstatsfiles, assemblylabels, titleval="Assembly indel errors for minimap2/winnowmap2", legendposxindex=11, legendposy=25)
dev.off()

### Phasing errors plot

phaseswitchfile <- "CompareAlignerOutput/assembly_switch_rates.txt"

assembly_switchrate_plot <- function(switchfile, plottitle="Phase switch errors per megabase") {
  switchrates <- read.table(switchfile, sep="\t", header=FALSE)
  names(switchrates) <- c("Assembly", "Hap", "SwitchRate")
  #switchrates <- switchrates[switchrates$Assembly != "Ash1v2.0", ]
  assemblylabelswithhap <- sapply(seq(1, length(switchrates$Assembly)), function(i) {label = paste(c(switchrates[i, "Assembly"], "\n(", switchrates[i, "Hap"], ")"), sep="", collapse=""); return(label)})
  #assemblycolorswithdups <- assemblycolors[find(assemblynames==switchrates$Assembly)]
  assemblycolorswithdups <- sapply(switchrates$Assembly, function(x) {assemblycolors[which(assemblylabels==x)]})
  
  # bottom, left, top, right  
  defaultmargin <- c(5.1, 4.1, 4.1, 2.1)
  par(mar = c(5.1, 7.1, 4.1, 4.1))
  #par(mar = c(5.1, 10.1, 4.1, 4.1))
  options(scipen=8)
  
  out <- barplot(log10(rev(switchrates$SwitchRate)), names.arg=rev(assemblylabelswithhap), las=1, horiz=TRUE, xaxt='n', col=rev(assemblycolorswithdups), xlab=c("Errors per megabase"), main=plottitle)
  xval <- c(1, 10, 100, 1000)
  xpos <- log10(xval)
  axis(1, xpos, xval, las=1)
  formattedrates <- as.integer(switchrates$SwitchRate*10)/10
  text(log10(rev(switchrates$SwitchRate)), out, rev(formattedrates), pos=4, xpd=NA)
  par(mar=defaultmargin)
  
}

pdf("CompareAlignerOutput/PhaseSwitchRates.pdf")  
assembly_switchrate_plot(phaseswitchfile, plottitle="Phase switch errors per megabase for minimap2, winnowmap2")
dev.off()

### Substitution rates

substitutionstatsfiles <- sapply(assemblynames, function(x) {file=paste(c("CompareAlignerOutput/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})

# Plot substitution rate-by-type histogram

typeorder <- c("A_C", "A_G", "A_T", "T_C", "T_G", "T_A", "G_A", "G_T", "G_C", "C_A", "C_T", "C_G")
titv <- c("tv", "ti", "tv", "ti", "tv", "tv", "ti", "tv", "tv", "tv", "ti", "tv" )

assembly_substitutions_plot <- function(assemblysubsfiles, assemblylabels, xlabval="Assembly", ylabval="Substitutions per mb", legendposxindex=NA, legendposy=NA, titleval="", ymax=NA, legend=TRUE) {

  firsthist <- read.table(assemblysubsfiles[1], sep="\t")
  names(firsthist) <- c("errortype", "errorcount", "errorspermbaligned")
  typeorderindex <- sapply(firsthist$errortype, function(x) {which(typeorder==x)})
  firsttis <- sum(firsthist[titv[typeorderindex]=="ti", "errorspermbaligned"])
  firsttvs <- sum(firsthist[titv[typeorderindex]=="tv", "errorspermbaligned"])
  tivals <- c(firsttis)
  tvvals <- c(firsttvs)
  assemblylabelswithtitv <- sapply(assemblylabels, function(x) {label = paste(c(x, "\nTi/Tv"), sep="", collapse=""); return(label)})
  
  for (i in seq(2, length(assemblysubsfiles))) {
    subshist <- read.table(assemblysubsfiles[i], sep="\t")
    names(subshist) <- c("errortype", "errorcount", "errorspermbaligned")
    typeorderindex <- sapply(subshist$errortype, function(x) {which(typeorder==x)})
    
    assemblytis <- sum(subshist[titv[typeorderindex]=="ti", "errorspermbaligned"])
    assemblytvs <- sum(subshist[titv[typeorderindex]=="tv", "errorspermbaligned"])
    
    tivals <- append(tivals, assemblytis)
    tvvals <- append(tvvals, assemblytvs)
  }
  
  barcolors <- sapply(assemblycolors, function(x) {c(darken(x), lighten(x, 0.3))})
  
  if (is.na(ymax)) {
    out <- barplot(rbind(tivals, tvvals), names.arg=assemblylabelswithtitv, beside=TRUE, col=barcolors, main=titleval, xlab=xlabval, ylab=ylabval)
  }
  else {
    out <- barplot(rbind(tivals, tvvals), names.arg=assemblylabelswithtitv, beside=TRUE, col=typecolors, main=titleval, xlab=xlabval, ylab=ylabval, ylim=c(0,ymax))
  }
  if (legend) {
    if (is.na(legendposxindex)) {
      legendposx = out[5]
    }
    else {
      legendposx = out[legendposxindex]
    }
    if (is.na(legendposy)) {
      legendposy = 13.9
    }
    
    legend(legendposx, legendposy, assemblylabels, col=assemblycolors, pch=15)
  }
  formattedrates <- as.integer(rbind(tivals, tvvals)*10)/10
  text(out, rbind(tivals, tvvals), formattedrates, pos=3, xpd=NA)
}

# Make substitutions plot:

pdf("CompareAlignerOutput/SubstitutionRates.pdf", width=16, height=8)  
assembly_substitutions_plot(substitutionstatsfiles, assemblylabels, titleval="Assembly substitution errors for minimap2/winnowmap2", legendposxindex=11, legendposy=180)
dev.off()

### Quality value (QV) scores for different assemblies by different methods
qvfile <- "CompareAlignerOutput/AssemblyQVScores.txt"

assembly_qv_plot <-function(assembly_stats_file, assemblylabels, titleval="Assembly QV by different methods") {

  assembly_stats <- read.table(assembly_stats_file, sep="\t", header=TRUE)
  names(assembly_stats) <- c("Assembly", "MerqQV", "YakQV", "BenchQV")
  
  minqv <- min(rbind(assembly_stats$MerqQV, assembly_stats$YakQV, assembly_stats$BenchQV), na.rm=TRUE)
  maxqv <- max(rbind(assembly_stats$MerqQV, assembly_stats$YakQV, assembly_stats$BenchQV), na.rm=TRUE)
  
  if (length(assemblylabels)==0) {
    assemblylabels=assembly_stats$Assembly
  }
  barplot(t(cbind(assembly_stats$MerqQV, assembly_stats$YakQV, assembly_stats$BenchQV)), beside=TRUE, names.arg=assembly_stats$Assembly, col=methodcolors, xlab="Assembly", ylim=c(0, maxqv+20), main=titleval)

  legend("topright", c("Merqury", "Yak", "Q100 Benchmark"), col=methodcolors, pch=15)
}

pdf("CompareAlignerOutput/QVMeasures.pdf", width=13, height=8.5)  
assembly_qv_plot(qvfile, assemblylabels)
dev.off()

### Plot them all together:
pdf("CompareAlignerMultiplot.pdf", width=13, height=9)
par(mfrow=c(2,3))
assembly_ngax_plot(assemblysizefiles, assemblylabels=assemblylabels, plottitle="NGAx for different assemblies")
#assembly_qv_plot(qvfile, assemblylabels)
assembly_mononucqv_plot(mnstatsfiles, assemblylabels, pointcex=1.5)
assembly_substitutions_plot(substitutionstatsfiles, assemblylabels, titleval="Substitution discrepancies in assemblies", legendposxindex=11, legendposy=180)
assembly_indels_plot(indelstatsfiles, assemblylabels, titleval="Indel discrepancies in assemblies", legendposxindex=11, legendposy=25)
assembly_switchrate_plot(phaseswitchfile, plottitle="Phase switch errors per mb")
dev.off()


