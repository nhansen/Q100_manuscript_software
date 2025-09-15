
setwd("/Users/nhansen/HG002_diploid_benchmark/plots/ngaplots")
args = commandArgs(trailingOnly=TRUE)

clustersizefile <- ifelse(!( is.na(args[1])), args[1], "v2_trio.alignclusterlengths.txt")
contigsizefile <- ifelse(!( is.na(args[2])), args[2], "v2_trio.contiglengths.txt")
scaffoldsizefile <- ifelse(!( is.na(args[3])), args[3], "v2_trio.scaffoldlengths.txt")
genomename <- ifelse(!( is.na(args[2])), args[2], "v2_trio")
benchname <- ifelse(!( is.na(args[3])), args[3], "v1.0.1")
outputdir <- ifelse(!( is.na(args[4])), args[4], ".")

readlengths <- function(sizefile) {
  clusterlengths <- read.table(sizefile, sep="\t", header=FALSE)
  names(clusterlengths) <- c("perc", "clusterlength", "totallength", "chromosome")
  clusterlengths$clusterlength <- clusterlengths$clusterlength/1000000

  return(clusterlengths)
}

plotclusterlengths <- function(clusterlengths, color="red", title="NGAx", dashed=FALSE) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  ltyval <- ifelse(dashed, 2, 1)
  plot(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,100), ylim=c(0,250), xlab="Percent of Haploid Genome", ylab="Aligned length", main=title)
}

addclusterlengths <- function(clusterlengths, color="blue", dashed=FALSE) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  ltyval <- ifelse(dashed, 2, 1)
  lines(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,100))
}

plotsingleassemblyresults <- function(clusterfile, contigfile, scaffoldfile, ideal=TRUE, haplotype="MAT", plottitle="") {
  cluster_lengths <- readlengths(clusterfile)
  contig_lengths <- readlengths(contigfile)
  scaffold_lengths <- readlengths(scaffoldfile)

  if (haplotype=="MAT") {
    ideallengths <- readlengths("v1.0.1.mat.alignclusterlengths.txt")
  }
  else {
    ideallengths <- readlengths("v1.0.1.pat.alignclusterlengths.txt")
  }

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

plottriomatforserge <- function(ideal=FALSE, plottitle="NGAx with trio data, maternal haplotype") {
  hifiasm_hap2_trio <- readlengths("hifiasm_v19_trio.hap2.alignclusterlengths.txt")
  #hifiasm_hap2_ngx <- readlengths("hifiasm_v19_trio.hap2.scaffoldlengths.txt")
  hifiasm_hap2_ngx <- readlengths("hifiasm_v19_trio.hap2.contiglengths.txt")
  matideal <- readlengths("v1.0.1.mat.alignclusterlengths.txt")
  v2_hap1_trio <- readlengths("v2_trio.hap1.alignclusterlengths.txt")
  #v2_hap1_ngx <- readlengths("v2_trio.hap1.scaffoldlengths.txt")
  v2_hap1_ngx <- readlengths("v2_trio.hap1.contiglengths.txt")
  
  plotclusterlengths(hifiasm_hap2_trio, col="red", title=plottitle)
  addclusterlengths(v2_hap1_trio, col="blue")
  if (ideal) {
    addclusterlengths(matideal, col="black")
  }
  addclusterlengths(hifiasm_hap2_ngx, col="red", dashed=TRUE)
  addclusterlengths(v2_hap1_ngx, col="blue", dashed=TRUE)

  if (ideal) {
    legend("topright", c("ideal", "hifiasm NGAx", "hifiasm NGx", "verkkov2.0 NGAx", "verkkov2.0 NGx"), col=c("black", "red", "red", "blue", "blue"), lty=c(1, 1, 2, 1, 2))
  }
  else {
    legend("topright", c("hifiasm NGAx", "hifiasm NGx", "verkkov2.0 NGAx", "verkkov2.0 NGx"), col=c("red", "red", "blue", "blue"), lty=c(1, 2, 1, 2))
  }
}

plottriopatforserge <- function(ideal=FALSE, ngax="contig") {
  if (ngax=="contig") {
    hifiasm_hap1_ngx <- readlengths("hifiasm_v19_trio.hap1.contiglengths.txt")
  }
  else {
    hifiasm_hap1_ngx <- readlengths("hifiasm_v19_trio.hap1.scaffoldlengths.txt")
  }
  hifiasm_hap1_trio <- readlengths("hifiasm_v19_trio.hap1.alignclusterlengths.txt")
  patideal <- readlengths("v1.0.1.pat.alignclusterlengths.txt")
  v2_hap2_trio <- readlengths("v2_trio.hap2.alignclusterlengths.txt")
  if (ngax=="contig") {
    v2_hap2_ngx <- readlengths("v2_trio.hap2.contiglengths.txt")
  }
  else {
    v2_hap2_ngx <- readlengths("v2_trio.hap2.scaffoldlengths.txt")
  }
  
  plotclusterlengths(hifiasm_hap1_trio, col="red", title="NGAx with trio data, paternal haplotype")
  addclusterlengths(v2_hap2_trio, col="blue")
  if (ideal) {
    addclusterlengths(patideal, col="black")
  }
  addclusterlengths(hifiasm_hap1_ngx, col="red", dashed=TRUE)
  addclusterlengths(v2_hap2_ngx, col="blue", dashed=TRUE)

  if (ideal) {
    legend("topright", c("ideal", "hifiasm NGAx", "hifiasm NGx", "verkkov2.0 NGAx", "verkkov2.0 NGx"), col=c("black", "red", "red", "blue", "blue"), lty=c(1, 1, 2, 1, 2))
  }
  else {
    legend("topright", c("hifiasm NGAx", "hifiasm NGx", "verkkov2.0 NGAx", "verkkov2.0 NGx"), col=c("red", "red", "blue", "blue"), lty=c(1, 2, 1, 2))
  }
}

plothicmatforserge <- function(ideal=FALSE, ngax="contig") {
  #hifiasm_hap2_ngx <- readlengths("hifiasm_v19_hic.hap2.scaffoldlengths.txt")
  hifiasm_hap2_ngx <- readlengths("hifiasm_v19_hic.hap2.contiglengths.txt")
  hifiasm_hap2_hic <- readlengths("hifiasm_v19_hic.hap2.alignclusterlengths.txt")
  matideal <- readlengths("v1.0.1.mat.alignclusterlengths.txt")
  #v2_hap1_ngx <- readlengths("v2_hic.hap1.scaffoldlengths.txt")
  v2_hap1_ngx <- readlengths("v2_hic.hap1.contiglengths.txt")
  v2_hap1_hic <- readlengths("v2_hic.hap1.alignclusterlengths.txt")
  
  plotclusterlengths(hifiasm_hap2_hic, col="red", title="NGAx with hi-C data, maternal haplotype")
  addclusterlengths(v2_hap1_hic, col="blue")
  #addclusterlengths(matideal, col="black")
  addclusterlengths(hifiasm_hap2_ngx, col="red", dashed=TRUE)
  addclusterlengths(v2_hap1_ngx, col="blue", dashed=TRUE)

  #legend("topright", c("ideal", "hifiasm NGAx", "hifiasm NGx", "verkkov2.0 NGAx", "verkkov2.0 NGx"), col=c("black", "red", "red", "blue", "blue"), lty=c(1, 1, 2, 1, 2))
  legend("topright", c("hifiasm NGAx", "hifiasm NGx", "verkkov2.0 NGAx", "verkkov2.0 NGx"), col=c("red", "red", "blue", "blue"), lty=c(1, 2, 1, 2))
}

plothicpatforserge <- function(ideal=FALSE, ngax="contig") {
  if (ngax=="contig") {
    hifiasm_hap1_ngx <- readlengths("hifiasm_v19_hic.hap1.contiglengths.txt")
  }
  else {
    hifiasm_hap1_ngx <- readlengths("hifiasm_v19_hic.hap1.scaffoldlengths.txt")
  }
  hifiasm_hap1_hic <- readlengths("hifiasm_v19_hic.hap1.alignclusterlengths.txt")
  patideal <- readlengths("v1.0.1.pat.alignclusterlengths.txt")
  if (ngax=="contig") {
    v2_hap2_ngx <- readlengths("v2_hic.hap2.contiglengths.txt")
  }
  else {
    v2_hap2_ngx <- readlengths("v2_hic.hap2.scaffoldlengths.txt")
  }
  v2_hap2_hic <- readlengths("v2_hic.hap2.alignclusterlengths.txt")
  
  plotclusterlengths(hifiasm_hap1_hic, col="red", title="NGAx with hi-C data, paternal haplotype")
  addclusterlengths(v2_hap2_hic, col="blue")
  if (ideal) {
    addclusterlengths(patideal, col="black")
  }
  addclusterlengths(hifiasm_hap1_ngx, col="red", dashed=TRUE)
  addclusterlengths(v2_hap2_ngx, col="blue", dashed=TRUE)

  if (ideal) {
    legend("topright", c("ideal", "hifiasm", "verkkov2.0"), col=c("black", "red", "blue"), pch=15)
  }
  else {
    legend("topright", c("hifiasm NGAx", "hifiasm NGx", "verkkov2.0 NGAx", "verkkov2.0 NGx"), col=c("red", "red", "blue", "blue"), lty=c(1, 2, 1, 2))
  }
}

plotbogposterfigure <- function(ideal=FALSE, ngax="contig", haplotype="maternal") {
  if (ngax=="contig") {
    year1mat_ngx <- readlengths("year1mat.contiglengths.txt")
    year1pat_ngx <- readlengths("year1pat.contiglengths.txt")
  }
  else {
    year1mat_ngx <- readlengths("year1mat.scaffoldlengths.txt")
    year1pat_ngx <- readlengths("year1pat.scaffoldlengths.txt")
  }
  patideal <- readlengths("v1.0.1.pat.alignclusterlengths.txt")
  matideal <- readlengths("v1.0.1.mat.alignclusterlengths.txt")
  
  if (ngax=="contig") {
    v2_hap2_ngx <- readlengths("v2_trio.hap2.contiglengths.txt")
    v2_hap1_ngx <- readlengths("v2_trio.hap1.contiglengths.txt")
  }
  else {
    v2_hap2_ngx <- readlengths("v2_trio.hap2.scaffoldlengths.txt")
    v2_hap1_ngx <- readlengths("v2_trio.hap1.scaffoldlengths.txt")
  }

  if (haplotype=="maternal") {
    plotclusterlengths(year1mat_ngx, col="red", title="NGAx for HG002 Assemblies (Maternal haplotype)")
    addclusterlengths(v2_hap1_ngx, col="blue")
    if (ideal) {
      addclusterlengths(matideal, col="black")
    }
    
    if (ideal) {
      legend("topright", c("Year1 HPRC", "Verkkov2.0", "Ideal"), col=c("red", "blue", "black"), pch=15)
    }
    else {
      legend("topright", c("Year1 HPRC", "Verkkov2.0"), col=c("red", "blue"), pch=15)
    }    
  }
  else {
    plotclusterlengths(year1pat_ngx, col="red", title="NGAx for HG002 Assemblies (Paternal haplotype)")
    addclusterlengths(v2_hap2_ngx, col="blue")
    if (ideal) {
      addclusterlengths(patideal, col="black")
    }
    
    if (ideal) {
      legend("topright", c("Year1 HPRC", "Verkkov2.0", "Ideal"), col=c("red", "blue", "black"), pch=15)
    }
    else {
      legend("topright", c("Year1 HPRC", "Verkkov2.0"), col=c("red", "blue"), pch=15)
    }   
  }
}

plotbogposterfigure(ideal=TRUE, haplotype="maternal")

pdf("ngax_hic_pat.pdf")
plothicpatforserge()
dev.off()
pdf("ngax_hic_mat.pdf")
plothicmatforserge()
dev.off()
pdf("ngax_trio_pat.pdf")
plottriopatforserge()
dev.off()
pdf("ngax_trio_mat.pdf")
plottriomatforserge()
dev.off()


