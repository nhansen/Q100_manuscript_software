setwd("/Users/nhansen/HG002_diploid_benchmark/plots/readerrors")

library(stringr)
args = commandArgs(trailingOnly=TRUE)

readsetname <- ifelse(!( is.na(args[1])), args[1], "elementins1000mat")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.0.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")

# Plot indel histogram

plotindellengths <- function(readsetname, outputdir, maxsize=10, xlabval="Length difference", ylabval="Indel Errors per mb", titleval="", ymax=NA, legend=TRUE) {
  plottitle <- ifelse(!( is.na(args[4])), args[4], titleval)
  indelfile <- paste(c(outputdir, "/", readsetname, ".indelerrorstats.txt"), sep="", collapse="")
  
  indellengthhist <- read.table(indelfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  insertionrates <- sapply(seq(1, maxsize), function(i) {indellengthhist[indellengthhist$indellength==i, "indelspermbaligned"]})
  deletionrates <- sapply(seq(-1, -1*maxsize, -1), function(i) {indellengthhist[indellengthhist$indellength==i, "indelspermbaligned"]})

  if (is.na(ymax)) {
    barplot(rbind(deletionrates, insertionrates), beside=TRUE, names.arg=seq(1,maxsize), col=c("blue", "red"), main=plottitle, xlab=xlabval, ylab=ylabval)
  }
  else {
    barplot(rbind(deletionrates, insertionrates), beside=TRUE, names.arg=seq(1,maxsize), col=c("blue", "red"), main=plottitle, xlab=xlabval, ylab=ylabval, ylim=c(0,ymax))
  }
  if (legend) {
    legend("topright", c("Deletions", "Insertions"), col=c("blue", "red"), pch=15)
  }
}

plotindellengths("ont_r10_duplex_split", outputdir)
plotindellengths("ont_q28_corrected", outputdir)
plotindellengths("ont_epi2me_q28", outputdir)
plotindellengths("hifidcv1.1", outputdir)
plotindellengths("elementins1000mat", outputdir)
plotindellengths("elementins1000pat", outputdir)
plotindellengths("element_avitilongmat", outputdir)

par(mfrow=c(2,2))
plotindellengths("element_avitilongmat", outputdir, titleval="Element Aviti (27 errors/Mb)")
plotindellengths("hifidcv1.1", outputdir, titleval="HiFi/DCv1.1 (1,557 errors/Mb)")
plotindellengths("hifirevio3cell", outputdir, titleval="HiFi Revio (1,335 errors/Mb)")
plotindellengths("ont_epi2me_q28", outputdir, titleval="ONT Q28 UL (1,515 errors/Mb)")
#plotindellengths("ont_r10_duplex_split", outputdir, titleval="ONT R10 Duplex (1,119/Mb)")

#currentpars <- par()
#par(mar=c(5.1, 4.1, 4.1, 2.1))
setwd("/Users/nhansen/HG002_diploid_polishing_validation/polishing_update_slides/comparison_may_2024/plots/benchresults")
plotindellengths("m64011_190830_dcv1.1_nonexc", outputdir, maxsize=5, titleval="HiFi SequelII/DCv1.1 (1,491 Indels/Mb)", ymax=800)
plotindellengths("hifi_revio_3cell_nonexc", outputdir, maxsize=5, titleval="HiFi Revio (1,333 Indels/Mb)", ymax=800, legend=FALSE)
plotindellengths("ONT_Q28", outputdir, maxsize=5, titleval="Epi2me Q28 (1,526 Indels/Mb)", ymax=800, legend=FALSE)
plotindellengths("ONT_R10_duplex_nonexc", outputdir, maxsize=5, titleval="R10 Duplex (1,450 Indels/Mb)", ymax=800, legend=FALSE)

plotindellengths("element_avitilongmat", outputdir, titleval="Element Aviti (23 errors/Mb)")
plotindellengths("m64011_190830_dcv1.1_nonexc", outputdir, titleval="HiFi/DCv1.1 (1,491 errors/Mb)")
plotindellengths("hifi_revio_3cell_nonexc", outputdir, titleval="HiFi Revio (1,333 errors/Mb)")
plotindellengths("ONT_Q28", outputdir, titleval="ONT Q28 UL (1,527 errors/Mb)")

plotindellengths("m64011_190830_dcv1.1_nonexc", outputdir, maxsize=5, titleval="HiFi SequelII/DCv1.1 (1,491 Indels/Mb)", ymax=800)
