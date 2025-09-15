setwd("/Users/nhansen/HG002_diploid_benchmark/PaperFigures")

################
### FIGURE 1 ###
### Adam's description: ###
### Overview of the complete diploid HG002 genome. ###
### [Diploid ideogram of all HG002 chromosomes, with maternal/paternal homologs ###
### adjacent to one another and each chromosome labeled with things like ###
### satellite annotation, het SNV density, het SVs, suspicious regions, gaps, etc.] ###
################

library(stringr)
library(karyoploteR)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

autosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX/Y")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

strandcolors <- c("#332288", "#88CCEE")

censat_color_palette <- c("#662F90",
                          "#FA99FF",
                          "#AC33C7",
                          "#00CCCC",
                          "#990000",
                          "#FF6600",
                          "#FF9200",
                          "#FFCC99",
                          "#00DE60",
                          "#1B998B",
                          "#0080FA",
                          "#335189",
                          #"#E0E0E0",
                          "#000000",
                          "#FFCC00",
                          "#CC0000",
                          "#78A1BB",
                          "#FF6600",
                          "#FFCC99",
                          "#990000")
censat_labels <- c("rDNA",
  "bsat",
  "gsat",
  "censat",
  "activeHOR",
  "inactiveHOR",
  "dHOR",
  "mon",
  "HSat1A",
  "HSat1B",
  "HSat2",
  "HSat3",
  #"ct",
  "GAP",
  "subTerm",
  "mixedAlpha",
  "HSat2_3", 
  "hor(SF3-5)",
  "mon/hor(SF4)",
  "active_hor")

genomefile <- "KaryoploteRFiles/v1.1.karyotype.txt"
matgenomefile <- "KaryoploteRFiles/v1.1.karyotype.mat.txt"
patgenomefile <- "KaryoploteRFiles/v1.1.karyotype.pat.txt"
benchgenome <- toGRanges(genomefile)
matbenchgenome <- toGRanges(matgenomefile)
patbenchgenome <- toGRanges(patgenomefile)

matgiabconfidentfile <-"KaryoploteRFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.liftedtov1.1.mat.bed"
patgiabconfidentfile <-"KaryoploteRFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.liftedtov1.1.pat.bed"
matgiabconf <- toGRanges(matgiabconfidentfile)
patgiabconf <- toGRanges(patgiabconfidentfile)

matcensatfile <- "KaryoploteRFiles/hg002v1.1.cenSatv2.0.MAT.bed"
matcensats <- toGRanges(matcensatfile)
matcensatdf <- read.table(matcensatfile, sep="\t", header=FALSE)
matrdnasdf <- matcensatdf[matcensatdf$V4=="rDNA",]
matrdnas <- toGRanges(matrdnasdf)

patcensatfile <- "KaryoploteRFiles/hg002v1.1.cenSatv2.0.PAT.bed"
patcensats <- toGRanges(patcensatfile)
patcensatdf <- read.table(patcensatfile, sep="\t", header=FALSE)
patrdnasdf <- patcensatdf[patcensatdf$V4=="rDNA",]
patrdnas <- toGRanges(patrdnasdf)

allcensats <- c(matcensats, patcensats)

hetfile <- "KaryoploteRFiles/v1.1.hets.non1to1.200k.hapaligned.windows.withheatmap.bed"
heterozygosity <- toGRanges(hetfile)

mathetfile <- "KaryoploteRFiles/v1.1.hets.non1to1.200k.hapaligned.windows.withheatmap.mat.bed"
pathetfile <- "KaryoploteRFiles/v1.1.hets.non1to1.200k.hapaligned.windows.withheatmap.pat.bed"
matheterozygosity <- toGRanges(mathetfile)
patheterozygosity <- toGRanges(pathetfile)

excludedfile <- "KaryoploteRFiles/v1.1.excluded_regions.bed"
excluded <- toGRanges(excludedfile)

matexcludedfile <- "KaryoploteRFiles/v1.1.excluded_regions.mat.bed"
matexcluded <- toGRanges(matexcludedfile)

patexcludedfile <- "KaryoploteRFiles/v1.1.excluded_regions.pat.bed"
patexcluded <- toGRanges(patexcludedfile)

all_aligns <- "KaryoploteRFiles/v1.1.hets.non1to1.200k.hapaligned.bed"
all_align_df <- read.table(all_aligns, sep="\t", header=FALSE)
names(all_align_df) <- c("chrom1", "start1", "end1", "name", "score", "strand")
longall_align_df <- all_align_df[all_align_df$end1 - all_align_df$start1 >=100000, ]
longall_align_df$chrom2 <- gsub("([^:]+):.*", "\\1", longall_align_df$name, perl = TRUE)
longall_align_df$start2 <- as.integer(gsub("[^:]+:(\\d+).*", "\\1", longall_align_df$name, perl = TRUE))
longall_align_df$end2 <- as.integer(gsub("[^:]+:\\d+\\-(\\d+).*", "\\1", longall_align_df$name, perl = TRUE))
longall_align_df$strand <- ifelse(longall_align_df$strand=="F", '+', '-')

matlongaligndf <- longall_align_df[str_detect(longall_align_df$chrom1, '^.*MATERNAL.*$' ), ]
patlongaligndf <- longall_align_df[str_detect(longall_align_df$chrom1, '^.*PATERNAL.*$' ), ]

plot_chrom_pairs_het_chroms <- function(genome=benchgenome, matgenome=matbenchgenome, patgenome=patbenchgenome, linkpanel1=1, linkpanel2=2, aligndf=matlongaligndf, plotchroms=NA, chromlinkname=NA, mathet=matheterozygosity, pathet=patheterozygosity, matcentro=matcensats, patcentro=patcensats, centerlow=0.2, centerhigh=0.4, ideoheight=200, matloqual=matexcluded, patloqual=patexcluded, loquallow=0.1, loqualhigh=0.4, invert=FALSE, labelsat=NA, labelwidths=NA, chromlabels=NA, plottitle="hg002v1.1") {
  allchroms <- rev(seqlevels(genome))
  if (length(plotchroms)>1 || !is.na(plotchroms)) {
    plotchroms <- plotchroms
  }
  else {
    plotchroms <- allchroms
  }

  if (!is.na(chromlinkname)) {
    regexstring <- paste(c("^", chromlinkname, "_"), sep="", collapse="")
    if (chromlinkname=="chrX") {
      linkchroms <- c("chrX_MATERNAL", "chrY_PATERNAL")
    }
    else {
      linkchroms <- allchroms[grep(regexstring, allchroms)]
    }
  }
  pp <- getDefaultPlotParams(plot.type=2)
  pp$ideogramheight <- ideoheight
  pp$data1height <- 2000
  pp$data2height <- 2000
  #pp$topmargin <- 650
  
  kp <- plotKaryotype(genome=genome, chromosomes=plotchroms, plot.type=2, lwd=0.0, ideogram.plotter=NULL, labels.plotter=NULL, cex=0.4, plot.params=pp)

  # plot the het tracks between the ideograms (panel 1 for mat, panel 2 for pat)
  kpRect(kp, data=matgenome, col="gray", border="white", data.panel=1, lwd=0.0, y0=rep(centerlow, length(matgenome)), y1=rep(centerhigh, length(matgenome)))
  kpRect(kp, data=patgenome, col="gray", border="white", data.panel=2, lwd=0.0, y0=rep(centerlow, length(patgenome)), y1=rep(centerhigh, length(patgenome)))
  kpRect(kp, data=mathet, col=mathet$itemRgb, data.panel=1, lwd=0.0, border="white", y0=rep(centerlow, length(mathet)), y1=rep(centerhigh, length(mathet)))
  kpRect(kp, data=pathet, col=pathet$itemRgb, data.panel=2, border="white", lwd=0.0, y0=rep(centerlow, length(pathet)), y1=rep(centerhigh, length(pathet)))
  kpRect(kp, data=matrdnas, col="black", data.panel=1, lwd=0.0, border="white", y0=rep(centerlow, length(matrdnas)), y1=rep(centerhigh, length(matrdnas)))
  kpRect(kp, data=patrdnas, col="black", data.panel=2, border="white", lwd=0.0, y0=rep(centerlow, length(patrdnas)), y1=rep(centerhigh, length(patrdnas)))
  
  # if possible, plot alignment bands between chromosomes
  if (length(aligndf) > 1 || !is.na(aligndf)) {
    alignstarts <-toGRanges(aligndf[,c("chrom1", "start1", "end1")])
    alignends <- toGRanges(aligndf[,c("chrom2", "start2", "end2")])
    strand(alignends) <- aligndf$strand
    forwardcolor <- paste(strandcolors[2], 'CC', sep="", collapse="")
    reversecolor <- paste(strandcolors[1], 'CC', sep="", collapse="")
    aligncolors <- ifelse(aligndf$strand=="+", forwardcolor, reversecolor)
    #kpRect(kp, data=matgenome, col="#E0E0E0", data.panel=1, lwd=0.0, y0=rep(centerhigh, length(matgenome)), y1=rep(1, length(matgenome)))
    #kpRect(kp, data=patgenome, col="#E0E0E0", data.panel=2, lwd=0.0, y0=rep(centerhigh, length(patgenome)), y1=rep(1, length(patgenome)))
    kpPlotLinks(kp, data=alignstarts, data2=alignends, col="white", lwd=0.05, border="black", r0=centerhigh, r1=centerhigh, data.panel=linkpanel1, data.panel2=linkpanel2, y=0)
  }
  
  ##GIAB regions plotted on ideograms:
  #kpRect(kp, data=matgiabconf, col="black", border="white", data.panel="ideogram", lwd=0.0, y0=rep(0, length(matgiabconf)), y1=rep(1, length(matgiabconf)))
  #kpRect(kp, data=patgiabconf, col="black", border="white", data.panel="ideogram", lwd=0.0, y0=rep(0, length(patgiabconf)), y1=rep(1, length(patgiabconf)))

  ##cenSat annotations plotted outside ideograms:
  #kpRect(kp, data=matcensats, col=matcensats$itemRgb, data.panel=2, lwd=0.0, border="white", y0=rep(centerlow, length(matcensats)), y1=rep(centerhigh, length(matcensats)))
  #kpRect(kp, data=patcensats, col=patcensats$itemRgb, data.panel=1, border="white", lwd=0.0, y0=rep(centerlow, length(patcensats)), y1=rep(centerhigh, length(patcensats)))

  #cenSat annotations plotted on ideograms:  
  kpRect(kp, data=matcensats, col=matcensats$itemRgb, data.panel="ideogram", lwd=0.0, y0=rep(0, length(matcensats)), y1=rep(1, length(matcensats)))
  kpRect(kp, data=patcensats, col=patcensats$itemRgb, data.panel="ideogram", lwd=0.0, y0=rep(0, length(patcensats)), y1=rep(1, length(patcensats)))
  mtext("hg002v1.1", side = 2, outer = TRUE)
  #mtext(paste(chroms_to_plot, sep=" ", collapse="   "), side=2, outer=FALSE, cex=0.4)

  if (!is.na(labelsat) && !is.na(labelwidths)) {
     labelpositions <- seq(0, length(chromlabels)-1)*labelwidths + labelsat
  }
  if (length(chromlabels)>1 || !is.na(chromlabels)) {
    mtext(chromlabels, side=2, at=labelpositions, line=2, outer=FALSE, cex=0.4)
  }
  mtext(plottitle, side = 2, outer = TRUE, line=-1)
}

plot_chrom_pairs_censat_chroms <- function(genome=benchgenome, matgenome=matbenchgenome, patgenome=patbenchgenome, linkpanel1=1, linkpanel2=2, aligndf=matlongaligndf, plotchroms=NA, chromlinkname=NA, mathet=matheterozygosity, pathet=patheterozygosity, matcentro=matcensats, patcentro=patcensats, centerlow=0.2, centerhigh=0.4, ideoheight=200, matloqual=matexcluded, patloqual=patexcluded, loquallow=0.1, loqualhigh=0.4, invert=FALSE, labelsat=NA, labelwidths=NA, chromlabels=NA, plottitle="hg002v1.1") {
  allchroms <- rev(seqlevels(genome))
  if (length(plotchroms)>1 || !is.na(plotchroms)) {
    plotchroms <- plotchroms
  }
  else {
    plotchroms <- allchroms
  }
  
  if (!is.na(chromlinkname)) {
    regexstring <- paste(c("^", chromlinkname, "_"), sep="", collapse="")
    if (chromlinkname=="chrX") {
      linkchroms <- c("chrX_MATERNAL", "chrY_PATERNAL")
    }
    else {
      linkchroms <- allchroms[grep(regexstring, allchroms)]
    }
  }
  pp <- getDefaultPlotParams(plot.type=2)
  pp$ideogramheight <- ideoheight
  #pp$data1height <- 2000
  #pp$data2height <- 2000
  #pp$topmargin <- 650
  
  kp <- plotKaryotype(genome=genome, chromosomes=plotchroms, plot.type=2, lwd=0.0, ideogram.plotter=NULL, labels.plotter=NULL, cex=0.4, plot.params=pp)
  
  # plot the het tracks between the ideograms (panel 1 for mat, panel 2 for pat)
  kpRect(kp, data=matgenome, col="gray", data.panel=1, lwd=0.0, y0=rep(centerlow, length(matgenome)), y1=rep(centerhigh, length(matgenome)))
  kpRect(kp, data=patgenome, col="gray", data.panel=2, lwd=0.0, y0=rep(centerlow, length(patgenome)), y1=rep(centerhigh, length(patgenome)))
  kpRect(kp, data=mathet, col=mathet$itemRgb, data.panel=1, lwd=0.0, y0=rep(centerlow, length(mathet)), y1=rep(centerhigh, length(mathet)))
  kpRect(kp, data=pathet, col=pathet$itemRgb, data.panel=2, lwd=0.0, y0=rep(centerlow, length(pathet)), y1=rep(centerhigh, length(pathet)))
  
  # if possible, plot alignment bands between chromosomes
  if (length(aligndf) > 1 || !is.na(aligndf)) {
    alignstarts <-toGRanges(aligndf[,c("chrom1", "start1", "end1")])
    alignends <- toGRanges(aligndf[,c("chrom2", "start2", "end2")])
    strand(alignends) <- aligndf$strand
    forwardcolor <- paste(strandcolors[2], 'CC', sep="", collapse="")
    reversecolor <- paste(strandcolors[1], 'CC', sep="", collapse="")
    aligncolors <- ifelse(aligndf$strand=="+", forwardcolor, reversecolor)
    kpPlotLinks(kp, data=alignstarts, data2=alignends, col=NA, border="black", r0=centerhigh, r1=centerhigh, data.panel=linkpanel1, data.panel2=linkpanel2, y=0)
  }  
  # plot the centromere track on top of the ideograms:
  #kpRect(kp, data=matgenome, col="grey", border="black", data.panel="ideogram", lwd=0.01, y0=rep(0, length(matgenome)), y1=rep(1, length(matgenome)))
  #kpRect(kp, data=patgenome, col="grey", border="black", data.panel="ideogram", lwd=0.01, y0=rep(0, length(patgenome)), y1=rep(1, length(patgenome)))
  kpRect(kp, data=matcensats, col=matcensats$itemRgb, data.panel="ideogram", lwd=0.0, y0=rep(0, length(matcensats)), y1=rep(1, length(matcensats)))
  kpRect(kp, data=patcensats, col=patcensats$itemRgb, data.panel="ideogram", lwd=0.0, y0=rep(0, length(patcensats)), y1=rep(1, length(patcensats)))
  mtext("hg002v1.1", side = 2, outer = TRUE)
  #mtext(paste(chroms_to_plot, sep=" ", collapse="   "), side=2, outer=FALSE, cex=0.4)
  
  if (!is.na(labelsat) && !is.na(labelwidths)) {
    labelpositions <- seq(0, length(chromlabels)-1)*labelwidths + labelsat
  }
  if (length(chromlabels)>1 || !is.na(chromlabels)) {
    mtext(chromlabels, side=2, at=labelpositions, line=1, outer=FALSE, cex=0.4)
  }
  mtext(plottitle, side = 2, outer = TRUE, line=-1)
}

pdf('Figure1.083024.pdf')
plot_chrom_pairs_het_chroms(benchgenome, mathet=matheterozygosity, pathet=patheterozygosity, matcentro=matcensats, patcentro=patcensats, centerlow=0.1, centerhigh=0.5, linkpanel1=1, linkpanel2=2, ideoheight=800, chromlabels=autosomes, labelsat=-0.175, labelwidth=0.058)

dev.off()

plot_karyoplot_legends <- function() {
  dataframe <- data.frame(X=rep(0, length(censat_labels)), Y=rep(1, length(censat_labels)), CentromereSequence=factor(censat_labels, levels = censat_labels), itemRgb=censat_color_palette)
  gplot <- ggplot(dataframe, aes(X, Y, fill=CentromereSequence)) + geom_point(size = 7, shape=22) + scale_fill_manual(values=censat_color_palette) + scale_color_manual(values=censat_color_palette) + guides(fill=guide_legend(ncol=4))
  legend <- get_legend(gplot)                     
  
  # Create new plot window 
  grid.newpage()                               
  
  # Draw Only legend  
  grid.draw(legend)  
}

pdf('Figure1Legend.083124.pdf')
plot_karyoplot_legends()
dev.off()

plot_het_heatmap_legend <- function() {
  col_fun = colorRamp2(c(0.0001, 0.0015, 0.003), c("yellow", "orange", "red")) 
  lgd = Legend(col_fun = col_fun)
  #draw(lgd)
  grid.newpage()
  #grid.draw()
  draw(lgd, x = unit(0.9, "npc"), y = unit(0.5, "npc"))
  #legend(x="right", legend=c("min", "med", "max"), fill=heat.colors(3))
}

#plot_chrom_pairs(benchgenome, mathet=matheterozygosity, pathet=patheterozygosity, matcentro=matcensats, patcentro=patcensats)

# Can we make a neat plot of the chromosome 8 inverted region?

plot_chromosome_alignments <- function(genome, chromname) {
  regexstring <- paste(c("^", chromname, "_"), sep="", collapse="")
  allchroms <- seqlevels(genome)
  if (chromname=="chrX") {
    chroms <- c("chrX_MATERNAL", "chrY_PATERNAL")
  }
  else {
    chroms <- allchroms[grep(regexstring, allchroms)]
  }
  
  pp <- getDefaultPlotParams(plot.type=1)
  pp$ideogramheight <- 0
  pp$data1height <- 0
  pp$data2height <- 0
  pp$data1inmargin <- 0
  
  chrom_align_df <- longall_align_df[grep(regexstring, longall_align_df$matchrom), ] 

  kp <- plotKaryotype(genome=genome, chromosomes=chroms, plot.type=1, cex=0.4, plot.params=pp)
  kpAddBaseNumbers(kp, tick.dist = 20000000, cex=0.3, minor.ticks=FALSE, tick.len=1)
  inv1start <-toGRanges(chrom_align_df[,c("matchrom", "matstart", "matend")])
  inv1end <- toGRanges(chrom_align_df[,c("patchrom", "patstart", "patend")])
  strand(inv1end) <- chrom_align_df$strand
  forwardcolor <- paste(strandcolors[2], 'CC', sep="", collapse="")
  reversecolor <- paste(strandcolors[1], 'CC', sep="", collapse="")
  aligncolors <- ifelse(chrom_align_df$strand=="+", forwardcolor, reversecolor)
  kpPlotLinks(kp, data=inv1start, data2=inv1end, col=aligncolors, r0=0, r1=0)
}

plot_chromosome_alignments(benchgenome, "chr8")

######
#' kpPlotLinks
#' 
#' @description 
#' 
#' Given 2 \code{GRanges} objects, plot lines or ribbons between region pairs
#' 
#' @details 
#'  
#'  This is one of the high-level, or specialized, plotting functions of karyoploteR.
#'  It takes two \code{GRanges} objects (or a single specially crafted one) and
#'  plots links (either lines or ribbons) between region pairs. Links are 
#'  plotted bewteen the first region of both objects, between the second one, etc...
#'  and therefore both objects need to have the same length. Specifying a region
#'  as negative strand, will "flip" it, so the the start of a region can be 
#'  linked to the end of its pair.
#'  
#' @note 
#'   For a link to be plotted BOTH ends must be visible in the karyoplot. In 
#'   particular, if a chromosome is not included in the plot (due to not
#'   being specified in \code{chromosomes}, for example) any link with an end
#'   on it will NOT be plotted. The same is true for zoomed in plots, where only
#'   intrachromosomal links will be visible. No warning or message will be
#'   generated.
#' 
#'
#' @usage kpPlotLinks(karyoplot, data, data2=NULL, y=0, arch.height=NULL, data.panel=1, r0=NULL, r1=NULL, ymin=NULL, ymax=NULL, col="#8e87eb", border=NULL, clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object with link start regions. If data2 is NULL, mcols(data) should be a bed-like structure with "link.chr", "link.start", "link.end" and optionally a "link.strand" columns. The first thee columns can have any name and the strand information will be extracted from the first column with "strand" in its name.
#' @param data2  (a \code{GRanges}) A GRanges object with the link end regions. If null, the end of the regions will be extracted from mcols(data). (Defaults to NULL)
#' @param y    (numeric) The y value where the origin and end of the links should be plotted (Defaults to 0)
#' @param arch.height    (numeric) The approximate arch height in links in the same chromosome in "y" scale. If NULL, it defaults to the whole span of the data panel.Also affects the curvature of links between chromosomes (Defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ymin    (numeric) The minimum value to be plotted on the data panel. If NULL, it is set to 0. (deafults to NULL)
#' @param ymax    (numeric) The maximum value to be plotted on the data.panel. If NULL the maximum density is used. (defaults to NULL)
#' @param col    (color) The background color of the links. If NULL and border is specified, it defaults to a lighter version of border.
#' @param border  (color) The border color of the links. If NULL and col is specified, it defaults to a darker version of col.
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
#' 
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpPlotRibbon}}, \code{\link{kpSegments}}
#' 
#' @examples
#'  
#'  
#'  set.seed(222)
#'  
#'  starts <- sort(createRandomRegions(nregions = 15))
#'  ends <- sort(createRandomRegions(nregions = 15))
#'  
#'  kp <- plotKaryotype()
#'  kpPlotLinks(kp, data=starts, data2=ends)
#'  
#'  #Create larger regions, so they look like ribbons
#'  starts <- sort(createRandomRegions(nregions = 15, length.mean = 8e6, length.sd = 5e6))
#'  ends <- sort(createRandomRegions(nregions = 15, length.mean = 8e6, length.sd = 5e6))
#'  
#'  kp <- plotKaryotype()
#'  kpPlotLinks(kp, data=starts, data2=ends)
#'  
#'  #flip some of them to represent inversions
#'  strand(ends) <- sample(c("+", "-"), length(ends), replace = TRUE)
#'  
#'  kp <- plotKaryotype()
#'  kpPlotLinks(kp, data=starts, data2=ends)
#'  
#'
'@importFrom bezier bezier
#'  
#'@export kpPlotLinks


kpPlotLinks <- function(karyoplot, data, data2=NULL, y=0, arch.height=NULL, data.panel=1, data.panel2=1, r0=NULL, r1=NULL, ymin=NULL, ymax=NULL, col="#8e87eb", border=NULL, clipping=TRUE, ...) {
  #Check parameters
  #karyoplot
  if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
  if(missing(data)) stop("The parameter 'data' is required")
  if(!methods::is(data, "GRanges") && !methods::is(data, "SimpleRleList")) stop("'data' must be a GRanges object or a SimpleRleList")
  
  #Assign complex default values
  #arc.height
  if(is.null(arch.height)) {
    arch.height <- karyoplot$plot.params[[paste0("data", data.panel, "max")]] - karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  }
  
  #colors
  prep.cols <- karyoploteR:::preprocessColors(col=col, border=border)
  col <- prep.cols$col
  border <- prep.cols$border
  
  
  #Define the link ends
  if(!is.null(data2) && !any(is.na(data2))) {
    if(!methods::is(data2, "GRanges")) {
      stop("If present, data2 must be a GRanges object")
    } else {
      if(length(data) != length(data2)) {
        stop("data and data2 must have the same length")
      }
    }
  } else {
    #if data2 is not defined, try to define it from the mcols of data
    tryCatch(data2 <- toGRanges(data.frame(mcols(data))),
             error=function(e) {stop("It was not possible to create data2 from mcols(data). ", e)})
    strand.col <- which(grepl(pattern = "strand", names(mcols(data2))))[1]
    if(length(strand.col)>0) {
      strand(data2) <- mcols(data2)[[strand.col]]
    }
  }
  
  
  #remove any links with at least one end out of the plotted regions
  to.keep <- overlapsAny(data, karyoplot$genome) & overlapsAny(data2, karyoplot$genome)
  if(any(!to.keep)) {
    data <- data[to.keep]
    data2 <- data2[to.keep]
  }
  
  
  if(length(data)==0) invisible(karyoplot) #return fast if no links will be plotted
  
  #Filter the additional arguments
  #NOTE: in this case we do not use the filter returned by prepareParameters
  #because we have a specific filtering logic due to links having two genomic 
  #regions
  dots <- filterParams(list(...), to.keep, length(to.keep))
  col <- filterParams(col, to.keep, length(to.keep))
  border <- filterParams(border, to.keep, length(to.keep))
  y <- filterParams(y, to.keep, length(to.keep))
  arch.height <- filterParams(arch.height, to.keep, length(to.keep))
  
  
  
  
  #Prepare the coordinates
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  ccf <- karyoplot$coord.change.function
  
  #Transform the coordinates of the starts
  pp.start <- prepareParameters4("kpPlotLinks", karyoplot=karyoplot, data=data, chr=NULL, x0=NULL, x1=NULL,
                                 y0=y, y1=y, ymin=ymin, ymax=ymax, r0=r0, r1=r1,
                                 data.panel=data.panel, ...)
  
  #Transform the coordinates of the starts
  pp.end <- prepareParameters4("kpPlotLinks", karyoplot=karyoplot, data=data2, chr=NULL, x0=NULL, x1=NULL,
                               y0=y, y1=y, ymin=ymin, ymax=ymax, r0=r0, r1=r1,
                               data.panel=data.panel2, ...)
  
  x0.start <- ccf(chr=pp.start$chr, x=pp.start$x0, data.panel=data.panel)$x
  x1.start <- ccf(chr=pp.start$chr, x=pp.start$x1, data.panel=data.panel)$x
  y.start <- ccf(chr=pp.start$chr, y=pp.start$y0, data.panel=data.panel)$y
  
  #swap the order of the regions in the negative strand
  aux <- numeric(length(data))
  neg.strand.start <- which(as.logical(strand(data)=="-"))
  aux[neg.strand.start] <- x0.start[neg.strand.start]
  x0.start[neg.strand.start] <- x1.start[neg.strand.start]
  x1.start[neg.strand.start] <- aux[neg.strand.start]
  
  x0.end <- ccf(chr=pp.end$chr, x=pp.end$x0, data.panel=data.panel2)$x
  x1.end <- ccf(chr=pp.end$chr, x=pp.end$x1, data.panel=data.panel2)$x
  y.end <- ccf(chr=pp.end$chr, y=pp.end$y0, data.panel=data.panel2)$y
  
  #swap the order of the regions in the negative strand
  neg.strand.end <- which(as.logical(strand(data2)=="-"))
  aux[neg.strand.end] <- x0.end[neg.strand.end]
  x0.end[neg.strand.end] <- x1.end[neg.strand.end]
  x1.end[neg.strand.end] <- aux[neg.strand.end]
  
  
  #transform the arch height to the plot coords using ccf and prepare parameters2 to take into account r0 and r1, ymin and ymax, etc
  y.min <- prepareParameters2("kpPlotLinks", karyoplot = karyoplot, data=NULL, chr=karyoplot$chromosomes[1], x=0, y=0, r0=r0, r1=r1, ymin=ymin, ymax=ymax, data.panel=data.panel)$y
  y.max <- prepareParameters2("kpPlotLinks", karyoplot = karyoplot, data=NULL, chr=karyoplot$chromosomes[1], x=0, y=arch.height, r0=r0, r1=r1, ymin=ymin, ymax=ymax, data.panel=data.panel)$y
  
  y.ctrl <- abs(ccf(chr=rep(karyoplot$chromosomes[1], length(y.min)), x=0, y=y.min, data.panel=data.panel)$y - ccf(chr=rep(karyoplot$chromosomes[1], length(y.max)), x=0, y=y.max, data.panel=data.panel)$y)
  
  #recycle parameters
  col <- karyoploteR:::recycle.first(col, data)
  border <- karyoploteR:::recycle.first(border, data)
  y.ctrl <- karyoploteR:::recycle.first(y.ctrl, data)
  
  #TODO: Make this 50 a parameter (num.bezier.segments=50) and document
  t <- seq(0, 1, length=50)
  for(i in seq_along(data)) {
    #The position of the control points (above or below the start and end points) depend on the relative position of start and end and in some cases on the data.panel
    if(y.start[i]>y.end[i]) { #if the start is above the end
      bezier_points_1 <- bezier:::bezier(t=t, p=list(c(x0.start[i],x0.start[i], x0.end[i], x0.end[i]), c(y.start[i],y.start[i]-y.ctrl[i],y.end[i]+y.ctrl[i],y.end[i])))
      bezier_points_2 <- bezier:::bezier(t=t, p=list(c(x1.start[i],x1.start[i], x1.end[i], x1.end[i]), c(y.start[i],y.start[i]-y.ctrl[i],y.end[i]+y.ctrl[i],y.end[i])))
    } else if(y.start[i]<y.end[i]) { #if start is below the end
      bezier_points_1 <- bezier:::bezier(t=t, p=list(c(x0.start[i],x0.start[i], x0.end[i], x0.end[i]), c(y.start[i],y.start[i]+y.ctrl[i],y.end[i]-y.ctrl[i],y.end[i])))
      bezier_points_2 <- bezier:::bezier(t=t, p=list(c(x1.start[i],x1.start[i], x1.end[i], x1.end[i]), c(y.start[i],y.start[i]+y.ctrl[i],y.end[i]-y.ctrl[i],y.end[i])))
    } else {
      #if they are in the same chromosome, detect whether we are in an upward or downward data panel
      if(ccf(chr=as.character(seqnames(data[1])), x=0, y=0, data.panel=data.panel)$y <
         ccf(chr=as.character(seqnames(data[1])), x=0, y=1, data.panel=data.panel)$y) {
        #it's upwards
        bezier_points_1 <- bezier:::bezier(t=t, p=list(c(x0.start[i],x0.start[i], x0.end[i], x0.end[i]), c(y.start[i],y.start[i]+y.ctrl[i],y.end[i]+y.ctrl[i],y.end[i])))
        bezier_points_2 <- bezier:::bezier(t=t, p=list(c(x1.start[i],x1.start[i], x1.end[i], x1.end[i]), c(y.start[i],y.start[i]+y.ctrl[i],y.end[i]+y.ctrl[i],y.end[i])))
      } else {
        #it's downwards
        bezier_points_1 <- bezier:::bezier(t=t, p=list(c(x0.start[i],x0.start[i], x0.end[i], x0.end[i]), c(y.start[i],y.start[i]-y.ctrl[i],y.end[i]-y.ctrl[i],y.end[i])))
        bezier_points_2 <- bezier:::bezier(t=t, p=list(c(x1.start[i],x1.start[i], x1.end[i], x1.end[i]), c(y.start[i],y.start[i]-y.ctrl[i],y.end[i]-y.ctrl[i],y.end[i])))
      }
    }
    
    processClipping(karyoplot=karyoplot, clipping=clipping, data.panel=data.panel)  
    graphics::polygon(x = c(bezier_points_1[,1], rev(bezier_points_2[,1])), y=c(bezier_points_1[,2], rev(bezier_points_2[,2])), col=col[i], border=NA, ...)
    if(!is.na(border[i])) {
      graphics::lines(bezier_points_1, col=border[i], ...)
      graphics::lines(bezier_points_2, col=border[i], ...)
    }
  }
  
  
  invisible(karyoplot)
}

