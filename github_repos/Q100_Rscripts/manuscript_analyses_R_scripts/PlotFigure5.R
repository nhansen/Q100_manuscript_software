setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/PaperFigures/Figures")

library(colorspace)
library(Hmisc)

source("/Users/nhansen/OneDrive/HG002_diploid_benchmark/Q100_Rscripts/manuscript_analyses_R_scripts/ReadBenchComparisonPlotFunctions.R")

################
### FIGURE 5 ###
### Sequencing coverage and accuracy of various platforms compared against ###'
### the HG002 genome benchmark. [Essential read quality statistics from ###
### Nanopore, PacBio, Illumina, Element, etc. illustrating subtle base ###
### call and coverage biases not evident from typically reported QV scores. ###
### Things like real vs. reported QV values, indel biases, coverage ###
### uniformity, transition/transversion biases ###
################

outdir <- "Figure4Output"

# currently plotted read sets (alternatives below)
readsetnames <- c("ont_epi2me_q28", "hifi_revio_pbmay24", "element_ultraq_jun2024", "illumina_googlepcrfree")
platformlabels <- c("ONT Q28", "HiFi Revio", "Element", "Illumina")

plusreadsetnames <- c("ont_epi2me_q28", "hifi_revio_pbmay24", "element_ultraq_jun2024", "illumina_googlepcrplus")
plusplatformlabels <- c("ONT Q28", "HiFi Revio", "Element", "Illumina (PCR Plus)")

arangsreadsetnames <- c("ont_epi2me_q28", "hifi_revio_pbmay24", "element_ultraq_jun2024", "illumina_googlepcrfree", "Ultima_2025apr")
arangsplatformlabels <- c("ONT Q28", "HiFi Revio", "Element", "Illumina", "Ultima")

#doublereadsetnames <- c("ont_epi2me_q28", "lc24_apk", "hifi_revio_pbmay24", "sprq_revio_30hr", "element_ultraq_jun2024", "element_avitistddip", "illumina_2x250", "illumina_googlepcrfree")
#doubleplatformlabels <- c("ONT Q28", "ONT APK", "HiFi Revio", "Revio SPRQ", "Element UltraQ", "Element Aviti", "Illumina Old", "Ill PCR Free")
#doubleplatformcolors <- c("#6699CC", "#6699CC", "#01541F", "#01541F", "#332288", "#332288", "#661100", "#661100")

doublereadsetnames <- c("ont_epi2me_q28", "lc24_apk", "hifi_revio_pbmay24", "sprq_revio_30hr", "element_ultraq_jun2024", "element_avitistddip", "illumina_googlepcrplus", "illumina_googlepcrfree")
doubleplatformlabels <- c("ONT Q28", "ONT APK", "HiFi Revio", "Revio SPRQ", "Element UltraQ", "Element Aviti", "Ill PCR Plus", "Ill PCR Free")
doubleplatformcolors <- c("#6699CC", "#6699CC", "#01541F", "#01541F", "#332288", "#332288", "#661100", "#661100")

#readsetnames <- c("ont_epi2me_q28", "lc24_apk", "hifi_revio_pbmay24", "NIST_onso_2024Q1", "element_ultraq_jun2024", "illumina_2x250")
#readsetnames <- c("ont_epi2me_q28", "lc24_apk", "hifi_revio_pbmay24", "element_ultraq_jun2024", "illumina_2x250")
#readsetnames <- c("ont_epi2me_q28", "lc24_apk", "hifi_revio_pbmay24", "element_ultraq_jun2024", "illumina_googlepcrfree_ds0.1")
#platformlabels <- c("ONT Q28", "ONT APK", "HiFi Revio", "Onso", "Element", "Illumina")
#platformlabels <- c("ONT Q28", "ONT APK", "HiFi Revio", "Element", "Illumina")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
readplatformcolors <- c("#6699CC", "#01541F", "#332288", "#661100", "#5D5D5D")
platformlinetype <- c(1, 2, 3, 4, 5, 6)
platformpchvals <- c(0, 1, 2, 5, 6, 8)
filledplatformpchvals <- c(15, 16, 17, 18, 23, 25, 8)
multiplevector <- c(1.0, 1.0, 1.0, 1.4, 1.0, 1.0, 1.0, 1.4, 1.4, 1.4, 1.4, 1.9)

title <- ""

pdf("ReadQVAccuracy.pdf")
read_qv_plot(readsetnames, platformlabels)
dev.off()

pdf("ArangsReadQVAccuracy.pdf")
read_qv_plot(arangsreadsetnames, arangsplatformlabels)
dev.off()

### Make mononucleotide accuracy plot ###

pdf("ReadMononucAccuracy.pdf")  
read_mononucacc_plot(readsetnames, platformlabels)
dev.off()

pdf("ArangsReadMononucAccuracy.pdf")  
read_mononucacc_plot(arangsreadsetnames, arangsplatformlabels)
dev.off()

pdf("ReadMononucErrorRate.pdf")  
read_mononucerror_plot(readsetnames, platformlabels)
dev.off()

pdf("ArangsReadMononucErrorRate.pdf")  
read_mononucerror_plot(arangsreadsetnames, arangsplatformlabels)
dev.off()

pdf("ReadMononucQVScores.pdf", width=7, height=7)  
read_mononucqvscore_plot(readsetnames, platformlabels, ymax=25, legendypos=23, filledpoints=TRUE)
dev.off()

pdf("ReadMononucQVScoresLines.pdf", width=7, height=7)  
read_mononucqvscore_plot(readsetnames, platformlabels, ymax=25, legendypos=23, plotlines=TRUE, linetype=2, errorbars=TRUE)
dev.off()

# Plot substitution rate-by-type histogram

pdf("ReadSubstitutionRates.pdf", width=7, height=7)
read_substitutions_plot(readsetnames, platformlabels, outputdir="Figure4Output", legend=TRUE, shiftnumbers=0.06)
dev.off()

pdf("ArangsReadSubstitutionRates.pdf", width=7, height=7)
read_substitutions_plot(arangsreadsetnames, arangsplatformlabels, outputdir="Figure4Output", legend=TRUE)
dev.off()

### Indel error plot ###

pdf("ReadIndelRates.pdf", width=7, height=7)  
read_indels_plot(readsetnames, platformlabels, shiftnumbers=0.07)
dev.off()

pdf("ArangsReadIndelRates.pdf", width=7, height=7)  
read_indels_plot(arangsreadsetnames, arangsplatformlabels)
dev.off()

### Plot the full figure together:
pdf("Figure4Output/Figure4Multiplot.pdf", width=7, height=7)
par(mfrow=c(2,3))
read_mononucqvscore_plot(readsetnames, platformlabels, strtype='mononuc', xlabel=c("Run length"), ymax=25)
#read_mononucqvscore_plot(readsetnames, platformlabels, strtype='dinuc', ymax=25, minlength=25, xlabel=c("Run length"), plottitle='Accuracy of dinucleotide runs')
#read_mononucqvscore_plot(readsetnames, platformlabels, strtype='trinuc', ymax=25, minlength=25, xlabel=c("Run length"), plottitle='Accuracy trinucleotide runs')
read_qv_plot(readsetnames, platformlabels)
read_substitutions_plot(readsetnames, platformlabels, outputdir="Figure4Output")
read_indels_plot(readsetnames, platformlabels)
dev.off()

# compare Illumina old and new
illuminareadsetnames <- c("illumina_2x250", "illumina_googlepcrplus", "illumina_googlepcrfree", "NIST_onso_2024Q1")
illuminaplatformlabels <- c("2x250 (2016)", "Google PCR Plus", "Google PCR Free", "NIST/PB Onso (2024)")

pdf("IlluminaSubsIndelComparison.pdf", width=20, height=11)
par(mfrow=c(1,2))
read_substitutions_plot(illuminareadsetnames, illuminaplatformlabels, outputdir="Figure4Output", legend=TRUE, legendposx=7)
read_indels_plot(illuminareadsetnames, illuminaplatformlabels, legendposx=7, legendposy=25)
dev.off()

pdf("IlluminaMononucAccuracy.pdf", width=11, height=11)
read_mononucqvscore_plot(illuminareadsetnames, illuminaplatformlabels, errorbars=TRUE, plottitle='Illumina Accuracy of Mononucleotide Runs')
dev.off()

# compare Illumina old and new
#shortreadsetnames <- c("illumina_2x250", "illumina_googlepcrplus_ds0.1", "illumina_googlepcrfree_ds0.1", "NIST_onso_2024Q1", "element_ultraq_jun2024")
#shortreadsetnames <- c("illumina_2x250", "illumina_googlepcrplus", "illumina_googlepcrfree", "NIST_onso_2024Q1", "element_ultraq_jun2024")
shortreadsetnames <- c("illumina_2x250", "illumina_googlepcrplus", "illumina_googlepcrfree", "element_avitistddip", "element_ultraq_jun2024")
#shortplatformlabels <- c("Illumina 2x250 (2016)", "Illumina PCR Plus (2020)", "Illumina PCR Free (2020)", "NIST/PB Onso (2024)", "Element Ultraq (2024)")
shortplatformlabels <- c("Illumina 2x250 (2016)", "Illumina PCR Plus (2020)", "Illumina PCR Free (2020)", "Element Aviti Std", "Element Ultraq (2024)")
shortreadcolors <- c("#88CCEE", "#AA4499", "#661100", "#44AA99", "#332288")


pdf("ShortReadSubsIndelComparison.pdf", width=23, height=11)
par(mfrow=c(1,2))
read_substitutions_plot(shortreadsetnames, shortplatformlabels, platformcolors=shortreadcolors, outputdir="Figure4Output", legend=TRUE, legendcex=1.2, titlecex=1.3, legendposx=10, legendposy=3700)
read_indels_plot(shortreadsetnames, shortplatformlabels, platformcolors=shortreadcolors, legendcex=1.2, titlecex=1.3, legendposx=10, legendposy=28)
dev.off()

pdf("ShortReadMononucAccuracy.pdf", width=11, height=11)
#read_mononucqvscore_plot(shortreadsetnames, shortplatformlabels, errorbars=TRUE, plottitle='Short Read Accuracy of Mononucleotide Runs', titlecex=1.3, legendcex=1.2)
read_mononucqvscore_plot(shortreadsetnames, shortplatformlabels, errorbars=TRUE, plottitle='', legendcex=1.2)
dev.off()

#pdf("CumulativeBinnedCoverage.pdf", width=7, height=7)
#read_cum_coverage_plot(readsetnames, platformlabels, subdir="binnedcov")
#dev.off()

pdf("CumulativeBinnedCoverage.pdf", width=7, height=7)
read_cum_coverage_plot(readsetnames, platformlabels, subdir="bincovsinglefile")
dev.off()

pdf("CumulativeBinnedCoverageWithPCRPlus.pdf", width=7, height=7)
read_cum_coverage_plot(readsetnameswithpcrplus, platformlabelswithpcrplus, subdir="binnedcov")
dev.off()

read_cum_coverage_plot(readsetnames, platformlabels, subdir = "binnedcovnocensats")

extremekmerdf <- read_extreme_kmer_counts(paste0(outdir, "/k40_kmer_ratios.txt"))
plot_readset_relative_ratios(extremekmerdf, doublereadsetnames, barcolors=doubleplatformcolors, legend=TRUE, legendlabels=doubleplatformlabels)

extremekmernocensatdf <- read_extreme_kmer_counts(paste0(outdir, "/k40nocensat_kmer_ratios.txt"))
pdf("ReadDinuc40merNoCensatCoverage.pdf", width=7, height=7)
plot_readset_relative_ratios(extremekmernocensatdf, readsetnames, barcolors=readplatformcolors, legend=TRUE, legendlabels=platformlabels)
dev.off()

pdf("ReadDinuc40merNoCensatWithPCRPlusCoverage.pdf", width=7, height=7)
plot_readset_relative_ratios(extremekmernocensatdf, readsetnameswithpcrplus, barcolors=readplatformcolors, legend=TRUE, legendlabels=platformlabelswithpcrplus)
dev.off()

pdf("ReadDinuc40merCoverage.pdf", width=7, height=7)
plot_readset_relative_ratios(extremekmerdf, readsetnames, barcolors=readplatformcolors, legend=TRUE, legendlabels=platformlabels)
dev.off()

# 3-mer composition in aligned read sequence compared to benchmark 3-mer composition 
threemer_df <- read_threemer_counts(paste0(outdir, "/k3_kmer_ratios.txt"))
pdf("Read3merCoverage.pdf", width=7, height=7)
plot_readset_threemer_percentages(threemer_df, readsetnames, barcolors=readplatformcolors, legend=TRUE, legendlabels=platformlabels)
dev.off()
pdf("Read3merCoverage10threemers.pdf", width=7, height=7)
plot_readset_threemer_percentages(threemer_df, readsetnames, barcolors=readplatformcolors, legend=TRUE, legendlabels=platformlabels, threemerlimit=10)
dev.off()

# 3-mer composition in aligned read sequence, but exclude censats:
threemer_nocensat_df <- read_threemer_counts(paste0(outdir, "/k3nocensat_kmer_ratios.txt"))
pdf("Read3merCoverage.nocensat.pdf", width=7, height=7)
plot_readset_threemer_percentages(threemer_nocensat_df, readsetnames, barcolors=readplatformcolors, legend=TRUE, legendlabels=platformlabels, title="3-mer coverage in read sets excluding censat")
dev.off()
pdf("Read3merCoverage10threemers.nocensat.pdf", width=7, height=7)
plot_readset_threemer_percentages(threemer_nocensat_df, readsetnames, barcolors=readplatformcolors, legend=TRUE, legendlabels=platformlabels, threemerlimit=10, title="3-mer coverage in read sets excluding censat")
dev.off()

# Trying to figure out Illumina read coverage bias:
plot_readset_threemer_percentages(threemer_df, doublereadsetnames, barcolors=doubleplatformcolors, legend=TRUE, legendlabels=doubleplatformlabels, threemerlimit=10, ymin=-0.2, ymax=0.2)
plot_readset_threemer_percentages(threemer_df, plusreadsetnames, barcolors=readplatformcolors, legend=TRUE, legendlabels=plusplatformlabels, threemerlimit=10)

plot_readset_threemer_percentages(threemer_nocensat_df, readsetnames, barcolors=readplatformcolors, legend=TRUE, legendlabels=platformlabels, threemerlimit=64, title="3-mer coverage in read sets excluding censat")

pdf("Bin100CoverageVsGC.pdf", width=7, height=7)
read_bin100_gccov_plot(readsetnames[1:4], platformlabels[1:4], readplatformcolors[1:4], subdir="bin100nocensat", errorbars=TRUE, ymax=50, xpercbin=5, cexmultipliers=c(0.8, 0.8, 0.7, 1.1), pchvals=filledplatformpchvals)
dev.off()

pdf("Bin100CoverageDistribution.pdf", width=7, height=7)
read_bin_covdist_plot(readsetnames[1:4], platformlabels[1:4], readplatformcolors[1:4], subdir="bincovsinglefile")
dev.off()
