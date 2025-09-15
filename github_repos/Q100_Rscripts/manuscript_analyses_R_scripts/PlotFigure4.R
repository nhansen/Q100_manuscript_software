setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/PaperFigures/Figures")

library(colorspace)
library(Hmisc)
library(plotrix)

source("/Users/nhansen/OneDrive/HG002_diploid_benchmark/Q100_Rscripts/manuscript_analyses_R_scripts/AssemblyBenchComparisonPlotFunctions.R")

################
### FIGURE 3 ###
### Adam's description: ###
### Historical assemblies evaluated against the HG002 genome benchmark. ###
### [Comparison of historical assemblies against the benchmark to include ###
### things like Ash1, HPRCv1, new Verkko, new Hifiasm, and whatever older ###
### HG002 assemblies we can find in GenBank to show progress towards ###
### personalized genomes.] ###
################

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
assemblycolors <- c("#44AA99", "#332288", "#882255", "#888888", "#661100","#999933","#88CCEE", "#CC6677", "#DDCC77")
qvmethodcolors <- c("#DDCC77", "#661100", "#6699CC")
assemblynames <- c("ash1v2", "hprc_year1", "hprc_year2_polished", "hifi_q28_trio_hic")
assemblylabels <- c("Ash1v2 2020", "Yr1 HPRC 2023", "Yr2_polished HPRC 2024", "Revio/Q28 2024")
fiveassemblynames <- c("ash1v2", "hifiasm_2021", "hprc_year1", "hprc_year2_polished", "hifi_q28_trio_hic")
fiveassemblylabels <- c("Ash1v2 2020", "Hifiasm 2021", "Yr1 HPRC 2023", "Yr2 HPRC 2024", "Revio/Q28 2024")
fivewithv2trionames <- c("ash1v2", "hifiasm_2021", "hprc_year1", "hprc_year2_polished", "v2_trio")
fivewithv2triolabels <- c("Ash1v2 2020", "Hifiasm 2021", "Yr1 HPRC 2023", "Yr2 HPRC 2024", "Verkko2 Trio 2025")
ksassemblynames <- c("ash1v2", "hg002v0.1", "hifiasm_2021", "hprc_year1", "hprc_curated", "hprc_year2_polished", "lc24_medaka_6b4_test", "v2_trio", "hprc_year2_v2_thic")
ksassemblylabels <- c("Ash1v2 2020", "HG002v0.1", "Hifiasm 2021", "Yr1 HPRC 2023", "Jarvis HPRC Curated", "HPRC Release2 HiFiasm", "LC24 ONT Medaka", "Verkko2 Trio 2025", "HPRC Release2 Verkko")
finalassemblynames <- c("ash1v2", "hifiasm_2021", "hprc_year1", "verkko_2023", "lc24_medaka_6b4")
finalassemblylabels <- c("Ash1 2020", "Hifiasm 2021", "HPRCv1 2022", "Verkko 2023", "ONT LC24 2024")
finalshortassemblylabels <- c("Ash1v2", "Hifiasm", "HPRCv1", "Verkko v1", "ONT/LC24")
# for serge's friday talk Jun 20 2025:
#skassemblynames <- c("verkko_2023", "v2_trio", "lc24_nopolish", "lc24_medaka_6b4_test", "hifiasm_ontonly_2025")
#skassemblylabels <- c("Verkko 2023", "Verkko2 Trio 2025", "LC24 Unpolished", "LC24 ONT Medaka", "Hifiasm ONTonly 2025")
#skshortassemblylabels <- c("Verkko", "Verkko2 Trio", "LC24 Unpolished", "LC24 Medaka", "Hifiasm ONTonly")
# for serge's south korea talk:
#skassemblynames <- c("verkko_2023", "lc24_nopolish", "lc24_medaka_6b4_test", "hifiasm_ontonly_2025", "verkko_hg002_sequel")
#skassemblylabels <- c("Verkko 2023", "LC24 Unpolished", "LC24 ONT Medaka", "Hifiasm ONTonly 2025", "Verkko 2.3")
#skshortassemblylabels <- c("Verkko", "LC24 Unpolished", "LC24 Medaka", "Hifiasm ONTonly", "Verkko 2.3")

# examining continuity of v3_issue-masked hg002v1.1 (vs. other assemblies from figure 3 in the paper):
maskedassemblynames <- c("ash1v2", "hifiasm_2021", "hprc_year1", "verkko_2023", "lc24_medaka_6b4", "v1.1.v3issuesmasked")
maskedassemblylabels <- c("Ash1 2020", "Hifiasm 2021", "HPRCv1 2022", "Verkko 2023", "ONT LC24 2024", "Issue-masked v1.1")
maskedshortassemblylabels <- c("Ash1v2", "Hifiasm", "HPRCv1", "Verkko v1", "ONT/LC24", "Issue-masked")

skassemblynames <- c("verkko_24hr_hic_2.2", "verkko_24hr_hic_2.3")
skassemblylabels <- c("Verkko 2.2", "Verkko 2.3")
skshortassemblylabels <- c("Verkko 2.2", "Verkko 2.3")

test1assemblynames <- c("lc24_nopolish")
test1assemblylabels <- c("LC24 Unpolished")

### Plot the final figure in separate pdfs so we can break the axis in the substitutions plot:
finalassemblysizefiles <- sapply(finalassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
finalmnstatsfiles <- sapply(finalassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})
finalsubstitutionstatsfiles <- sapply(finalassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})
finalindelstatsfiles <- sapply(finalassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})
assemblyqvfile <- "Figure3Output/gqc_assembly_qvs.txt"
yakqvfile <- "Figure3Output/yak_sprq_element_std_hybrid_qvs.txt"
merqqvfile <- "Figure3Output/merqury_qvs.txt"
phaseswitchfile <- "Figure3Output/switchrates.txt"
missingfile <- "Figure3Output/missingbases.txt"

pdf("Figure3NGAx.pdf", width=7, height=7)
assembly_ngax_plot(finalassemblysizefiles, assemblylabels=finalassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
dev.off()
pdf("Figure3Subs.pdf", width=7, height=7)
# Substitution panel values added manually in adobe illustrator:
figure3subsvalues <- assembly_substitutions_plot(finalsubstitutionstatsfiles, finalshortassemblylabels, cexnames=0.8, legendlabels=finalassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies", ybreak=c(30, 135), ymax=158)
dev.off()
pdf("Figure3Indels.pdf", width=7, height=7)
assembly_indels_plot(finalindelstatsfiles, finalshortassemblylabels, titleval="Indel discrepancies in assemblies", cexnames=0.8, legendlabels=finalassemblylabels, legendypos=20.0)
dev.off()
pdf("Figure3HPRuns.pdf", width=7, height=7)
assembly_mononucqv_plot(finalmnstatsfiles, finalassemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
dev.off()
pdf("Figure3QVs.pdf", width=7, height=7)
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, finalassemblynames, finalshortassemblylabels)
dev.off()
pdf("Figure3SwitchMissing.pdf", width=7, height=7)
par(mfrow=c(2,1))
assembly_switchrate_plot(finalassemblynames, finalassemblylabels, phaseswitchfile)
assembly_missingness_plot(finalassemblynames, finalassemblylabels, missingfile)
dev.off()

# NGA plots
assemblysizefiles <- sapply(assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
fiveassemblysizefiles <- sapply(fiveassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
fivewithv2triosizefiles <- sapply(fivewithv2trionames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
ksassemblysizefiles <- sapply(ksassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
skassemblysizefiles <- sapply(skassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
maskedassemblysizefiles <- sapply(maskedassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})

# Make NGAx plots:

#pdf("Figure3Output/OriginalFourNGAxPlot.pdf")
#assembly_ngax_plot(assemblysizefiles, assemblylabels=assemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
#dev.off()

#pdf("Figure3Output/FiveAssemblyNGAxPlot.pdf")
#assembly_ngax_plot(fiveassemblysizefiles, assemblylabels=fiveassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
#dev.off()

#pdf("Figure3Output/FiveWithVerkko2TrioNGAxPlot.pdf")
#assembly_ngax_plot(fivewithv2triosizefiles, assemblylabels=fiveassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
#dev.off()

#pdf("Figure3Output/KitchenSinkNGAxPlot.pdf")
#assembly_ngax_plot(ksassemblysizefiles, assemblylabels=ksassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
#dev.off()

#pdf("Figure3Output/AdamsNGAxPlot.pdf")
#assembly_ngax_plot(finalassemblysizefiles, assemblylabels=finalassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
#dev.off()

pdf("Figure3Output/V1.1IssueMaskedNGAxPlot.pdf")
assembly_ngax_plot(maskedassemblysizefiles, assemblylabels=maskedassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
dev.off()

# Mononucleotide accuracy:
mnstatsfiles <- sapply(assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})
ksmnstatsfiles <- sapply(ksassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})
skmnstatsfiles <- sapply(skassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})

# Make mononucleotide accuracy/QV plots

pdf("Figure3Output/AssemblyMononucAccuracy.pdf")  
assembly_mononucacc_plot(mnstatsfiles, assemblylabels)
dev.off()

pdf("Figure3Output/KitchenSinkMononucAccuracy.pdf")  
assembly_mononucacc_plot(ksmnstatsfiles, ksassemblylabels)
dev.off()

pdf("Figure3Output/KitchenSinkAssemblyMononucErrorQVs.pdf")  
assembly_mononucqv_plot(ksmnstatsfiles, ksassemblylabels)
dev.off()

pdf("Figure3Output/AdamsMononucErrorQVs.pdf")  
assembly_mononucqv_plot(finalmnstatsfiles, finalassemblylabels)
dev.off()

# Indel errors:
indelstatsfiles <- sapply(assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})
ksindelstatsfiles <- sapply(ksassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})
skindelstatsfiles <- sapply(skassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})

# Make indels plots:

pdf("Figure3Output/IndelRates.pdf", width=9, height=8)  
assembly_indels_plot(indelstatsfiles, assemblylabels, titleval="Indel discrepancies in assemblies")
dev.off()

pdf("Figure3Output/KitchenSinkIndelRates.pdf", width=9, height=8)  
assembly_indels_plot(ksindelstatsfiles, ksassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=18.0)
dev.off()

pdf("Figure3Output/AdamsRates.pdf", width=9, height=8)  
assembly_indels_plot(finalindelstatsfiles, finalassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=18.0)
dev.off()


### Phasing errors plot

#pdf("Figure3Output/PhaseSwitchRates.pdf")  
#assembly_switchrate_plot(assemblynames, assemblylabels, phaseswitchfile)
#dev.off()

pdf("Figure3Output/KitchenSinkPhaseSwitchRates.pdf")  
assembly_switchrate_plot(ksassemblynames, ksassemblylabels, phaseswitchfile)
dev.off()

pdf("Figure3Output/AdamsPhaseSwitchRates.pdf")  
assembly_switchrate_plot(finalassemblynames, finalassemblylabels, phaseswitchfile)
dev.off()

### Substitution rates

substitutionstatsfiles <- sapply(assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})
kssubstitutionstatsfiles <- sapply(ksassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})
sksubstitutionstatsfiles <- sapply(skassemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})

# Make substitutions plot:

pdf("Figure3Output/SubstitutionRates.pdf", width=9, height=8)  
assembly_substitutions_plot(substitutionstatsfiles, assemblylabels, titleval="Substitution discrepancies in assemblies")
dev.off()

pdf("Figure3Output/KitchenSinkSubstitutionRates.pdf", width=9, height=8)  
assembly_substitutions_plot(kssubstitutionstatsfiles, ksassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
dev.off()

pdf("Figure3Output/AdamsSubstitutionRates.pdf", width=9, height=8)  
assembly_substitutions_plot(finalsubstitutionstatsfiles, finalassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
dev.off()

### Quality value (QV) scores for different assemblies by different methods

#pdf("Figure3Output/QVMeasures.pdf", width=13, height=8.5)  
#assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, assemblynames, assemblylabels)
#dev.off()

pdf("Figure3Output/KitchenSinkQVMeasures.pdf", width=13, height=8.5)  
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, ksassemblynames, ksassemblylabels)
dev.off()

pdf("Figure3Output/AdamsQVMeasures.pdf", width=13, height=8.5)  
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, finalassemblynames, finalassemblylabels)
dev.off()

#### Plot everything all together:
#pdf("Figure3KitchenSinkMultiplot.pdf", width=13, height=9)
#par(mfrow=c(2,3))
#assembly_ngax_plot(ksassemblysizefiles, assemblylabels=ksassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
#assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, ksassemblynames, ksassemblylabels)
#assembly_mononucqv_plot(ksmnstatsfiles, ksassemblylabels)
#assembly_substitutions_plot(kssubstitutionstatsfiles, ksassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
#assembly_indels_plot(ksindelstatsfiles, ksassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=18.0)
#assembly_switchrate_plot(ksassemblynames, ksassemblylabels, phaseswitchfile)
#dev.off()
#### Plot the final figure in separate pdfs so we can break the axis in the substitutions plot:
#### Plot them all together:
#pdf("Figure3AdamsMultiplot.pdf", width=10, height=7)
#par(mfrow=c(2,3))
#assembly_ngax_plot(adamsassemblysizefiles, assemblylabels=adamsassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
#assembly_substitutions_plot(adamssubstitutionstatsfiles, adamsshortassemblylabels, cexnames=0.8, legendlabels=adamsassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
#assembly_indels_plot(adamsindelstatsfiles, adamsshortassemblylabels, titleval="Indel discrepancies in assemblies", cexnames=0.8, legendlabels=adamsassemblylabels, legendypos=20.0)
##assembly_mononucqv_plot(adamsmnstatsfiles, adamsassemblylabels, pointcex=1.0)
#assembly_mononucqv_plot(adamsmnstatsfiles, adamsassemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
#assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, adamsassemblynames, adamsshortassemblylabels)
#assembly_switchrate_plot(adamsassemblynames, adamsassemblylabels, phaseswitchfile)
#dev.off()

#pdf("Figure3AdamsMultiplotLinesBars.pdf", width=16, height=10)
#par(mfrow=c(2,3))
#assembly_ngax_plot(adamsassemblysizefiles, assemblylabels=adamsassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
#assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, adamsassemblynames, adamsassemblylabels)
#assembly_mononucqv_plot(adamsmnstatsfiles, adamsassemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
##assembly_substitutions_plot(adamssubstitutionstatsfiles, adamsassemblylabels, legendypos=10, titleval="Substitution discrepancies in assemblies", ybreak=c(40, 120), ymax=160, spanbreak=c(9,12))
#assembly_substitutions_plot(adamssubstitutionstatsfiles, adamsassemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies" )
#assembly_indels_plot(adamsindelstatsfiles, adamsassemblylabels, titleval="Indel discrepancies in assemblies", legendypos=20.0)
#assembly_switchrate_plot(adamsassemblynames, adamsassemblylabels, phaseswitchfile)
#dev.off()

### Plot them all together:
##pdf("Figure3SergesKitchenSinkMultiplot.pdf", width=13, height=9)
#pdf("SergesSouthKoreaTalkMultiplot.pdf", width=16, height=9)
pdf("SergesBackPocketMultiplot.pdf", width=16, height=9)
par(mfrow=c(2,3))
assembly_ngax_plot(skassemblysizefiles, assemblylabels=skshortassemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
assembly_qv_plot(assemblyqvfile, yakqvfile, merqqvfile, skassemblynames, skassemblylabels)
assembly_mononucqv_plot(skmnstatsfiles, skshortassemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
assembly_substitutions_plot(sksubstitutionstatsfiles, skshortassemblylabels, legendypos=2.65, titleval="Substitution discrepancies in assemblies")
assembly_indels_plot(skindelstatsfiles, skshortassemblylabels, titleval="Indel discrepancies in assemblies", legendxpos=1, legendypos=62.0)
assembly_switchrate_plot(skassemblynames, skshortassemblylabels, phaseswitchfile)
dev.off()

# see how it looks to plot just one assembly's data:
#pdf("Test1AssemblyMultiplot.pdf", width=9, height=9)
#par(mfrow=c(2,2))
#test1assemblysizefiles <- sapply(test1assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
#test1substitutionstatsfiles <- sapply(test1assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})
#test1indelstatsfiles <- sapply(test1assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".indelerrorstats.txt"), sep="", collapse=""); return(file)})
#test1mnstatsfiles <- sapply(test1assemblynames, function(x) {file=paste(c("Figure3Output/", x, ".mononucstats.txt"), sep="", collapse=""); return(file)})
#assembly_ngax_plot(test1assemblysizefiles, assemblylabels=test1assemblylabels, ideal=TRUE, plottitle="NGAx for different assemblies")
#assembly_mononucqv_plot(test1mnstatsfiles, test1assemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
#assembly_substitutions_plot(test1substitutionstatsfiles, test1assemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
#assembly_indels_plot(test1indelstatsfiles, test1assemblylabels, titleval="Indel discrepancies in assemblies", legendypos=62.0)
#dev.off()


