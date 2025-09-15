setwd("/Users/nhansen/OneDrive/HG002_diploid_polishing_validation/polishing_update_slides/plotsforadamstalk")

correctionlength_hist <- read.table("correction_length.hist", header=FALSE, sep="\t")
falsehetlength_hist <- read.table("correction_length.falsehets.hist", header=FALSE, sep="\t")
falsehomlength_hist <- read.table("correction_length.falsehoms.hist", header=FALSE, sep="\t")
hetwitherrorlength_hist <- read.table("correction_length.hetwitherror.hist", header=FALSE, sep="\t")

correctionlength_hist_notail <- correctionlength_hist[correctionlength_hist$V1>-6 & correctionlength_hist$V1<6,]
falsehetlength_hist_notail <- falsehetlength_hist[falsehetlength_hist$V1>-6 & falsehetlength_hist$V1<6,]
falsehomlength_hist_notail <- falsehomlength_hist[falsehomlength_hist$V1>-6 & falsehomlength_hist$V1<6,]
hetwitherrorlength_hist_notail <- hetwitherrorlength_hist[hetwitherrorlength_hist$V1>-6 & hetwitherrorlength_hist$V1<6,]

barplot(correctionlength_hist_notail$V2, names.arg=correctionlength_hist_notail$V1, main="HG002v0.7 correction lengths", xlab="length(corrected)-length(assembly)", ylab="Number of corrections")

barplot(rbind(falsehetlength_hist_notail$V2, falsehomlength_hist_notail$V2, hetwitherrorlength_hist_notail$V2), names.arg=hetwitherrorlength_hist_notail$V1, xlab="Correction length (bases)", ylab="# Corrections", col=c("brown", "blue", "darkgreen"), cex.lab=1.5, cex.axis=1.1, main="Length of corrections applied to HG002v0.7 to create HG002v0.9")
legend("topright", c("FalseHet", "FalseHom", "HetError"), col=c("brown", "blue", "darkgreen"), pch=15, cex=1.4)
