setwd("/Users/nhansen/HG002_diploid_benchmark/plots/benchhapmerdistances")

hprc_curated_hapmer_distances <- read.table("hprc_curated_distances_between_hapmers.sort.txt")
names(hprc_curated_hapmer_distances) <- c("distances")
hprc_curated_hist_data <- hist(hprc_curated_hapmer_distances$distances, plot=FALSE)

options(scipen=999)
plot(hprc_curated_hist_data$mids, hprc_curated_hist_data$count, col='darkgreen', log="y", type='h', xlab='Distance between adjacent bench hapmers', ylab='Count', lwd=10, lend=2, main="HPRC curated assembly inter-hap distances", cex.lab=0.8, cex.axis=0.7)

v1.1_hapmer_distances <- read.table("v1.1.k40.distances_between_hapmers.nordnagaps.sort.txt")
names(v1.1_hapmer_distances) <- c("distances")
v1.1_hist_data <- hist(v1.1_hapmer_distances$distances, plot=FALSE)

options(scipen=999)
plot(v1.1_hist_data$mids, v1.1_hist_data$count, col='darkred', log="y", type='h', xlab='Distance between adjacent bench hapmers', ylab='Count', lwd=10, lend=2, main="hg002v1.1 inter-hapmer distances", cex.lab=0.8, cex.axis=0.7)

