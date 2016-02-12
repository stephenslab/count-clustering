

## plot hierarchical and admixture F for GTEX data
setwd("/Users/kushal/Documents/count-clustering/src/")
hierarchy_F <- read.table("../internal_data/hierarchy_F_thin_0.txt")
admixture_F <- read.table("../internal_data/admixture_F_thin_0.txt")

#hierarchy_F <- read.table("../internal_data/hierarchy_F_clus_5_thin_1.txt")
#admixture_F <- read.table("../internal_data/admixture_F_clus_5_thin_1.txt")

hierarchy_F <- hierarchy_F[1:48,1:48];
admixture_F <- admixture_F[1:48,1:48];

for(i in 1:dim(hierarchy_F)[1])
{ 
  for(j in 1:(i-1))
  {
    hierarchy_F[j,i] <- hierarchy_F[i,j];
    admixture_F[j,i] <- admixture_F[i,j];
  }
}

hierarchy_F[hierarchy_F < 0.8] = 0;
admixture_F[admixture_F < 0.8] = 0;

heatmap(as.matrix(hierarchy_F))
heatmap(as.matrix(admixture_F))

samples_id=read.table("../external_data/GTEX_V6/samples_id.txt");


unique_tissues <- unique(samples_id[,3]);
unique_tissues <- unique_tissues[-which(as.numeric(table(samples_id[,3])) < 50)];

colnames(hierarchy_F) <- unique_tissues
rownames(hierarchy_F) <- unique_tissues
colnames(admixture_F) <- unique_tissues
rownames(admixture_F) <- unique_tissues

library(cape)

png(filename="../plots/hierarchy_F_thin_clus_5_1.png", res=100)
myImagePlot(as.matrix(hierarchy_F))
dev.off()

png(filename="../plots/admixture_F_thin_clus_5_1.png", res=100)
myImagePlot(as.matrix(admixture_F))
dev.off()




