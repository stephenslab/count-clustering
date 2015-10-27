

## GTEX V6 structure plots

setwd('/Users/kushal/Documents/count-clustering/src') 

docweights <- as.matrix(read.table("../data/GTEX_V6/admix_out_GTEX_V6/omega_cis_genes_0_1_2.txt"));
K=dim(docweights)[2];
color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");


samples_id=read.table("/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/samples_id.txt")[,3];

docweights_per_tissue_mean <- apply(docweights,2,function(x) tapply(x,samples_id,mean));
ordering=heatmap(docweights_per_tissue_mean)$rowInd;

unique_samples_id_ordered = unique(samples_id)[ordering];

clus_ordered =unlist(lapply(1:53, function(x) which(samples_id == unique_samples_id_ordered[x])));


docweights_ordering = docweights[clus_ordered,];
samples_id_ordered = samples_id[clus_ordered];
#indices1 <- which(samples_id_ordered =="Artery - Coronary");
#indices2 <- 1:7667;
#indices3 <- 7668:8227;
#indices4 <- 8361:8555;

#indices <- c(indices2, indices1, indices3, indices4);
#samples_id_ordered <- samples_id_ordered[indices];
#docweights_ordering <- docweights_ordering[indices,]

png(filename=paste0('../plots/GTEX_V6_thin_',0,'.png'),width=700,height=300)
par(mar=c(15,2,2,1))
barplot(t(docweights_ordering),col=color[1:K],axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=0.7,cex.main=1)

labels = match(unique(samples_id_ordered), samples_id_ordered);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(samples_id_ordered),las=2, cex.axis=0.8);
dev.off()