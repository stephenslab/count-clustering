

##  Structure plots for the paper

setwd('/Users/kushal/Documents/count-clustering/src')
docweights_thinned_version_2 <- as.matrix(read.table("../internal_data/gtex_thinned/topics_omega_clus_12_version_1.txt"));
K=12;
#docweights <- docweights_thinned_version_2;
#omega <- as.matrix(read.table("../internal_data/gtex_thinned/topics_omega_clus_12_version_2.txt"));


color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");


samples_id=read.table("/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/samples_id.txt")[,3];

docweights_per_tissue_mean <- apply(docweights_thinned_version_2,2,function(x) tapply(x,samples_id,mean));
#ordering=heatmap(docweights_per_tissue_mean)$rowInd;

unique_samples_id_ordered = unique(samples_id)[ordering];

clus_ordered =unlist(lapply(1:53, function(x) which(samples_id == unique_samples_id_ordered[x])));


docweights_ordering = docweights_thinned_version_2[clus_ordered,];
samples_id_ordered = samples_id[clus_ordered];
png(filename='../plots/gtex_thinned_clus_12_version_1.png',width=700,height=300)
par(mar=c(14,2.5,1,1))
barplot(t(docweights_ordering),col=color[1:K],axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=1.2,cex.main=1.4)

labels = match(unique(samples_id_ordered), samples_id_ordered);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(samples_id_ordered),las=2);
dev.off()

