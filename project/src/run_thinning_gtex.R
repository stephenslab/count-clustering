
## The gtex structure plots under thinning

setwd('/Users/kushal/Documents/count-clustering/src') 
for(run in 1:3)
{
  docweights <- as.matrix(read.table(paste0("../internal_data/multiple_runs/run_",run,"_omega.txt")));
  K=dim(docweights)[2];
  color=c("red","blue","cornflowerblue","black","cyan","darkblue",
          "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
          "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");
  
  
  samples_id=read.table("/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/samples_id.txt")[,3];
  
  docweights_per_tissue_mean <- apply(docweights,2,function(x) tapply(x,samples_id,mean));
 # ordering=heatmap(docweights_per_tissue_mean)$rowInd;
  
  unique_samples_id_ordered = unique(samples_id)[ordering];
  
  clus_ordered =unlist(lapply(1:53, function(x) which(samples_id == unique_samples_id_ordered[x])));
  
  
  docweights_ordering = docweights[clus_ordered,];
  samples_id_ordered = samples_id[clus_ordered];
  png(filename=paste0('../plots/multiple_run_',run,'.png'),width=700,height=300)
  par(mar=c(14,2.5,1,1))
  barplot(t(docweights_ordering),col=color[1:K],axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=1.2,cex.main=1.4)
  
  labels = match(unique(samples_id_ordered), samples_id_ordered);
  abline(v=labels)
  
  labels_low=labels;
  labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
  mid_point=labels_low +0.5*(labels_up-labels_low);
  
  axis(1,at=mid_point, unique(samples_id_ordered),las=2);
  dev.off()
  
  
}

