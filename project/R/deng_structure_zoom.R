
## Deng micro structure plot (Zoom in)

rm(list=ls())
library(data.table)
#install_github('kkdey/maptpx') 
library(maptpx)
library(CountClust)
library(data.table)

reads <- data.frame(fread('../external_data/Deng_Data/Deng_cell_data.txt'),row.names=1);
files <- list.files("../external_data/Deng_Data/Deng_files/");
cell_meta <- unlist(lapply(files, function(x) strsplit(x,"_")[[1]][2]));
embryo_id <- unlist(lapply(files, function(x) strsplit(strsplit(x,"_")[[1]][3],"-")[[1]][1]));
cell_meta[grep("zy",cell_meta)]="zy";
cell_meta[grep("smartseq2", files)]="8cell_nd";
cell_meta[grep("8cell_2pooled", files)]="8cell_nd";
cell_meta[grep("8cell_split", files)]="8cell_nd";
cell_meta[grep("16cell_2pooled", files)]="16cell_nd";
cell_meta[grep("16cell_split", files)]="16cell_nd";
indices_not_reqd <- which(cell_meta=="BXC"   | cell_meta=="C57twocell" | cell_meta=="fibroblast" | cell_meta =="8cell_nd" | cell_meta == "16cell_nd");
cell_meta <- cell_meta[-indices_not_reqd];
embryo_id <- embryo_id[-indices_not_reqd];
embryo_id[which(embryo_id == "expression.txt")]="."
cell_embryo <- paste0(cell_meta,"_", embryo_id);
reads <- reads[,-indices_not_reqd];
cell_meta_unique <- c("zy","early2cell","mid2cell","late2cell","4cell","8cell","16cell","earlyblast","midblast","lateblast") ;
order_of_development <- order(match(cell_meta,cell_meta_unique))
reads <- reads[,order_of_development];
cell_meta <- cell_meta[order_of_development]
cell_embryo <- cell_embryo[order_of_development];
colnames(reads) <- cell_meta;


samp_metadata <- cbind.data.frame(cell_meta, cell_embryo);
counts <- t(reads);
colnames(samp_metadata) <- c("dev_phase", "dev.phase_embryo");

nclus <- 6

if(file.exists("../rdas/deng_topic_fit.rda")) {
  deng_topics <- get(load("../rdas/deng_topic_fit.rda"));
} else {
  StructureObj(as.matrix(counts),nclus_vec,samp_metadata = samp_metadata, tol=10, batch_lab = NULL, path_rda="../rdas/deng_topic_fit.rda",partition=c('TRUE'),path_struct = "../plots/deng_structure");
  deng_topics <- get(load("../rdas/deng_topic_fit.rda"));
}

index_of_zoom <- which(cell_meta=="4cell" | cell_meta=="8cell" | cell_meta=="16cell");

omega <- deng_topics[[5]]$omega;
omega <- omega[index_of_zoom,];

samp_metadata <- samp_metadata[index_of_zoom,];
colnames(samp_metadata) <- c("dev_phase_zoom", "dev.phase_embryo_zoom");

docweights <- omega
num <- 2
metadata <- samp_metadata[,num];

path_struct = "../plots/deng_structure"

if(!dir.exists(paste0(path_struct,"/clus_",dim(docweights)[2])))
    dir.create(paste0(path_struct,"/clus_",dim(docweights)[2]))
  
  control.default <- list(struct.width=600, struct.height=400, cex.axis=1, cex.main=1.5, las=2, lwd=3,mar.bottom =14, mar.left=2, mar.top=2, mar.right=2,color=2:(dim(docweights)[2]+1));
  control <- control.default;
  struct.width <- control$struct.width;
  struct.height <- control$struct.height;
  cex.axis <- control$cex.axis;
  cex.main <- control$cex.main;
  las <- control$las;
  lwd <- control$lwd;
  mar.bottom <- control$mar.bottom;
  mar.left <- control$mar.left;
  mar.top <- control$mar.top;
  mar.right <- control$mar.right;
  color <- control$color;
  
  png(filename=paste0(path_struct,'/clus_',nclus,'/struct_clus_',nclus,'_',colnames(samp_metadata)[num],'.png'),width=struct.width, height=struct.height);
  par(mar=c(mar.bottom,mar.left, mar.top,mar.right))
  barplot(t(docweights),col=color,axisnames=F,space=0,border=NA,
          las=las,ylim=c(0,1),ylab="admix prop", xlab=paste0(colnames(samp_metadata)[num]),
          cex.axis=cex.axis,cex.main=cex.main);
  labels = match(unique(metadata), metadata);
  abline(v=labels-1, lty=1, lwd=lwd)
  
  labels_low=labels-1;
  labels_up=c(labels_low[2:length(labels_low)],dim(docweights)[1]);
  mid_point <- labels_low +0.5*(labels_up-labels_low);
  axis(1,at=mid_point, unique(metadata),las=las,cex.axis=cex.axis,lwd=lwd);
  dev.off()

