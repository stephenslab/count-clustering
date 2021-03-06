

## GTEX V6 structure plots

#setwd('/Users/kushal/Documents/count-clustering/src')

docweights <- as.matrix(read.table("../external_data/GTEX_V6/admix_out_GTEX_V6/omega_cis_genes_0_1_2.txt"));
K=dim(docweights)[2];
color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");


samples_id=read.table("../external_data/GTEX_V6/samples_id.txt")[,3];
samples_id <- as.character(samples_id)
docweights_per_tissue_mean <- apply(docweights,2,function(x) tapply(x,samples_id,mean));
ordering=heatmap(docweights_per_tissue_mean)$rowInd;
unique_samples_id_ordered = unique(samples_id)[ordering];
clus_ordered =unlist(lapply(1:53, function(x) which(samples_id == unique_samples_id_ordered[x])));


samples_id[grep("Nucleus", samples_id)] = "Brain -N. accumbens"
samples_id[grep("Putamen", samples_id)] = "Brain -Putamen"
samples_id[grep("Caudate", samples_id)] = "Brain -Caudate"
samples_id[grep("Gastroe", samples_id)] = "Esophagus -Gastroesophageal Jn."
samples_id[grep("cingulate", samples_id)] = "Brain - Anterior cortex (BA24)."
samples_id[grep("EBV", samples_id)] = "Cells -EBV-lymphocytes"
samples_id[grep("Suprapubic", samples_id)] = "Skin - Unexposed (Suprapubic)"
samples_id[grep("Lower Leg", samples_id)] = "Skin - Exposed (Lower Leg)"


samples_id <- as.factor(samples_id)


docweights_ordering = docweights[clus_ordered,];
samples_id_ordered = samples_id[clus_ordered];
indices1 <- which(samples_id_ordered =="Artery - Coronary");
indices2 <- 1:7667;
indices3 <- 7668:8227;
indices4 <- 8361:8555;

indices <- c(indices2, indices1, indices3, indices4);
samples_id_ordered <- samples_id_ordered[indices];
samples_id_ordered_all <- samples_id_ordered;
docweights_ordering <- docweights_ordering[indices,]

#png(paste0('../plots/fig1_paper/temp1.png'),width=700,height=300)
par(mar=c(15,3,2,1))
barplot(t(docweights_ordering),col=color[1:K],axisnames=F,space=0,border=NA,main=paste("Unthinned data: K=",K),las=1,ylim=c(0,1),cex.axis=0.9,cex.main=1)

labels = match(unique(samples_id_ordered), samples_id_ordered);
abline(v=labels)
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(samples_id_ordered),las=2, cex.axis=1, font.axis=2);
#dev.off()

clus_brain <- clus_ordered[grep("Brain",samples_id_ordered)]


docweights <- as.matrix(read.table('../external_data/GTEX_V6/admix_out_GTEX_V6/omega_cis_genes_brain.txt'));
K=dim(docweights)[2];
color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");
ordering <- clus_brain - min(clus_brain)+1;
samples_id_ordered <- samples_id[clus_brain];
docweights_ordering <- docweights[ordering,];
#png(filename='../plots/fig1_paper/temp2.png',width=700,height=300)
par(mar=c(13,4,2,2))
barplot(t(docweights_ordering),col=2:(K+1),axisnames=F,space=0,border=NA,
        main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=1,cex.main=2)

labels = match(unique(samples_id_ordered), samples_id_ordered);
abline(v=labels,lwd=3)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);
unique_brain_id_ordered <- as.character(unique(samples_id_ordered));
unique_brain_id_ordered[grep("Nucleus",unique(samples_id_ordered))]="Brain -N. accumbens (basal ganglia)"
axis(1,at=mid_point, factor(unique_brain_id_ordered),las=2, cex.axis=0.8, font.axis=2);
#dev.off()

docweights <- as.matrix(read.table("../external_data/GTEX_V6/admix_out_GTEX_V6/omega_cis_genes_0_0001.txt"));
docweights_ordering = docweights[clus_ordered,];
docweights_ordering <- docweights_ordering[indices,]
K=dim(docweights)[2];
samples_id_ordered <- samples_id_ordered_all;
color1 <- color;
color1[15] <- color[14];
color1[14] <- color[13];
color1[10] <- color[10];
color1[3] <- color[4];
color1[7] <- color[7];
color1[8] <- color[9];
color1[11] <- color[12];
color1[12] <- color[11];
color1[1] <- color[1];
color1[2] <- color[2];
color1[4] <- color[3];
color1[5] <- color[8];
color1[6] <- color[6];
color1[9] <- color[5];
color1[13] <- color[15];

#png(paste0('../plots/fig1_paper/temp3.png'),width=700,height=300)
par(mar=c(15,3,2,1))
barplot(t(docweights_ordering),col=color1[1:K],axisnames=F,space=0,border=NA,main=paste("Thinned data: pthin=0.001, K=",K),las=1,ylim=c(0,1),cex.axis=0.9,cex.main=1)

labels = match(unique(samples_id_ordered), samples_id_ordered);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(samples_id_ordered),las=2, cex.axis=1, font.axis=2);
#dev.off()


docweights <- as.matrix(read.table("../external_data/GTEX_V6/admix_out_GTEX_V6/omega_cis_genes_0_001.txt"));
docweights_ordering = docweights[clus_ordered,];
docweights_ordering <- docweights_ordering[indices,]
K=dim(docweights)[2];
samples_id_ordered <- samples_id_ordered_all;
color2 <- color1;
color2[12] <- color1[10];
color2[15] <- color1[15];
color2[14] <- color1[14];
color2[10] <- color1[13];
color2[3] <- color1[3];
color2[7] <- color1[7];
color2[8] <- color1[8];
color2[11] <- color1[12];
color2[1] <- color1[1];
color2[2] <- color1[2];
color2[4] <- color1[4];
color2[5] <- color1[5];
color2[6] <- color1[6];
color2[9] <- color1[9];
color2[13] <- color1[11];

#png(paste0('../plots/fig1_paper/temp3.png'),width=700,height=300)
par(mar=c(15,3,2,1))
barplot(t(docweights_ordering),col=color2[1:K],axisnames=F,space=0,border=NA,main=paste("Thinned data: pthin=0.01, K=",K),las=1,ylim=c(0,1),cex.axis=0.9,cex.main=1)

labels = match(unique(samples_id_ordered), samples_id_ordered);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(samples_id_ordered),las=2, cex.axis=1, font.axis=2);
#dev.off()

docweights <- as.matrix(read.table("../external_data/GTEX_V6/admix_out_GTEX_V6/omega_cis_genes_0_1.txt"));
docweights_ordering = docweights[clus_ordered,];
docweights_ordering <- docweights_ordering[indices,]
K=dim(docweights)[2];
samples_id_ordered <- samples_id_ordered_all;
color3 <- color2;
color3[12] <- color2[9];
color3[15] <- color2[10];
color3[14] <- color2[14];
color3[3] <- color2[3];
color3[7] <- color2[7];
color3[8] <- color2[4];
color3[11] <- color2[12];
color3[1] <- color2[1];
color3[2] <- color2[2];
color3[4] <- color2[5];
color3[5] <- color2[6];
color3[6] <- color2[15];
color3[9] <- color2[8];
color3[10] <- color2[13];
color3[13] <- color2[11];


#png(paste0('../plots/fig1_paper/temp3.png'),width=700,height=300)
par(mar=c(15,3,2,1))
barplot(t(docweights_ordering),col=color3[1:K],axisnames=F,space=0,border=NA,main=paste("Thinned data: pthin=0.001, K=",K),las=1,ylim=c(0,1),cex.axis=0.9,cex.main=1)

labels = match(unique(samples_id_ordered), samples_id_ordered);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(samples_id_ordered),las=2, cex.axis=1, font.axis=2);
#dev.off()
