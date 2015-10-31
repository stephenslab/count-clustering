###  Brain GTEx V6 data

data <- data.frame(fread('../data/GTEX_V6/cis_gene_expression.txt'));
matdata <- data[,-(1:2)];
samples_id=read.table("/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/samples_id.txt")[,3];
brain_indices <- grep("Brain", samples_id);

brain_data <- matdata[,brain_indices];
colnames(brain_data) <- samples_id[brain_indices];
brain_data_frame <- cbind.data.frame(data[,2],brain_data);

write.table(brain_data_frame, '../data/GTEX_V6/cis_gene_expression_brain.txt');

gene_names <- cbind.data.frame(data[,2]);
colnames(gene_names) <- "cis_gene_names"
write.table(gene_names, "../data/GTEX_V6/gene_names_GTEX_V6.txt");

docweights <- as.matrix(read.table('../data/GTEX_V6/admix_out_GTEX_V6/omega_cis_genes_brain.txt'));
K=dim(docweights)[2];
color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");
brain_ids <- samples_id[brain_indices];
ordering <- order(brain_ids);
samples_id_ordered <- brain_ids[ordering];
docweights_ordering <- docweights[ordering,];
png(filename=paste0('../plots/GTEX_V6_brain_thin_',0,'.png'),width=700,height=300)
par(mar=c(14,2,2,1))
barplot(t(docweights_ordering),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=1,cex.main=1)

labels = match(unique(samples_id_ordered), samples_id_ordered);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights_ordering)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(samples_id_ordered),las=2, cex.axis=0.8);
dev.off()

