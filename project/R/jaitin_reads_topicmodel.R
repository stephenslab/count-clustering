
### gtex reads topic model
library(data.table)
library(maptpx)
setwd('/Users/kushal/Documents/count-clustering/src')
data <- t(data.frame(fread('../data/Jaitin_data/GSE54006_umitab.txt'),row.names=1));
gene_names=colnames(data);

Exp_details=read.table(file="../data/Jaitin_data/jaitin_etal_experimental_design.txt",fill=T,header=T);

indices=read.table("../data/Jaitin_data/indices_1041_cells.txt");
indices=as.vector(as.matrix(indices));

gene_names=read.csv("../data/Jaitin_data/gene_names.csv")[,2];
ERCC_genes=grep("ERCC",gene_names);
data_non_ERCC=data[,-ERCC_genes];
counts <- data_non_ERCC

batch_well_ID=rownames(data_non_ERCC)
bwlist=strsplit(batch_well_ID,"_");
index_ID = as.numeric(matrix(unlist(bwlist),nrow=2)[2,])  

filterlist=c("M34473","abParts","M13680","Tmsb4x","S100a4","B2m","Atpase6","Rpl23","Rps18","Rpl13","Rps19","H2-Ab1","Rplp1","Rpl4","Rps26","EF437368") ;
fcounts = counts[,-match(filterlist,gene_names[-ERCC_genes])];

fcounts_reduced=fcounts[match(indices,index_ID),];

Topic_Clus=topics(fcounts_reduced,7,kill=0,tol=0.005);
docweights=Topic_Clus$omega;

exp_labels <- match(indices,Exp_details$index);

Exp_details_reduced <- Exp_details[exp_labels,];

omega <- as.matrix(docweights)

write.table(omega,"../internal_data/jaitin_omega_7.txt");

omega <- as.matrix(read.table("../internal_data/jaitin_omega_7.txt"))

library(CountClust)

if(!dir.exists("../internal_data/Structure")) dir.create("../internal_data/Structure")

samp_metadata <- cbind.data.frame(Exp_details_reduced$sequencing_batch, Exp_details_reduced$amplification_batch);
colnames(samp_metadata) = c("seq.batch", "amp.batch");

library(CountClust)
obj <- StructureObj_omega(omega, samp_metadata = samp_metadata, batch_lab = NULL, path="../internal_data/Structure",partition = c("TRUE","TRUE"),
                          control=list(cex.axis=1));



