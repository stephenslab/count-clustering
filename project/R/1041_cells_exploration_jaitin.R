

######  In this code, we try to replicate the CAP projection plot code sent

######  by Effi and others and also try to do the batchwise and cluster label 

######   wise distruct plots, and also do the PCA analysis

library(data.table)
Data=data.frame(fread("../external/GSE54006_umitab.txt"));
gene_names=as.matrix(Data[,1]);
Exp_details=read.table(file="D:/Matthew_Stephens_Project/Jaitin_single_cell/Data/jaitin_etal_experimental_design.txt",fill=T,header=T);




#################   replicating the CAP  projection plot  ####################


indices=read.table("D:/Matthew_Stephens_Project/Jaitin_single_cell/Data/indices_1041_cells.txt");
indices=as.vector(as.matrix(indices));

loglikelihood=read.table("D:/Matthew_Stephens_Project/Jaitin_single_cell/Data/GSEloglikelihoodData.txt");
rownames_loglik=rownames(loglikelihood);
bwlist=strsplit(rownames_loglik,"_");
index_loglik = as.numeric(matrix(unlist(bwlist),nrow=2)[2,])


###############  extracting the loglikelihood table for the 1041 cells

mm=match(indices,index_loglik);
reduced_loglikelihood_table=loglikelihood[mm,];

plot_CAP(reduced_loglikelihood_table)
mtext("CAP plot for Jaitin loglikelihood data",side=1,line=3);

gene_names=read.csv("D:/Matthew_Stephens_Project/Jaitin_single_cell/Data/gene_names.csv");
Exp_details=read.table(file="D:/Matthew_Stephens_Project/Jaitin_single_cell/Data/jaitin_etal_experimental_design.txt",fill=T,header=T);
gene_names=gene_names[,2];
ERCC_genes=grep("ERCC",gene_names);
Data_non_ERCC=Data[-ERCC_genes,];
counts =t(as.matrix(Data_non_ERCC[,-1]))

batch_well_ID=colnames(Data_non_ERCC)[-1]
bwlist=strsplit(batch_well_ID,"_");
index_ID = as.numeric(matrix(unlist(bwlist),nrow=2)[2,])  

filterlist=c("M34473","abParts","M13680","Tmsb4x","S100a4","B2m","Atpase6","Rpl23","Rps18","Rpl13","Rps19","H2-Ab1","Rplp1","Rpl4","Rps26","EF437368") ;
fcounts = counts[,-match(filterlist,gene_names[-ERCC_genes])];



dim(fcounts)



fcounts_reduced=fcounts[match(indices,index_ID),];

###  fcounts_reduced is the counts matrix of the 1041 cells for which the CAP plot was made and 
###  the 200091 genes of interest (removing controls)




library(maptpx); library(slam);
Topic_Clus=topics(fcounts_reduced,7,kill=0,tol=30);

docweights=Topic_Clus$omega;
write.table(docweights,"D:/Matthew_Stephens_Project/Jaitin_single_cell/Data/docweights_1041_cells_cluster_7",sep="\t");

rowSums(docweights); # check that all values are equal to 1 or not
cluster_index=array(0,dim(docweights)[1]);
for(i in 1:dim(docweights)[1])
{
	cluster_index[i]=which(docweights[i,]==max(docweights[i,]),arr.ind=TRUE);
}

poster=t(docweights);
Nseed=7;

loglikelihood_topics=matrix(0,dim(docweights)[1],dim(docweights)[2]);

for (j in 1:dim(docweights)[1])
{
	loglikelihood_topics[j,]=reduced_loglikelihood_table[j,1]+
						log(docweights[j,]/docweights[j,1]);
}
loglikelihood_topics=data.frame(loglikelihood_topics)
names(loglikelihood_topics)=names(reduced_loglikelihood_table)
rownames(loglikelihood_topics)=rownames(reduced_loglikelihood_table)

plot_CAP(loglikelihood_topics)

Labels=numeric(0);
rowMax=apply(docweights,1,max);
length_vec=array(0,dim(docweights)[2])
for (idx in 1:dim(docweights)[2])
{
	labels=which(docweights[,idx]==rowMax);
	length_vec[idx]=length(labels);
	Labels=c(Labels,labels);
}
dat_sort=docweights[Labels,];


#################  go to Structurelikeplots.R and run the first section ####

names_reduced=rownames(reduced_loglikelihood_table);
bwlist=strsplit(names_reduced,"_");
batch = as.numeric(matrix(unlist(bwlist),nrow=2)[1,])

#################  go to Structurelieplots.R and run the second section ####


Princomp.analysis=prcomp(docweights,center=TRUE);
eigenPrincomp=Princomp.analysis$sdev;
rotPrincomp=Princomp.analysis$rotation;
barplot(eigenPrincomp,col="red", 
		ylab="eigenvalues of Q matrix",
		xlab="Labels of topics or clusters",xaxt="n");
axis(1,at=0.7+1.2*(0:6),c("PC1","PC2","PC3","PC4","PC5","PC6","PC7"));

projection=t(rotPrincomp[,-dim(rotPrincomp)[2]])%*%t(docweights);

plot(projection[1,],projection[2,],main="FirstPC vs SecondPC",xlab="First PC",
	ylab="Second PC",col=cluster_index);

plot(projection[1,],projection[3,],main="FirstPC vs ThirdPC",xlab="First PC",
	ylab="Third PC",col=cluster_index);


################################  The End  ###############################






