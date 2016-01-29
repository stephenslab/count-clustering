
## Zeisel single cell structure plots 

library(data.table)
library(CountClust)
data=suppressWarnings(suppressMessages(data.frame(fread("../data/expression_mRNA_17-Aug-2014.txt"))));

counts_data=as.matrix(data[-(1:10),-(1:2)]);

counts_data <- apply(counts_data,2,as.numeric)
counts_data <- t(counts_data);

cell_type = unlist(lapply(strsplit(rownames(counts_data),"\\."),function(x) x[1]));

nclus_vec <- c(2,5,7,10);

if(!dir.exists("Structure")) dir.create("Structure")
if(!dir.exists("Structure/batch_uncorrected")) dir.create("Structure/batch_uncorrected")

finer_subtype = data[9,-(1:2)];

for(num in 1:length(nclus_vec))
{
  if(!dir.exists(paste0("Structure/batch_uncorrected/clus_",nclus_vec[num]))) dir.create(paste0("Structure/batch_uncorrected/clus_",nclus_vec[num]))
  omega <- as.matrix(read.table(paste0("Structure/batch_uncorrected/clus_",nclus_vec[num],"/omega_mat.txt")));
  omega <- omega[which(finer_subtype!="(none)"),];
  finer_subtype1 <- as.vector(as.matrix(finer_subtype[which(finer_subtype!="(none)")]));
  cell_type1 <- cell_type[which(finer_subtype!="(none)")];
  samp_metadata <- cbind.data.frame(cell_type1,finer_subtype1); colnames(samp_metadata) <- c("cell_type","cell_subtype");
  obj <- structureObj_omega(omega,samp_metadata = samp_metadata, tol=0.001, 
                            batch_lab = NULL, path=paste0("Structure/batch_uncorrected/clus_",nclus_vec[num]),
                            partition=rep('TRUE',ncol(samp_metadata)),
                            control <- list(struct.width=800, struct.height=300,cex.axis=1,lwd=3,mar.bottom=10,mar.left=2.2));
}
