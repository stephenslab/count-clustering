

#######  Hierarchical vs Admixture: Separation Strength  ###############

library(limma)
library(data.table)

counts <- data.frame(data.table::fread("cis_gene_expression.txt"))
matdata <- counts[,-(1:2)]
voom_out  <- limma::voom(matdata)
voom_weights <- voom_out$weights
voom_data <- t(voom_out$E)

samples_id <- read.table("samples_id.txt")
tissue_labels <- samples_id[,3]
library(DESeq2)

library(edgeR)

cpm_data <- cpm(matdata, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)

tt <- table(tissue_labels)
tissue_names <- names(tt)

tissues_to_consider <- tissue_names[which(as.numeric(table(tissue_labels))>60)]

admix_prop <- matrix(0, length(tissues_to_consider), length(tissues_to_consider))
hierarchy_prop <- matrix(0, length(tissues_to_consider), length(tissues_to_consider))

library(maptpx)
library(slam)

source("topics.R")
source("tpx.R")
source("count.R")


for(m in 2:length(tissues_to_consider)){
  for(n in 1:(m-1)){
    tissue1 <- tissues_to_consider[m];
    tissue2 <- tissues_to_consider[n];
    index1 <- which(tissue_labels==tissue1)
    index2 <- which(tissue_labels==tissue2)
    
    matdata1 <- matdata[,index1[1:50]]
    matdata2 <- matdata[,index2[1:50]]
    pooled_data <- cbind(matdata1, matdata2)
    
    cpm_data1 <- cpm_data[,index1[1:50]]
    cpm_data2 <- cpm_data[,index2[1:50]]
    cpm_data_pooled <- cbind(cpm_data1, cpm_data2)
    
    factor1 <- c(rep(1,50), rep(2,50))
    dd <- dist(t(cpm_data_pooled))
    
    hclusters <- cutree(hclust(dd),2)
    tab <- xtabs(~hclusters + factor1)
    misclass1 <- (min(tab[1,1]+tab[2,2], tab[1,2]+ tab[2,1]))/100;
    cat("Hierarchical clustering: misclass", misclass1, "\n")
    
    suppressWarnings(topics_fit <- topics(t(pooled_data), K=2, tol=10000))
    tclusters <- apply(topics_fit$omega, 1, function(x) return(which.max(x)))
    levels(tclusters) <- c("1","2")
    tab2 <- xtabs(~tclusters + factor1)
    if(dim(tab2)[1]==2){
     misclass2 <- (min(tab2[1,1]+tab2[2,2], tab2[1,2]+ tab2[2,1]))/100;
    }else{
      misclass2 <- 0.5
    }
    cat("Admixture clustering: misclass", misclass2, "\n")
    
    admix_prop[m,n] <- misclass2;
    hierarchy_prop[m,n] <- misclass1;
    admix_prop[n,m] <- admix_prop[m,n]
    hierarchy_prop[n,m] <- hierarchy_prop[m,n]
  }
  cat("We are at tissue", m, "\n")
}

save(admix_prop, file="admix_prop_gtex.rda")
save(hierarchy_prop, file="hierarchy_prop_gtex.rda")

