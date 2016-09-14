

##############  multiple hierarchical model runs  #######################

counts <- data.frame(data.table::fread("../external_data/GTEX_V6/cis_gene_expression.txt"))
matdata <- counts[,-(1:2)]
voom_out  <- limma::voom(matdata)
voom_weights <- voom_out$weights
voom_data <- t(voom_out$E)

samples_id <- read.table("../external_data/GTEX_V6/samples_id.txt")
tissue_labels <- samples_id[,3]
library(DESeq2)

library(edgeR)

cpm_data <- cpm(matdata, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)


cpm_data_norm <- t(apply(cpm_data, 1, function(x) return((x-mean(x))/sd(x))))

counts_adjusted <- t(apply(matdata, 1, function(x) return(x/sqrt(mean(x)+1))))

counts_adjusted_norm <- t(apply(counts_adjusted, 1, function(x) return((x-mean(x))/sd(x))))



hierarchy_prop_1 <- matrix(0, length(tissues_to_consider), length(tissues_to_consider))
hierarchy_prop_2 <- matrix(0, length(tissues_to_consider), length(tissues_to_consider))
hierarchy_prop_3 <- matrix(0, length(tissues_to_consider), length(tissues_to_consider))
hierarchy_prop_4<- matrix(0, length(tissues_to_consider), length(tissues_to_consider))

library(maptpx)
library(slam)
source("../../../maptpx/R/count.R")
source("../../../maptpx/R/topics.R")
source("../../../maptpx/R/tpx.R")

for(m in 2:length(tissues_to_consider)){
  for(n in 1:(m-1)){
    tissue1 <- tissues_to_consider[m];
    tissue2 <- tissues_to_consider[n];
    index1 <- which(tissue_labels==tissue1)
    index2 <- which(tissue_labels==tissue2)
    
    cpm_data1 <- cpm_data[,index1[1:50]]
    cpm_data2 <- cpm_data[,index2[1:50]]
    cpm_data_pooled <- cbind(cpm_data1, cpm_data2)
    
    factor1 <- c(rep(1,50), rep(2,50))
    dd <- dist(t(cpm_data_pooled))
    
    hclusters <- cutree(hclust(dd),2)
    tab <- xtabs(~hclusters + factor1)
    misclass1 <- (min(tab[1,1]+tab[2,2], tab[1,2]+ tab[2,1]))/100;
    cat("Hierarchical clustering (cpm): misclass", misclass1, "\n")
    
    cpm_data_norm_1 <- cpm_data_norm[,index1[1:50]]
    cpm_data_norm_2 <- cpm_data_norm[,index2[1:50]]
    cpm_data_norm_pooled <- cbind(cpm_data_norm_1, cpm_data_norm_2)
    
    factor1 <- c(rep(1,50), rep(2,50))
    dd <- dist(t(cpm_data_norm_pooled))
    
    hclusters <- cutree(hclust(dd),2)
    tab <- xtabs(~hclusters + factor1)
    misclass2 <- (min(tab[1,1]+tab[2,2], tab[1,2]+ tab[2,1]))/100;
    cat("Hierarchical clustering (cpm norm): misclass", misclass2, "\n")
    
    counts_adjusted1 <- counts_adjusted[,index1[1:50]]
    counts_adjusted2 <- counts_adjusted[,index2[1:50]]
    counts_adjusted_pooled <- cbind(counts_adjusted1, counts_adjusted2)
    
    factor1 <- c(rep(1,50), rep(2,50))
    dd <- dist(t(counts_adjusted_pooled))
    
    hclusters <- cutree(hclust(dd),2)
    tab <- xtabs(~hclusters + factor1)
    misclass3 <- (min(tab[1,1]+tab[2,2], tab[1,2]+ tab[2,1]))/100;
    cat("Hierarchical clustering (counts adjusted): misclass", misclass3, "\n")
    
    counts_adjusted_norm1 <- counts_adjusted_norm[,index1[1:50]]
    counts_adjusted_norm2 <- counts_adjusted_norm[,index2[1:50]]
    counts_adjusted_norm_pooled <- cbind(counts_adjusted_norm1, counts_adjusted_norm2)
    
    factor1 <- c(rep(1,50), rep(2,50))
    dd <- dist(t(counts_adjusted_norm_pooled))
    
    hclusters <- cutree(hclust(dd),2)
    tab <- xtabs(~hclusters + factor1)
    misclass4 <- (min(tab[1,1]+tab[2,2], tab[1,2]+ tab[2,1]))/100;
    cat("Hierarchical clustering (counts adjusted norm) class", misclass4, "\n")
    
    hierarchy_prop_1[m,n] <- misclass1;
    hierarchy_prop_2[m,n] <- misclass2;
    hierarchy_prop_3[m,n] <- misclass3;
    hierarchy_prop_4[m,n] <- misclass4;
    
  }
  cat("We are at tissue", m, "\n")
}

save(hierarchy_prop_1, file="../rdas/hierarchy_prop_1.rda")
save(hierarchy_prop_2, file="../rdas/hierarchy_prop_2.rda")
save(hierarchy_prop_3, file="../rdas/hierarchy_prop_3.rda")
save(hierarchy_prop_4, file="../rdas/hierarchy_prop_4.rda")

