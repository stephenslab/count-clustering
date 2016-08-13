

#####  Deng and GTEx Euclidean distance matrices #######

library(singleCellRNASeqMouseDeng2014)
counts <- exprs(Deng2014MouseESC)
meta_data <- pData(Deng2014MouseESC)
gene_names <- rownames(counts)

###  Voom (log CPM ) transformation

voom_counts <- limma::voom(counts)$E;

euclidean_mat_deng <- as.matrix(dist(t(voom_counts), method="euclidean"))

euclidean_object_deng <- list("euclidean_mat"=euclidean_mat_deng,
                              "metadata"=meta_data)

save(euclidean_object_deng, file="../rdas/euclidean_object_deng.rda")

euclidean_object_deng <- get(load("../rdas/euclidean_object_deng.rda"))


##########  GTEx V6 log CPM PCs data  ###################

library(data.table)
data <- data.frame(fread("../external_data/GTEX_V6/cis_gene_expression.txt"))
matdata <- data[,-(1:2)];

voom_matdata <- limma::voom(matdata)$E;

pc_data <- prcomp(t(voom_matdata));
