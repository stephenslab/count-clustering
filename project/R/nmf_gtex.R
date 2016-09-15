

#######  Non negative matrix factorization GTEx tissues ###############

library(NMF)
library(data.table)

cis_expr_data <- data.frame(fread("../external_data/GTEX_V6/cis_gene_expression.txt"))
matdata <- cis_expr_data[,-(1:2)];
voom_matdata <- limma::voom(matdata)$E;

out <- nmf(matdata, rank=20)
