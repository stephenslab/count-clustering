library(data.table)
library(dendextend)
library(dendextendRcpp)
library(limma)

cis_expr_data <- data.frame(fread("cis_gene_expression.txt"))
matdata <- cis_expr_data[,-(1:2)];
voom_matdata <- voom(matdata)$E;

d_gtex <- dist(t(voom_matdata));
hc_gtex <- hclust(d_gtex, method = "complete")
save(hc_gtex, file="hclust_gtex.rda")


