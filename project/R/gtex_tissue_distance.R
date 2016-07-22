

### GTEx heatmp distance matrix

cis_expr_data <- data.frame(data.table::fread("../external_data/GTEX_V6/cis_gene_expression.txt"))
matdata <- cis_expr_data[,-(1:2)];
voom_expr_data <- t(limma::voom(matdata)$E)

pr_voom_expr <- prcomp(voom_expr_data)$x;
dist_gtex <- as.matrix(dist(voom_expr_data, method="euclidean"))