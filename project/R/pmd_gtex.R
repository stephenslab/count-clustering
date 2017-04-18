library(limma)
library(PMA)
library(data.table)
cis_expr_data <- data.frame(fread("cis_gene_expression.txt"))
matdata <- cis_expr_data[,-(1:2)];
voom_matdata <- voom(matdata)$E;
cat("We start the PMD")
pmd_out <- PMD(t(voom_matdata), K=15);
save(pmd_out, file="PMD_on_GTEX_15.rda")

