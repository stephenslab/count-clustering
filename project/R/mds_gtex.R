library(limma)
library(data.table)
cis_expr_data <- data.frame(fread("cis_gene_expression.txt"))
matdata <- cis_expr_data[,-(1:2)];
voom_matdata <- voom(matdata)$E;

distcor <- 1 - cor(voom_matdata)

fit <- cmdscale(distcor, eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

out <- data.frame(x=x, y=y)

save(out, file="mds_gtex_whole.rda")



