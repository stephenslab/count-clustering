

####  Multi-dimensional scaling (Deng et al)  #############


library(singleCellRNASeqMouseDeng2014)
library(CountClust)
library(ggplot2)
counts <- exprs(Deng2014MouseESC)
meta_data <- pData(Deng2014MouseESC)
gene_names <- rownames(counts)

voom_counts <- limma::voom(counts)$E;

distcor <- 1 - cor(voom_counts)

fit <- cmdscale(distcor, eig = TRUE, k = 2)
x <- fit$points[,1]
y <- fit$points[,2]

out <- data.frame(x=x, y=y)
save(out, file="../rdas/mds_deng.rda")

plot(out$x, out$y, col=meta_data$cell_type)


######  Multi-dimensional scaling (GTEx Brain) ################

library(GTExV6Brain)
library(ggplot2)
library(CountClust)
counts <- exprs(GTExV6Brain)
meta_data <- pData(GTExV6Brain)
gene_names <- rownames(counts)

voom_counts <- limma::voom(counts)$E;

distcor <- 1 - cor(voom_counts)

fit <- cmdscale(distcor, eig = TRUE, k = 2)
x <- fit$points[,1]
y <- fit$points[,2]

out <- data.frame(x=x, y=y)
save(out, file="../rdas/mds_gtex_brain.rda")

plot(out$x, out$y, col=factor(meta_data$tissue_subtype))

########   Multi-dimensional scaling (GTEx whole tissues) #################


