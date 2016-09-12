
########   PCA on GTEx Brain  #######################

library(GTExV6Brain)
library(ggplot2)
library(CountClust)
counts <- exprs(GTExV6Brain)
meta_data <- pData(GTExV6Brain)
gene_names <- rownames(counts)


voom_counts <- limma::voom(counts)$E

pca_brain <- prcomp(t(voom_counts))
pca_brain_x <- pca_brain$x;

out <- pca_brain_x[,1:2]
save(out, file="../rdas/pca_gtex_brain.rda")

tsne_out <- tsne::tsne(pca_brain_x[,1:100],2)
save(tsne_out, file="../rdas/tsne_gtex_brain.rda")
