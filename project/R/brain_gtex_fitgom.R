

## Fit GoM on Brain GTEx data

library(data.table)
data <- data.frame(fread('../external_data/GTEX_V6/cis_gene_expression.txt'));
matdata <- data[,-(1:2)];
samples_id=read.table("../external_data/GTEX_V6/samples_id.txt")[,3];
brain_indices <- grep("Brain", samples_id);

brain_data <- matdata[,brain_indices];
colnames(brain_data) <- samples_id[brain_indices];


source("../../../maptpx/R/topics.R")
source("../../../maptpx/R/tpx.R")
source("../../../maptpx/R/count.R")

library(slam)
library(maptpx)
Topic_clus <- topics(t(brain_data), K=5, tol=100, method_admix=1);
