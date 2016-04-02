
### Topic model fit ###

library(data.table)
## devtools::install_github('kkdey/maptpx')
library(maptpx)
gtexv6_data <- data.frame(fread("cis_gene_expression.txt"))


matdata <- t(gtexv6_data[,-c(1,2)]);

Topic_clus_1 <- topics(matdata, K=16, tol=100)
#Topic_clus_2 <- topics(matdata, K=18, tol=100)
#Topic_clus_3 <- topics(matdata, K=20, tol=100)

save(Topic_clus_1, file="gtexv6fit.k.16.rda")

