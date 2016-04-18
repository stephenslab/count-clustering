

## Loglikelihood comparisons of the GTEx V6 model fits for K=15,16 and 17

topic_fit_15_1 <- get(load("../external_data/GTEX_V6/gtexv6fit.k.15.rda"));
counts <- data.frame(data.table::fread("../external_data/GTEX_V6/cis_gene_expression.txt"));
matdata <- counts[,-(1:2)];
topic_fit_15_1_list <- list(topic_fit_15_1)
compgom_15 <- compGoM(matdata, topic_fit_15_1_list)

loglik_15 <- as.numeric(compgom_15[2,1]);


topic_fit_16_1 <- get(load("../external_data/GTEX_V6/gtexv6fit.k.16.rda"));
topic_fit_16_1_list <- list(topic_fit_16_1)
compgom_16 <- compGoM(matdata, topic_fit_16_1_list)
loglik_16 <- as.numeric(compgom_16[2,1]);


topic_fit_16_2 <- get(load("../external_data/GTEX_V6/gtexv6fit.k.16.part2.rda"));
topic_fit_16_2_list <- list(topic_fit_16_2)
compgom_16_2 <- compGoM(matdata, topic_fit_16_2_list)
loglik_16_2 <- as.numeric(compgom_16_2[2,1]);

topic_fit_17_1 <- get(load("../external_data/GTEX_V6/gtexv6fit.k.17.rda"));
topic_fit_17_1_list <- list(topic_fit_17_1)
compgom_17 <- compGoM(matdata, topic_fit_17_1_list)
loglik_17 <- as.numeric(compgom_17[2,1]);
