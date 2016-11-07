

###############   BA9 and BA24  GTEx brain comparisons ######################


data <- get(load("../external_data/GTEX_V6/ba_data_gtex.rda"))
            
ba_metadata <- rownames(data);

library(CountClust)

#topic_clus <- CountClust::FitGoM(data, K=2:3, tol=10);
#save(topic_clus, file="../rdas/topic_clus_ba_data_2_3.rda")

topic_clus <- get(load(file="../rdas/topic_clus_ba_data_2_3.rda"))

omega <- topic_clus$clust_2$omega;
rownames(omega) <- 1:NROW(omega)
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(ba_metadata))

CountClust::StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(9, "Set1"),
                yaxis_label = "Amplification batch",
                order_sample = FALSE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


omega <- topic_clus$clust_3$omega;
rownames(omega) <- 1:NROW(omega)
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(ba_metadata))

CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(9, "Set1"),
                            yaxis_label = "Amplification batch",
                            order_sample = FALSE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


##############  list of extracted important genes for clusters  ######################

gene_names <- ExtractTopFeatures(topic_clus$clust_3$theta, top_features = 50, method = "poisson", options = "min")

gene_names_file <- as.vector(read.table("../external_data/GTEX_V6/gene_names_GTEX_V6.txt"))
gene_names_extract <- as.character(gene_names_file[,1])
gene_names_extract <- sapply(gene_names_extract, function(l) substring(l,1,15))
top_genes_3 <- apply(gene_names, c(1,2), function(l) return(gene_names_extract[l]))

save(top_genes_3, file="../rdas/ba_gtex_top_genes_clus_3.rda")


gene_names <- ExtractTopFeatures(topic_clus$clust_2$theta, top_features = 50, method = "poisson", options = "min")

gene_names_file <- as.vector(read.table("../external_data/GTEX_V6/gene_names_GTEX_V6.txt"))
gene_names_extract <- as.character(gene_names_file[,1])
gene_names_extract <- sapply(gene_names_extract, function(l) substring(l,1,15))
top_genes_2 <- apply(gene_names, c(1,2), function(l) return(gene_names_extract[l]))

save(top_genes_2, file="../rdas/ba_gtex_top_genes_clus_2.rda")
