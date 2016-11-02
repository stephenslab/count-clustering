

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
