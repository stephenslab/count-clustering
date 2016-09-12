

###############  NMF model on the GTEx whole tissues data  ###############################

ll <- get(load("../rdas/nmf_gtex_K_20.rda"))

library(ggplot2)
omega <-  t(apply(t(coef(ll)), 1, function(x) return(x/sum(x))))

samples_id <- read.table("../external_data/GTEX_V6/samples_id.txt");
tissue_labels <- factor(as.character(samples_id[,3]));

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(as.character(tissue_labels),
                 levels=unique(tissue_labels)) )

rownames(omega) <- annotation$sample_id

cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3")[c(1,2,5,8,9)],
           RColorBrewer::brewer.pal(9, "Set1")[c(9,7)],
           RColorBrewer::brewer.pal(8, "Dark2")[c(3,4,8)])

CountClust::StructureGGplot(omega = omega,
                            annotation= annotation,
                            palette = cols1,
                            yaxis_label = "",
                            order_sample = TRUE,
                            split_line = list(split_lwd = .1,
                                              split_col = "white"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 5,
                                             axis_label_face="bold"))

