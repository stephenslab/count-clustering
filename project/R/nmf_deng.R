

#########  Non-negative matrix factorization on Deng et al #####################


library(singleCellRNASeqMouseDeng2014)
library(CountClust)
library(ggplot2)
counts <- exprs(Deng2014MouseESC)
meta_data <- pData(Deng2014MouseESC)
gene_names <- rownames(counts)

indices <- which(rowSums(counts)==0)
counts <- counts[-indices, ]


library(NMF)

out <- nmf(counts, rank=6)
save(out, file="../rdas/nmf_deng_k_6.rda")

out <- get(load("nmf_deng_k_6.rda"))

library(ggplot2)
omega <-  t(apply(t(coef(out)), 1, function(x) return(x/sum(x))))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(rownames(omega),
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
                                        "lateblast") ) ) )

rownames(omega) <- annotation$sample_id;

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(meta_data$cell_type,
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
                                        "lateblast") ) ) )
rownames(omega) <- annotation$sample_id; #(edited)



FactorGGStack(loadings = omega,
             annotation = annotation,
             palette = c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(4, "Spectral")),
             yaxis_label = "Development Phase",
             order_sample = TRUE,
             figure_title = "Factor Loadings Structure Plot",
             legend_labels = pve_percentage[-1],
             scale=TRUE,
             axis_tick = list(axis_ticks_length = .1,
                              axis_ticks_lwd_y = .1,
                              axis_ticks_lwd_x = .1,
                              axis_label_size = 7,
                              axis_label_face = "bold"))


StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))
