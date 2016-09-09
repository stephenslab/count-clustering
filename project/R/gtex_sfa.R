

###########  SFA on the GTEx 2013 data  #########################

library(flashr)

loadings_sfa <- read.table("../sfa_outputs/GTEX2013/voom_gtex/voom_gtex_sfa_lambda.out")
factors_sfa <- read.table("../sfa_outputs/GTEX2013/voom_gtex/voom_gtex_sfa_F.out")


library(CountClust)
samples_id <- read.table("../sfa_inputs/samples_id.txt");
tissue_labels <- factor(as.character(samples_id[,3]));

omega <- loadings_sfa

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  label = factor(as.character(tissue_labels),
                 levels=unique(tissue_labels)) )

rownames(omega) <- annotation$sample_id

cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3")[c(1,2,5,8,9)],
           RColorBrewer::brewer.pal(9, "Set1")[c(9,7)],
           RColorBrewer::brewer.pal(8, "Dark2")[c(3,4,8)])


FactorGGStack(loadings = omega,
              annotation = annotation,
              palette = cols1,
              yaxis_label = "Tissue Type",
              order_sample = TRUE,
              figure_title = "Factor Loadings Structure Plot (Sparse Loadings)",
              legend_labels = NULL,
              scale=TRUE,
              axis_tick = list(axis_ticks_length = .1,
                               axis_ticks_lwd_y = .1,
                               axis_ticks_lwd_x = .1,
                               axis_label_size = 5,
                               axis_label_face = "bold"))

loadings_sfa <- read.table("../sfa_outputs/GTEX2013_transpose/voom_gtex/gtex_voom_transpose_F.out")

omega <- t(loadings_sfa)

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  label = factor(as.character(tissue_labels),
                 levels=unique(tissue_labels)) )

rownames(omega) <- annotation$sample_id

cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3")[c(1,2,5,8,9)],
           RColorBrewer::brewer.pal(9, "Set1")[c(9,7)],
           RColorBrewer::brewer.pal(8, "Dark2")[c(3,4,8)])


FactorGGStack(loadings = omega,
              annotation = annotation,
              palette = cols1,
              yaxis_label = "Tissue Type",
              order_sample = TRUE,
              figure_title = "Factor Loadings Structure Plot (Sparse Factors)",
              legend_labels = NULL,
              scale=TRUE,
              axis_tick = list(axis_ticks_length = .1,
                               axis_ticks_lwd_y = .1,
                               axis_ticks_lwd_x = .1,
                               axis_label_size = 5,
                               axis_label_face = "bold"))
