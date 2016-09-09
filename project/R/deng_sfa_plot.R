
############  Deng SFA plot  ##############################

library(devtools)
library(flashr)

read.data1 = function() {
  x = tempfile()
  download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
  z = get(load((x)))
  return(z)
}

Deng2014MouseESC <- read.data1()

deng.counts <- Biobase::exprs(Deng2014MouseESC)
meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

lambda_out <- read.table("../sfa_outputs/Deng2014PPS/voom_deng_counts_lambda.out")

omega <- lambda_out

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  label = factor(meta_data$cell_type,
                 levels = c("zy", "early2cell",
                            "mid2cell", "late2cell","4cell", "8cell", "16cell","earlyblast","midblast","lateblast") ) ) 

rownames(omega) <- annotation$sample_id

FactorGGStack(loadings = omega[,1:6],
             annotation = annotation,
             palette = c(RColorBrewer::brewer.pal(8, "Accent"),
                         RColorBrewer::brewer.pal(4, "Spectral")),
             yaxis_label = "Development Phase",
             order_sample = TRUE,
             figure_title = "Factor Loadings Structure Plot (Sparse Loadings)",
             legend_labels = NULL,
             scale=TRUE,
             axis_tick = list(axis_ticks_length = .1,
                              axis_ticks_lwd_y = .1,
                              axis_ticks_lwd_x = .1,
                              axis_label_size = 7,
                              axis_label_face = "bold"))


lambda_out <- read.table("../sfa_outputs/Deng2014PPS/voom_deng_counts_transpose_F.out")

omega <- t(lambda_out)

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  label = factor(meta_data$cell_type,
                 levels = c("zy", "early2cell",
                            "mid2cell", "late2cell","4cell", "8cell", "16cell","earlyblast","midblast","lateblast") ) ) 

rownames(omega) <- annotation$sample_id

FactorGGStack(loadings = omega[,1:6],
              annotation = annotation,
              palette = c(RColorBrewer::brewer.pal(8, "Accent"),
                          RColorBrewer::brewer.pal(4, "Spectral")),
              yaxis_label = "Development Phase",
              order_sample = TRUE,
              figure_title = "Factor Loadings Structure Plot (Sparse Factors)",
              legend_labels = NULL,
              scale=TRUE,
              axis_tick = list(axis_ticks_length = .1,
                               axis_ticks_lwd_y = .1,
                               axis_ticks_lwd_x = .1,
                               axis_label_size = 7,
                               axis_label_face = "bold"))


