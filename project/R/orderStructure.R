#' Order samples by phenotype for Structure plot
#' 
orderStructure <- function(omega, annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            figure_title = "",
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            sample_order_decreasing = TRUE,
                            split_line = list(split_lwd = 1,
                                              split_col = "white"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 3) ) {
    
    
    # # rename columns to satisfy function input
    # df_mlt <- reshape2::melt(t(omega_ordered))
    # df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
    #                                            "Var2" = "document"))
    # df_mlt$document <- factor(df_mlt$document)
    # df_mlt$topic <- factor(df_mlt$topic)
    # head(df_mlt)
    # df_plot <- df_mlt
    # colnames(df_plot) <- c("score", "item", "value")
    # df_plot$family <- annotation$tissue_label[match(df_mlt$document, annotation$sample_id)]
    # df_plot$family <- factor(df_plot$family)
    # df_plot$score <- factor(df_plot$score)
    # df_plot$item <- factor(df_plot$item)
    # levels(df_plot$family) <- levels(annotation$tissue_label)
    