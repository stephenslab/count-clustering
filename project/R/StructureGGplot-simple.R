StructureGGplot_simple <- function (omega, annotation, palette = RColorBrewer::brewer.pal(8, 
                                                                "Accent"), figure_title = "", yaxis_label = "Tissue type", 
          order_sample = TRUE, sample_order_decreasing = TRUE, split_line = list(split_lwd = 1, 
                                                                                 split_col = "white"), plot_labels = TRUE, axis_tick = list(axis_ticks_length = 0.1, 
                                                                                                                                            axis_ticks_lwd_y = 0.1, axis_ticks_lwd_x = 0.1, axis_label_size = 3, 
                                                                                                                                            axis_label_face = "bold")) 
{
    if (dim(omega)[2] > length(palette)) {
        stop("Color choices is smaller than the number of clusters!")
    }
    if (length(unique(rownames(omega))) != NROW(omega)) {
        stop("omega rownames are not unique!")
    }
    if (!is.data.frame(annotation)) 
        stop("annotation must be a data.frame")
    if (!all.equal(colnames(annotation), c("sample_id", "tissue_label"))) {
        stop("annotation data.frame column names must be sample_id and tissue_label")
    }
    if (length(unique(annotation$sample_id)) != NROW(omega)) {
        stop("sample_id is not unique")
    }
    df_ord <- do.call(rbind, lapply(1:nlevels(annotation$tissue_label), 
                                    function(ii) {
                                        temp_label <- levels(annotation$tissue_label)[ii]
                                        temp_df <- omega[which(annotation$tissue_label == 
                                                                   temp_label), ]
                                        is_single_sample <- (length(temp_df) == nlevels(annotation$tissue_label) | 
                                                                 is.null(dim(temp_df)))
                                        if (is_single_sample) {
                                            each_sample_order <- which.max(temp_df)
                                        }
                                        else {
                                            each_sample_order <- apply(temp_df, 1, which.max)
                                        }
                                        sample_order <- as.numeric(attr(table(each_sample_order), 
                                                                        "name")[1])
                                        if (order_sample == TRUE & !is_single_sample) {
                                            temp_df_ord <- temp_df[order(temp_df[, sample_order], 
                                                                         decreasing = sample_order_decreasing), ]
                                        }
                                        else {
                                            temp_df_ord <- temp_df
                                        }
                                        temp_df_ord
                                    }))
    df_mlt <- reshape2::melt(t(df_ord))
    df_mlt <- plyr::rename(df_mlt, replace = c(Var1 = "topic", 
                                               Var2 = "document"))
    df_mlt$document <- factor(df_mlt$document)
    df_mlt$topic <- factor(df_mlt$topic)
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 12)) + ggplot2::theme_update(panel.grid.minor.x = ggplot2::element_blank(), 
                                                                                  panel.grid.minor.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), 
                                                                                  panel.grid.major.y = ggplot2::element_blank())
    value_ifl <- 10000
    ticks_number <- 6
    tissue_count <- table(droplevels(annotation$tissue_label))
    tissue_count_cumsum <- cumsum(table(droplevels(annotation$tissue_label)))
    tissue_names <- levels(droplevels(annotation$tissue_label))
    tissue_breaks <- sapply(1:length(tissue_count), function(i) {
        if (i == 1) {
            if (tissue_count[i] == 1) 
                bk <- 1
            if (tissue_count[i] > 1) 
                bk <- (tissue_count_cumsum[i] - 0)/2
            return(bk)
        }
        if (i > 1) {
            if (tissue_count[i] == 1) 
                bk_interval <- 1
            if (tissue_count[i] > 1) {
                bk_interval <- (tissue_count_cumsum[i] - tissue_count_cumsum[i - 
                                                                                 1])/2
            }
            bk <- tissue_count_cumsum[i - 1] + bk_interval
            return(bk)
        }
    })
    names(tissue_breaks) <- tissue_names
    a <- ggplot2::ggplot(df_mlt, ggplot2::aes(x = df_mlt$document, 
                                              y = df_mlt$value * 10000, fill = factor(df_mlt$topic))) + 
        ggplot2::xlab(yaxis_label) + ggplot2::ylab("") + ggplot2::scale_fill_manual(values = palette) + 
        ggplot2::theme(legend.position = "right", legend.key.size = ggplot2::unit(0.2, 
                                                                                  "cm"), legend.text = ggplot2::element_text(size = 5), 
                       axis.text = ggplot2::element_text(size = axis_tick$axis_label_size, 
                                                         face = axis_tick$axis_label_face), axis.ticks.y = ggplot2::element_line(size = axis_tick$axis_ticks_lwd_y), 
                       axis.ticks.x = ggplot2::element_line(size = axis_tick$axis_ticks_lwd_x), 
                       axis.ticks.length = ggplot2::unit(axis_tick$axis_ticks_length, 
                                                         "cm"), title = ggplot2::element_text(size = 6)) + 
        ggplot2::ggtitle(figure_title) + ggplot2::scale_y_continuous(breaks = seq(0, 
                                                                                  value_ifl, length.out = ticks_number), labels = seq(0, 
                                                                                                                                      1, 1/(ticks_number - 1))) + ggplot2::scale_x_discrete(breaks = as.character((levels(df_mlt$document)[round(tissue_breaks)])), 
                                                                                                                                                                                            labels = names(tissue_breaks)) + ggplot2::labs(fill = "Clusters") + 
        ggplot2::coord_flip()
    b <- a + ggplot2::geom_bar(stat = "identity", position = "stack", 
                               width = 1)
    b <- b + cowplot::panel_border(remove = TRUE)
    b <- b + ggplot2::geom_vline(xintercept = cumsum(table(droplevels(annotation$tissue_label)))[-length(table(droplevels(annotation$tissue_label)))] + 
                                     0.5, col = split_line$split_col, size = split_line$split_lwd)
    
    df_mlt
}