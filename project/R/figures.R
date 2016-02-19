

###--- Deng et al. 6 clusters ---###


# load the previously analyzed results
load("../../../count-clustering/project/rdas/deng_topic_fit.rda")

# extract the omega matrix: membership weights of each cell
names(Topic_clus_list)
str(Topic_clus_list$clust_6)
omega <- Topic_clus_list$clust_6$omega

library(ggplot2)
library(plotly)
library(reshape2)
library(cowplot)

# order omega so that the figure looks nicer...
cell_labels <- rownames(omega)
cell_labels_unique <- paste(cell_labels, 1:NROW(omega), sep = "_")
rownames(omega) <- cell_labels_unique

cell_labels_fact <- factor(cell_labels,
                          levels = c("zy", "early2cell", "mid2cell", "late2cell",
                                     "4cell", "8cell", "16cell", "earlyblast",
                                     "midblast", "lateblast"))

## Display for all cell types

# make the re-ordered dataframe
df_ord <- do.call(rbind,
    lapply(1:nlevels(cell_labels_fact), function(ii) {
    temp_label <- levels(cell_labels_fact)[ii]
    temp_df <- omega[which(cell_labels == temp_label), ]
    first_order <- order(temp_df[1, ], decreasing = FALSE)
    sample_order <- do.call(cbind,
                            lapply(1:length(first_order),
                                   function(ii) {
                                       which_column <- first_order[ii]
                                       order(temp_df[ ,which_column], 
                                             decreasing = FALSE)
                                   }) )
    temp_df_ord <- temp_df[order(sample_order[ ,1], 
                                 sample_order[ ,2],
                                 sample_order[ ,3], 
                                 sample_order[ ,4], 
                                 sample_order[ ,5], 
                                 sample_order[ ,6]), ]
    temp_df_ord
}) )
dim(df_ord)

df_mlt <- reshape2::melt(t(df_ord))
df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                           "Var2" = "document"))
head(df_mlt)

# set blank background
theme_set(theme_bw(base_size = 12)) +
    theme_update( panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_blank() )

value_ifl <- 10000
ticks_number <- 6
# Use RColorBrewer color
# http://bxhorn.com/rcolorbrewer-palettes/
a <- ggplot(df_mlt, 
            aes(x = document, y = value*10000, fill = factor(topic)) ) + 
        xlab("Cell types") + ylab("") +
        scale_fill_brewer(palette = "Set3") +
        theme(legend.position = "none",
              axis.text = element_text(size = 4),
              title = element_text(size = 6)) +
        ggtitle("STRUCTURE plot by developmental phase") + 
        scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                            labels = seq(0, 1, 1/(ticks_number -1 ) ) ) + 
        coord_flip() 

b <- a + geom_bar(stat = "identity", position = "stack")
b <- b + panel_border(remove = TRUE)
b <- ggdraw(switch_axis_position((b), axis = "y"))
cowplot::plot_grid(b)
save_plot("../../../count-clustering/project/plots/deng-figures/deng-ggplot.png", 
          b, base_height = 4, base_width = 2)
#ggplotly(b)

# make axis labels
# take a barchart of the frequency of cell labels
# then copy and paste
cell_labels_count <- table(cell_labels_fact)
png(file = "../../../count-clustering/project/plots/deng-figures/labels.png")
par(mfrow=c(1,2))
barplot(as.matrix(cell_labels_count), col = "white")
dev.off()




## Display for some early stages

# make the re-ordered dataframe
df_ord <- do.call(rbind,
                  lapply(5:7, function(ii) {
                      temp_label <- levels(cell_labels_fact)[ii]
                      temp_df <- omega[which(cell_labels == temp_label), ]
                      first_order <- order(temp_df[1, ], decreasing = FALSE)
                      sample_order <- do.call(cbind,
                                              lapply(1:length(first_order),
                                                     function(ii) {
                                                         which_column <- first_order[ii]
                                                         order(temp_df[ ,which_column], 
                                                               decreasing = FALSE)
                                                     }) )
                      temp_df_ord <- temp_df[order(sample_order[ ,1], 
                                                   sample_order[ ,2],
                                                   sample_order[ ,3], 
                                                   sample_order[ ,4], 
                                                   sample_order[ ,5], 
                                                   sample_order[ ,6]), ]
                      temp_df_ord
                  }) )
dim(df_ord)

df_mlt <- reshape2::melt(t(df_ord))
df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                           "Var2" = "document"))
head(df_mlt)

# set blank background
theme_set(theme_bw(base_size = 12)) +
    theme_update( panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_blank() )

value_ifl <- 10000
ticks_number <- 6
# Use RColorBrewer color
# http://bxhorn.com/rcolorbrewer-palettes/
a <- ggplot(df_mlt, 
            aes(x = document, y = value*10000, fill = factor(topic)) ) + 
    xlab("Cell types") + ylab("") +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "none",
          axis.text = element_text(size = 4),
          title = element_text(size = 6)) +
    ggtitle("STRUCTURE plot by developmental phase") + 
    scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                        labels = seq(0, 1, 1/(ticks_number -1 ) ) ) + 
    coord_flip() 

b <- a + geom_bar(stat = "identity", position = "stack")
b <- b + panel_border(remove = TRUE)
b <- ggdraw(switch_axis_position((b), axis = "y"))
cowplot::plot_grid(b)
save_plot("../../../count-clustering/project/plots/deng-figures/deng-ggplot-early.png", 
          b, base_height = 4, base_width = 2)
#ggplotly(b)

# make axis labels
# take a barchart of the frequency of cell labels
# then copy and paste
cell_labels_count <- table(cell_labels[cell_labels %in% c("4cell", "8cell", "16cell")])
png(file = "../../../count-clustering/project/plots/deng-figures/labels-early.png")
par(mfrow=c(1,2))
barplot(as.matrix(cell_labels_count), col = "white")
dev.off()

