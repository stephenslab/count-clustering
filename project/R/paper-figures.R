###--- GTEx ---###

sample_labels <- c("Adipose Subcutaneous: 111", "Adipose Visceral: 22", "Adrenal Gland: 47", 
                    "Artery Aorta: 74", "Artery Coronary: 40", "Artery Tibial: 124", 
                    "Blood: 168", "Brain Amygdala: 26", "Brain Caudate Basal Ganglia: ", 
                    "Brain Cerebellar Hemisphere", "Brain Cerebellum: 29",
                    "Brain Cortex: 24", "Brain Frontal Cortex: 27",
                    "Brain Hypthallamus: 28", "Brain Nucleus Accumbens: 30",
                    "Brain Putamen: 23", "Brain Substantia Nigr: 26", "Breast: 57",
                    "Cell Line EBV Lymphocytes: 42", "Colon Transverse: 53",
                    "Esophagus Mucosa: 95", "Esophagus MUscularis: 91",
                    "Heart Atrial Appendage: 27", "Heart Left Ventrical: 87",
                    "Liver: 32", "Lung: 124", "Muscle Skeletal: 143",
                    "Nerve Tibial: 102", "Ovary: 30", "Pancreas: 58",
                    "Pituitary: 22", "Prostate: 37", "Skin Not Sun Exposed: 31",
                    "Skin Sun Exposed: 114", "Spleen: 28", "Stomach: 76", 
                    "Testis: 54", "Thyroid: 112", "Uterus: 32", "Vagina: 31")

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

# making unique cell sample labels
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
        
        # find the dominant cluster in each sample
        each_sample_order <- apply(temp_df, 1, which.max)
        
        # find the dominant cluster across samples
        sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])
        
        # reorder the matrix
        temp_df_ord <- temp_df[order(temp_df[ , sample_order]), ]
        
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
        scale_fill_brewer(palette = "Accent") +
        theme(legend.position = "none",
              axis.text = element_text(size = 4),
              title = element_text(size = 6)) +
        ggtitle("STRUCTURE plot by developmental phase") + 
        scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                            labels = seq(0, 1, 1/(ticks_number -1 ) ) ) + 
        coord_flip() 

# width = 1: increase bar width and in turn remove space
# between bars
b <- a + geom_bar(stat = "identity", 
                  position = "stack", 
                  width = 1)
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

          # find the dominant cluster in each sample
          each_sample_order <- apply(temp_df, 1, which.max)
          
          # find the dominant cluster across samples
          sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])
          
          # reorder the matrix
          temp_df_ord <- temp_df[order(temp_df[ , sample_order]), ]
          
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
    scale_fill_brewer(palette = "Accent") +
    theme(legend.position = "none",
          axis.text = element_text(size = 4),
          title = element_text(size = 6)) +
    ggtitle("STRUCTURE plot by developmental phase") + 
    scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                        labels = seq(0, 1, 1/(ticks_number -1 ) ) ) + 
    coord_flip() 

# width = 1: increase bar width to 1 so that there's no 
# space between bars
b <- a + geom_bar(stat = "identity", 
                  position = "stack", 
                  width = 1)
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


## TBD
## Experiment with circle stacked bar chart
# polarBarChart(
#         df = df_mlt,
#         binSize=1,
#         spaceBar=0.05,
#         spaceItem=0.2,
#         spaceFamily=1.2,
#         innerRadius=0.3,
#         outerRadius=1,
#         nguides=3,
#         guides=pretty(range(c(0, df$value)), n=nguides, min.n=2),
#         alphaStart=-0.3,
#         circleProportion=0.8,
#         direction="inwards",
#         familyLabels=TRUE,
#         itemSize=3,
#         legLabels=NULL,
#         legTitle="Source")




###--- Jaitin et al. 7 clusters ---###
# see here for analysis steps
# http://stephenslab.github.io/count-clustering/project/src/jaitin_structure_genes.html

# load the previously analyzed results
topic_fit <- readRDS("rdas/MouseJaitinSpleen-topicFit.rds")

# extract the omega matrix: membership weights of each cell
names(topic_fit$clust_7)
omega <- topic_fit$clust_7$omega

# load phenotype information
library(singleCellRNASeqMouseJaitinSpleen)
data(MouseJaitinSpleen)
meta_data <- colData(jaitin2014)
gene_names <- rownames(counts)

# follow previous anaysis steps to extract
# amplification batch information
# 
# exclude ERCC genes 
ENSG_genes_index <- grep("ERCC", gene_names, invert = TRUE)

# expression matrix without ENSG genes
counts_ensg <- counts[ENSG_genes_index, ]
filter_genes <- c("M34473","abParts","M13680","Tmsb4x",
                  "S100a4","B2m","Atpase6","Rpl23","Rps18",
                  "Rpl13","Rps19","H2-Ab1","Rplp1","Rpl4",
                  "Rps26","EF437368") 
fcounts <- counts_ensg[ -match(filter_genes, rownames(counts_ensg)), ]
sample_counts <- colSums(fcounts)

filter_sample_index <- which(meta_data$number_of_cells == 1 & 
                                 meta_data$group_name == "CD11c+" & 
                                 sample_counts > 600)

# make filterd phenotype data
meta_data_filtered <- meta_data[filter_sample_index, ]
stopifnot(dim(meta_data_filtered)[1] == dim(omega)[1])

# load packages
library(ggplot2)
library(plotly)
library(reshape2)
library(cowplot)


## Order by amplification batch

# making sample labels
amp_batch <- as.numeric(meta_data_filtered[ , "amplification_batch"])
cell_labels <- amp_batch
cell_labels_unique <- paste(amp_batch, c(1:length(amp_batch)), sep = "_")
cell_labels_fact <- factor(cell_labels,
                           levels = sort(unique(cell_labels)))

# make the re-ordered dataframe
number_plot_subgroups <- length(unique(cell_labels))
df_ord <- do.call(rbind,
                  lapply(1:number_plot_subgroups, function(ii) {
                      temp_label <- levels(cell_labels_fact)[ii]
                      temp_df <- omega[which(cell_labels == temp_label), ]
                      
                      # find the dominant cluster in each sample
                      each_sample_order <- apply(temp_df, 1, which.max)
                      
                      # find the dominant cluster across samples
                      sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])
                      
                      # reorder the matrix
                      temp_df_ord <- temp_df[order(temp_df[ , sample_order]), ]
                      
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
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none",
          axis.text = element_text(size = 4),
          title = element_text(size = 6)) +
    ggtitle("STRUCTURE plot by developmental phase") + 
    scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                        labels = seq(0, 1, 1/(ticks_number -1 ) ) ) + 
    coord_flip() 

# width = 1: increase bar width and in turn remove space
# between bars
b <- a + geom_bar(stat = "identity", 
                  position = "stack", 
                  width = 1)
b <- b + panel_border(remove = TRUE)
b <- ggdraw(switch_axis_position((b), axis = "y"))
cowplot::plot_grid(b)
save_plot("plots/jaitin-figures/amplification.png", 
          b, base_height = 4, base_width = 2)


# make axis labels
# take a barchart of the frequency of cell labels
# then copy and paste
png(file = "plots/jaitin-figures/amplfication-labels.png")
par(mfrow=c(1,2))
barplot(as.matrix(cell_labels_count), col = "white")
dev.off()




## Order by sequencing batch

# making sample labels
seq_batch <- as.numeric(meta_data_filtered[ , "sequencing_batch"])
cell_labels <- seq_batch
cell_labels_unique <- paste(seq_batch, c(1:length(amp_batch)), sep = "_")
cell_labels_fact <- factor(cell_labels,
                           levels = sort(unique(cell_labels)))

# make the re-ordered dataframe
number_plot_subgroups <- length(unique(cell_labels))
df_ord <- do.call(rbind,
                        lapply(1:number_plot_subgroups, function(ii) {
                            temp_label <- levels(cell_labels_fact)[ii]
                            temp_df <- omega[which(cell_labels == temp_label), ]
                            
                            # find the dominant cluster in each sample
                            each_sample_order <- apply(temp_df, 1, which.max)
                            
                            # find the dominant cluster across samples
                            sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])
                            
                            # reorder the matrix
                            temp_df_ord <- temp_df[order(temp_df[ , sample_order]), ]
                            
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
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "none",
          axis.text = element_text(size = 4),
          title = element_text(size = 6)) +
    ggtitle("STRUCTURE plot by developmental phase") + 
    scale_y_continuous( breaks = seq(0, value_ifl, length.out = ticks_number),
                        labels = seq(0, 1, 1/(ticks_number -1 ) ) ) + 
    coord_flip() 

# width = 1: increase bar width and in turn remove space
# between bars
b <- a + geom_bar(stat = "identity", 
                  position = "stack", 
                  width = 1)
b <- b + panel_border(remove = TRUE)
b <- ggdraw(switch_axis_position((b), axis = "y"))
cowplot::plot_grid(b)
save_plot("plots/jaitin-figures/sequencing.png", 
          b, base_height = 4, base_width = 2)


# make axis labels
# take a barchart of the frequency of cell labels
# then copy and paste
png(file = "plots/jaitin-figures/sequencing-labels.png")
par(mfrow=c(1,2))
barplot(as.matrix(cell_labels_count), col = "white")
dev.off()

