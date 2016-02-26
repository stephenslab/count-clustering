###--- GTEx ---###

# prepare data

omega <- read.table("../project/rdas/omega_cis_genes_0_1_2.txt",
                    header = TRUE, sep = " ",
                    stringsAsFactors = FALSE)
dim(omega)
colnames(omega) <- c(1:NCOL(omega))


# make cell sample labels
# want a version consistent with majority of the literature
sample_labels <- read.table("../project/rdas/samples_id.txt", 
                            header = TRUE, sep = " ",
                            stringsAsFactors = FALSE)
tissue_labels <- vector("numeric", NROW(sample_labels))
tissue_labels <- sample_labels[ ,3]


# clean labels
tissue_labels[grep("Nucleus", tissue_labels)] <- "Brain -N. accumbens"
tissue_labels[grep("Putamen", tissue_labels)] <- "Brain -Putamen"
tissue_labels[grep("Caudate", tissue_labels)] <- "Brain -Caudate"
tissue_labels[grep("Gastroe", tissue_labels)] <- "Esophagus -Gastroesophageal Jn."
tissue_labels[grep("cingulate", tissue_labels)] <- "Brain - Anterior cortex (BA24)."
tissue_labels[grep("EBV", tissue_labels)] <- "Cells -EBV-lymphocytes"
tissue_labels[grep("Suprapubic", tissue_labels)] <- "Skin - Unexposed (Suprapubic)"
tissue_labels[grep("Lower Leg", tissue_labels)] <- "Skin - Sun Exposed (Lower Leg)"


# find sample orders in hierarchical clustering
tissue_label <- as.character(tissue_labels)
docweights_per_tissue_mean <- apply(omega, 2, 
        function(x) { tapply(x, tissue_labels, mean) })
ordering <- heatmap(docweights_per_tissue_mean)$rowInd


# # order samples by hierarchical clustering order
samples_order <- unlist(
    lapply(1:53, function(x) which(tissue_labels == unique(tissue_labels)[ordering][x])))
#tissue_labels <- as.factor(tissue_labels)
tissue_label <- factor(tissue_label,
                       levels = rownames(docweights_per_tissue_mean)[ordering])  
omega_reordered <- omega[samples_order, ]
tissue_labels_reordered <- tissue_labels[samples_order]

# manually move the samples to make plot look better
indices1 <- which(tissue_labels_ordered =="Artery - Coronary")
indices2 <- 1:7667
indices3 <- 7668:8227
indices4 <- 8361:8555
indices <- c(indices2, indices1, indices3, indices4)
tissue_labels_reordered <- tissue_labels_reordered[indices]
omega_reordered <- omega_reordered[indices, ]


# assign tissue labels
rownames(omega_reordered) <- paste0("X", 1:length(tissue_labels_reordered))
annotation <- data.frame(
    sample_id = paste0("X", 1:length(tissue_labels_reordered)),
    tissue_label = factor(tissue_labels_reordered,
                          levels = unique(tissue_labels_reordered)) )


# flip omega & annotation matrix so that the brain tissues get plotted on top
# omega_reordered <- omega_reordered[c(1:NROW(omega_reordered)), ]
# annotation <- annotation[c(1:NROW(annotation)), ]




###<----- Structure ggplot

source("R/StructureGGplot.R")

# define colors of the clusers
cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3"))

cols2 <- c("red", "blue", "cornflowerblue", "black", "cyan", "darkblue",
           "brown4", "burlywood", "darkgoldenrod1", "darkgray", "deepskyblue",
           "darkkhaki", "firebrick", "darkorchid", "hotpink", "green",
           "magenta", "yellow", "azure1", "azure4")

# Joyce's color scheme
gtex_ggplot_cols1 <- StructureGGplot(omega = omega_reordered, 
                                     annotation = annotation,
                                     palette = cols1,
                                     figure_title = "",
                                     yaxis_label = "", 
                                     order_sample = FALSE)
save_plot("plots/gtex-figures/main-barplot-cols13.pdf", 
          gtex_ggplot_cols1, base_height = 4, base_width = 2)

# Kushal's color scheme
gtex_ggplot_cols2 <- StructureGGplot(omega = omega_ordered, 
                               annotation = annotation,
                               palette = cols2,
                               figure_title = "",
                               yaxis_label = "",
                               order_sample = FALSE)
save_plot("plots/gtex-figures/main-barplot-cols2.pdf", 
          gtex_ggplot_cols2, base_height = 4, base_width = 2)


# print out a barplot of tissue type counts
png(file = "plots/gtex-figures/tissue-type-labels.png")
par(mfrow=c(1,2))
barplot(as.matrix(table(annotation$tissue_label)), col = "white")
dev.off()




# ##<----- Polar histogram of all tissues

source("R/polarHistogramStructure.R")

# rename columns to satisfy function input
df_mlt <- reshape2::melt(t(omega_ordered))
df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                           "Var2" = "document"))
df_mlt$document <- factor(df_mlt$document)
df_mlt$topic <- factor(df_mlt$topic)
head(df_mlt)
df_plot <- df_mlt
colnames(df_plot) <- c("score", "item", "value")
df_plot$family <- annotation$tissue_label[match(df_mlt$document, annotation$sample_id)]
df_plot$family <- factor(df_plot$family)
df_plot$score <- factor(df_plot$score)
df_plot$item <- factor(df_plot$item)
levels(df_plot$family) <- levels(annotation$tissue_label)


# being preparing for computing on midway
save(polarHistogramStructure,
     cols1, cols2,
     df_plot, file = "rdas/rda-for-midway-all.rda")

# the following code was run on midway to make the plot
# run large interactive job
# sinteractive --partition=bigmem --constraint=256G --ntasks=1 --cpus-per-task=1 --mem-per-cpu=128000
pdf("main-figure-polar-histogram-all.pdf",
    height = 12, width = 12)
polarHistogramStructure(df_plot, palette = cols1,
                       outerRadius = 1.8,
                       innerRadius = 0.3,
                       familyLabelDistance = 2, # same metric as outerRadius
                       binSize = 1,
                       spaceFamily = 4,
                       circleProportion = 0.90,
                       familyLabels = TRUE)
dev.off()





##<----- Polar histogram of brain tissues

source("R/polarHistogramStructure.R")

omega <- read.table("../project/rdas/omega_cis_genes_brain.txt",
                    header = TRUE, sep = " ",
                    stringsAsFactors = FALSE)
dim(omega)
colnames(omega) <- c(1:NCOL(omega))
head(omega)


# make cell sample labels
# want a version consistent with majority of the literature
sample_labels <- read.table("../project/rdas/samples_id.txt", 
                            header = TRUE, sep = " ",
                            stringsAsFactors = FALSE)
brain_labels <- sample_labels[grep("Brain", sample_labels[,3]), 3]

# assign tissue labels
rownames(omega) <- paste0("X", 1:length(brain_labels))
annotation <- data.frame(
    sample_id = paste0("X", 1:length(brain_labels)),
    tissue_label = factor(brain_labels,
                          levels = c("Brain - Cerebellar Hemisphere",
                                     "Brain - Cerebellum",
                                     "Brain - Spinal cord (cervical c-1)",
                                     "Brain - Anterior cingulate cortex (BA24)",
                                     "Brain - Frontal Cortex (BA9)",
                                     "Brain - Cortex",
                                     "Brain - Hippocampus",
                                     "Brain - Substantia nigra",
                                     "Brain - Amygdala",
                                     "Brain - Putamen (basal ganglia)",
                                     "Brain - Caudate (basal ganglia)",
                                     "Brain - Hypothalamus",
                                     "Brain - Nucleus accumbens (basal ganglia)") ) )
                                     
# define colors of the clusers
cols <- c("blue", "darkgoldenrod1", "cyan", "red")

# StructureGGplot(omega,
#                 annotation= annotation, palette = cols, order_sample = FALSE)

# rename columns to satisfy function input
df_mlt <- reshape2::melt(t(omega))
df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                           "Var2" = "document"))

df_mlt$document <- factor(df_mlt$document)
df_mlt$topic <- factor(df_mlt$topic)
head(df_mlt)
df_plot <- df_mlt
colnames(df_plot) <- c("score", "item", "value")
df_plot$family <- annotation$tissue_label[match(df_mlt$document, annotation$sample_id)]
df_plot$family <- factor(df_plot$family)
df_plot$score <- factor(df_plot$score)
df_plot$item <- factor(df_plot$item)
levels(df_plot$family) <- levels(annotation$tissue_label)


# preparing for computing on midway
save(polarHistogramStructure,
     cols1, cols2, cols,
     df_plot, file = "rdas/rda-for-midway.rda")


# the following code was run on midway to make the plot
# run large interactive job
# 
# sinteractive --partition=bigmem --constraint=256G --ntasks=1 --cpus-per-task=1 --mem-per-cpu=128000
pdf("main-figure-polar-histogram.pdf",
    height = 8, width = 8)
    polarHistogramStructure(
        df = df_plot, 
        palette = cols,
        outerRadius = 1.8,
        innerRadius = 0.3,
        familyLabelDistance = 2, # same metric as outerRadius
        binSize = 1,
        spaceFamily = 4,
        circleProportion = 0.90,
        familyLabels = TRUE)
dev.off()


 





###------ Deng et al. 6 clusters ------###

source("R/StructureGGplot.R")

# load the previously analyzed results
load("rdas/deng_topic_fit.rda")

# extract the omega matrix: membership weights of each cell
names(Topic_clus_list)
str(Topic_clus_list$clust_6)
omega <- Topic_clus_list$clust_6$omega

# make annotation matrix
annotation <- data.frame(
  sample_id = c(1:NROW(omega)),
  tissue_label = factor(rownames(omega),
                        levels = rev( c("zy", "early2cell", "mid2cell", "late2cell",
                                       "4cell", "8cell", "16cell", "earlyblast",
                                       "midblast", "lateblast") ) ) )

# after extracting tissue type of each sample
# recode each sample to have unique rownames
rownames(omega) <- paste0("X",annotation$sample_id)

# make the polar histogram
deng_ggplot <- StructureGGplot(omega = omega, 
                               annotation = annotation,
                               palette = RColorBrewer::brewer.pal(8, "Accent"),
                               figure_title = "",
                               yaxis_label = "Cell type")
    
cowplot::save_plot("plots/deng-figures/deng-ggplot.pdf", 
                   deng_ggplot, base_height = 4, base_width = 2)






##<----- Display for some early stages

# make the re-ordered dataframe
df_ord <- do.call(rbind,
      lapply(5:7, function(ii) {
          temp_label <- levels(annotation$tissue_label)[ii]
          temp_df <- omega[which(annotation$tissue_label == temp_label), ]

          # find the dominant cluster in each sample
          each_sample_order <- apply(temp_df, 1, which.max)
          
          # find the dominant cluster across samples
          sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])
          
          # reorder the matrix
          temp_df_ord <- temp_df[order(temp_df[ , sample_order]), ]
          
          temp_df_ord
      }) )
dim(df_ord)

#df_mlt <- reshape2::melt(t(df_ord))
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






##<----- Polar histogram of Deng et al., 

source("R/polarHistogramStructure.R")

dim(omega)
colnames(omega) <- c(1:NCOL(omega))
head(omega)


# make cell sample labels
# want a version consistent with majority of the literature
sample_labels <- read.table("../project/rdas/samples_id.txt", 
                            header = TRUE, sep = " ",
                            stringsAsFactors = FALSE)
brain_labels <- sample_labels[grep("Brain", sample_labels[,3]), 3]

# assign tissue labels
rownames(omega) <- paste0("X", 1:length(brain_labels))
annotation <- data.frame(
    sample_id = paste0("X", 1:length(brain_labels)),
    tissue_label = factor(brain_labels,
                          levels = c("Brain - Cerebellar Hemisphere",
                                     "Brain - Cerebellum",
                                     "Brain - Spinal cord (cervical c-1)",
                                     "Brain - Anterior cingulate cortex (BA24)",
                                     "Brain - Frontal Cortex (BA9)",
                                     "Brain - Cortex",
                                     "Brain - Hippocampus",
                                     "Brain - Substantia nigra",
                                     "Brain - Amygdala",
                                     "Brain - Putamen (basal ganglia)",
                                     "Brain - Caudate (basal ganglia)",
                                     "Brain - Hypothalamus",
                                     "Brain - Nucleus accumbens (basal ganglia)") ) )

# define colors of the clusers
cols <- c("blue", "darkgoldenrod1", "cyan", "red")

# StructureGGplot(omega,
#                 annotation= annotation, palette = cols, order_sample = FALSE)

# rename columns to satisfy function input
df_mlt <- reshape2::melt(t(omega))
df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                           "Var2" = "document"))

df_mlt$document <- factor(df_mlt$document)
df_mlt$topic <- factor(df_mlt$topic)
head(df_mlt)
df_plot <- df_mlt
colnames(df_plot) <- c("score", "item", "value")
df_plot$family <- annotation$tissue_label[match(df_mlt$document, annotation$sample_id)]
df_plot$family <- factor(df_plot$family)
df_plot$score <- factor(df_plot$score)
df_plot$item <- factor(df_plot$item)
levels(df_plot$family) <- levels(annotation$tissue_label)


# preparing for computing on midway
save(polarHistogramStructure,
     cols1, cols2, cols,
     df_plot, file = "rdas/rda-for-midway.rda")


# the following code was run on midway to make the plot
# run large interactive job
# 
# sinteractive --partition=bigmem --constraint=256G --ntasks=1 --cpus-per-task=1 --mem-per-cpu=128000
pdf("main-figure-polar-histogram.pdf",
    height = 8, width = 8)
polarHistogramStructure(
    df = df_plot, 
    palette = cols,
    outerRadius = 1.8,
    innerRadius = 0.3,
    familyLabelDistance = 2, # same metric as outerRadius
    binSize = 1,
    spaceFamily = 4,
    circleProportion = 0.90,
    familyLabels = TRUE)
dev.off()




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
cols <- "Set1"
# Use RColorBrewer color
# http://bxhorn.com/rcolorbrewer-palettes/
a <- ggplot(df_mlt, 
            aes(x = document, y = value*10000, fill = factor(topic)) ) + 
    xlab("Cell types") + ylab("") +
    scale_fill_brewer(palette = cols) +
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

