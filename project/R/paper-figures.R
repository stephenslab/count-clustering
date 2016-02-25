###--- GTEx ---###
omega <- read.table("../project/rdas/omega_cis_genes_0_1_2.txt",
                    header = TRUE, sep = " ",
                    stringsAsFactors = FALSE)
dim(omega)
head(omega)
colnames(omega) <- c(1:NCOL(omega))


library(ggplot2)
library(plotly)
library(reshape2)
library(cowplot)


# make cell sample labels
# want a version consistent with majority of the literature
sample_labels <- read.table("../project/rdas/samples_id.txt", 
                            header = TRUE, sep = " ",
                            stringsAsFactors = FALSE)
cell_labels <- vector("numeric", NROW(sample_labels))

cell_labels[sample_labels$SMTSD == "Adipose - Subcutaneous"] <- "Adipose Subcutaneous"
cell_labels[sample_labels$SMTSD == "Adipose - Visceral (Omentum)"] <- "Adipose Visceral"
cell_labels[sample_labels$SMTSD == "Adrenal Gland"] <- "Adrenal Gland"
cell_labels[sample_labels$SMTSD == "Artery - Aorta"] <- "Artery Aorta"
cell_labels[sample_labels$SMTSD == "Artery - Coronary"] <- "Artery Coronary"
cell_labels[sample_labels$SMTSD == "Artery - Tibial"] <- "Artery Tibial"
cell_labels[sample_labels$SMTSD == "Bladder"] <- "Bladder"
cell_labels[sample_labels$SMTSD == "Brain - Anterior cingulate cortex (BA24)"] <- "Brain Anterior Cingulate Cortex"
cell_labels[sample_labels$SMTSD == "Brain - Amygdala"] <- "Brain Amygdala"
cell_labels[sample_labels$SMTSD == "Brain - Caudate (basal ganglia)"] <- "Brain Caudate (basal ganglia)"
cell_labels[sample_labels$SMTSD == "Brain - Cerebellar Hemisphere"] <- "Brain Cerebellar Hemisphere"
cell_labels[sample_labels$SMTSD == "Brain - Cerebellum"] <- "Brain Cerebellum"
cell_labels[sample_labels$SMTSD == "Brain - Cortex"] <- "Brain Cortex"
cell_labels[sample_labels$SMTSD == "Brain - Frontal Cortex (BA9)"] <- "Brain Frontal Cortex"
cell_labels[sample_labels$SMTSD == "Brain - Hippocampus"] <- "Brain Hippocampus"
cell_labels[sample_labels$SMTSD == "Brain - Hypothalamus"] <- "Brain Hypothalamus"
cell_labels[sample_labels$SMTSD == "Brain - Nucleus accumbens (basal ganglia)"] <- "Brain Nucleus Accumbens (basal ganglia)"
cell_labels[sample_labels$SMTSD == "Brain - Putamen (basal ganglia)"] <- "Brain Putamen (basal ganglia)"
cell_labels[sample_labels$SMTSD == "Brain - Spinal cord (cervical c-1)"] <- "Brain Spinal cord"
cell_labels[sample_labels$SMTSD == "Brain - Substantia nigra"] <- "Brain - Substantia Nigr"
cell_labels[sample_labels$SMTSD == "Breast - Mammary Tissue"] <- "Breast"
cell_labels[sample_labels$SMTSD == "Cells - EBV-transformed lymphocytes"] <- "Cells Line EBV Lymphocytes"
cell_labels[sample_labels$SMTSD == "Cells - Transformed fibroblasts"] <- "Cells Transformed Fibroblasts"
cell_labels[sample_labels$SMTSD == "Cervix - Ectocervix"] <- "Cervix Ectocervix"
cell_labels[sample_labels$SMTSD == "Cervix - Endocervix"] <- "Cervix Endocervix"
cell_labels[sample_labels$SMTSD == "Colon - Sigmoid"] <- "Colon Sigmoid"
cell_labels[sample_labels$SMTSD == "Colon - Transverse"] <- "Colon Transverse"
cell_labels[sample_labels$SMTSD == "Esophagus - Gastroesophageal Junction"] <- "Esophagus Gastroesophageal Junction"
cell_labels[sample_labels$SMTSD == "Esophagus - Mucosa"] <- "Esophagus Mucosa"
cell_labels[sample_labels$SMTSD == "Esophagus - Muscularis"] <- "Esophagus Muscularis"
cell_labels[sample_labels$SMTSD == "Fallopian Tube"] <- "Fallopian Tube"
cell_labels[sample_labels$SMTSD == "Heart - Atrial Appendage"] <- "Heart Atrial Appendage"
cell_labels[sample_labels$SMTSD == "Heart - Left Ventricle"] <- "Heart Left Ventricle"
cell_labels[sample_labels$SMTSD == "Kidney - Cortex"] <- "Kidney Cortex"
cell_labels[sample_labels$SMTSD == "Liver"] <- "Liver"
cell_labels[sample_labels$SMTSD == "Lung"] <- "Lung"
cell_labels[sample_labels$SMTSD == "Minor Salivary Gland"] <- "Minor Salivary Gland"
cell_labels[sample_labels$SMTSD == "Muscle - Skeletal"] <- "Muscle Skeletal"
cell_labels[sample_labels$SMTSD == "Nerve - Tibial"] <- "Nerve Tibial"
cell_labels[sample_labels$SMTSD == "Ovary"] <- "Ovary"
cell_labels[sample_labels$SMTSD == "Pancreas"] <- "Pancreas"
cell_labels[sample_labels$SMTSD == "Pituitary"] <- "Pituitary"
cell_labels[sample_labels$SMTSD == "Prostate"] <- "Prostate"
cell_labels[sample_labels$SMTSD == "Skin - Not Sun Exposed (Suprapubic)"] <- "Skin Not Sun Exposed"
cell_labels[sample_labels$SMTSD == "Skin - Sun Exposed (Lower leg)"] <- "Skin Sun Exposed"
cell_labels[sample_labels$SMTSD == "Small Intestine - Terminal Ileum"] <- "Small Intestine Terminal Ileum"
cell_labels[sample_labels$SMTSD == "Spleen"] <- "Spleen"
cell_labels[sample_labels$SMTSD == "Stomach"] <- "Stomach"
cell_labels[sample_labels$SMTSD == "Testis"] <- "Testis"
cell_labels[sample_labels$SMTSD == "Thyroid"] <- "Thyroid"
cell_labels[sample_labels$SMTSD == "Uterus"] <- "Uterus"
cell_labels[sample_labels$SMTSD == "Vagina"] <- "Vagina"
cell_labels[sample_labels$SMTSD == "Whole Blood"] <- "Whole Blood"

# assign tissue labels
rownames(omega) <- paste(cell_labels, 1:length(cell_labels), sep = "_")
cell_labels_fact <- factor(cell_labels, 
                           levels = sort(unique(cell_labels)) )


# find tissue orders using hierarchical clustering
# this is needs to be done before re-ordering the data
docweights_per_family_mean <- apply(omega, 2, function(x) {
    tapply(x, cell_labels, mean) })
family_ordering <- heatmap(docweights_per_family_mean)$rowInd
levels(cell_labels_fact) <- levels(cell_labels_fact)[family_ordering]


# make the re-ordered dataframe
df_ord <- do.call(rbind,
                  lapply(1:nlevels(cell_labels_fact), function(ii) {
                      temp_label <- rev(levels(cell_labels_fact))[ii]
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

# transform for ggplot2 input
df_mlt <- reshape2::melt(t(df_ord))
df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                           "Var2" = "document"))

# start making polar histogram
source("R/polarHistogramStructure.R")

# rename columns to satisfy function input
df_plot <- df_mlt
colnames(df_plot) <- c("score", "item", "value")
df_plot$family <- sapply(df_mlt$document, function(x) {
    strsplit(as.character(x), split= "_")[[1]][1]
})
#df_mlt$value <- df2$value*10000
df_plot$family <- factor(df_plot$family)
df_plot$score <- factor(df_plot$score)
df_plot$item <- sapply(strsplit(as.character(df_plot$item), split = "_"), function(x) x[2])
df_plot$item <- factor(df_plot$item)
levels(df_plot$family) <- levels(cell_labels_fact)

# define colors of the clusers
cols1 <- c(RColorBrewer::brewer.pal(9, "Set1"),
          RColorBrewer::brewer.pal(12, "Set3"))

cols2 <- c("red", "blue", "cornflowerblue", "black", "cyan", "darkblue",
            "brown4", "burlywood", "darkgoldenrod1", "darkgray", "deepskyblue",
            "darkkhaki", "firebrick", "darkorchid", "hotpink", "green",
            "magenta", "yellow", "azure1", "azure4")

source("R/polarHistogramStructure.R")

# being preparing for computing on midway
save(polarHistogramStructure,
     cols1, cols2,
     df_plot, file = "rdas/rda-for-midway.rda")

# the following code was run on midway to make the plot
# run large interactive job
# sinteractive --partition=bigmem --constraint=256G --ntasks=1 --cpus-per-task=1 --mem-per-cpu=128000

png(file = "main-figure-cols1.png",
    height = 7, width = 7, units = "in", res = 300)
polarHistogramStructure(df2, palette = cols1,
                        outerRadius = 1.8,
                        innerRadius = 0.3,
                        familyLabelDistance = 2, # same metric as outerRadius
                        binSize = .5,
                        circleProportion = 0.97,
                        familyLabels = TRUE)
dev.off()

png(file = "main-figure-cols2.png",
    height = 7, width = 7, units = "in", res = 300)
polarHistogramStructure(df2, palette = cols2,
                        outerRadius = 1.8,
                        innerRadius = 0.3,
                        familyLabelDistance = 2, # same metric as outerRadius
                        binSize = .5,
                        circleProportion = 0.97,
                        familyLabels = TRUE)
dev.off()




###------ Deng et al. 6 clusters ------###


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
cols <- "Accent"
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

