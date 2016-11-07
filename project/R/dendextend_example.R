

########  Apply dendextend to the Deng et al and Iris data ################


#install.packages("dendextend")
#install.packages('dendextendRcpp')


library(dendextend)
library(dendextendRcpp)

iris <- datasets::iris
iris2 <- iris[,-5]
species_labels <- iris[,5]
library(colorspace) # get nice colors
species_col <- rev(rainbow_hcl(3))[as.numeric(species_labels)]


d_iris <- dist(iris2) # method="man" # is a bit better
hc_iris <- hclust(d_iris, method = "complete")
iris_species <- rev(levels(iris[,5]))

library(dendextend)
dend <- as.dendrogram(hc_iris)
species_col_adjusted <- species_col[hc_iris$order]
# order it the closest we can to the order of the observations:
#dend <- rotate(dend, 1:150)

# Color the branches based on the clusters:
#dend <- color_branches(dend, k=3) #, groupLabels=iris_species)

# Manually match the labels, as much as possible, to the real classification of the flowers:
#labels_colors(dend) <-
#  rainbow_hcl(3)[sort_levels_values(
#    as.numeric(iris[,5])[order.dendrogram(dend)]
#  )]

#labels(dend) <- paste(as.character(iris[,5])[order.dendrogram(dend)],
#                      "(",labels(dend),")",
#                      sep = "")
# We hang the dendrogram a bit:
labels(dend) <- rep("",length(species_col))

dend <- dend %>% set("leaves_pch", 15) %>% set("leaves_cex", 1) %>% set("leaves_col", species_col_adjusted)

#dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)
# And plot:
pdf("../plots/dendextend_iris.pdf")
par(mar = c(3,3,3,3))
plot(dend,
     main = "Clustered Iris data set
     (the labels give the true flower species)",
     horiz =  TRUE,  nodePar = list(cex = .007))
legend("topleft", legend = iris_species, fill = rainbow_hcl(3))
dev.off()

library(circlize)
circlize_dendrogram(dend)

pdf("../plots/circle_dendextend_iris.pdf")
par(mar = rep(0,4))
circlize_dendrogram(dend)
dev.off()


library(circlize)
rownames(iris) <- NULL
dend <- iris[,-5] %>% dist %>% hclust %>% as.dendrogram %>%
    set("labels_colors") %>% set("labels_cex", 0.5) %>%
    set("leaves_pch", 15) %>% set("leaves_cex", 1) %>% set("leaves_col", species_col_adjusted)
par(mar = rep(0,4))
circlize_dendrogram(dend, labels=FALSE)
legend("topleft", legend = iris_species, fill = rainbow_hcl(3))






#############  Deng et al   ########################

library(singleCellRNASeqMouseDeng2014)
library(CountClust)
library(ggplot2)
counts <- exprs(Deng2014MouseESC)
meta_data <- pData(Deng2014MouseESC)
gene_names <- rownames(counts)

voom_counts <- limma::voom(counts)$E;

celltype_labels <- meta_data$cell_type
celltype_labels <-   factor(celltype_labels, levels = c("zy",
                             "early2cell",
                             "mid2cell",
                             "late2cell",
                             "4cell",
                             "8cell",
                             "16cell",
                             "earlyblast",
                             "midblast",
                             "lateblast"))

numtypes <- length(unique(celltype_labels))
library(colorspace) # get nice colors
cols <- c( rev(c("darkblue", "blue", "cornflowerblue", "cadetblue2")),
           rev(c("darkgreen", "darkolivegreen4", "darkolivegreen3")),
           rev(c("coral4", "coral3", "coral")) )

species_col <- cols[as.numeric(celltype_labels)]

d_deng <- dist(t(voom_counts)) # method="man" # is a bit better
hc_deng <- hclust(d_deng, method = "complete")
species_col_adjusted <- species_col[hc_deng$order]

labels(hc_deng)
hc_deng$order
deng_celltypes <- levels(celltype_labels)

library(dendextend)
dend <- as.dendrogram(hc_deng)
# order it the closest we can to the order of the observations:
#dend <- rotate(dend, 1:259)

dend <- dend %>% set("leaves_pch", c(15)) %>% set("leaves_cex", c(2)) %>% set("leaves_col", species_col_adjusted)

#labels(dend) <- rep("", length(species_col))
# Color the branches based on the clusters:
#dend <- color_branches(dend, k=6) #, groupLabels=iris_species)

# Manually match the labels, as much as possible, to the real classification of the flowers:
#labels_colors(dend) <-
#  rainbow_hcl(numtypes)[sort_levels_values(
#    as.numeric(celltype_labels)[order.dendrogram(dend)]
#  )]

#labels(dend) <- paste(as.character(celltype_labels)[order.dendrogram(dend)],
#                     "(",labels(dend),")",
#                      sep = "")


dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)
# And plot:
pdf("../plots/dendextend_deng.pdf")
par(mar = c(3,3,3,3))
plot(dend,
     main = "Clustered Deng et al data set
     (Labels give true developmental stages)",
     horiz =  TRUE,  nodePar = list(cex = .007))
legend("topleft", legend = deng_celltypes, fill = rev(rainbow_hcl(numtypes)))
dev.off()

pdf("../plots/dendextend_deng_circle.pdf")
par(mar = c(3,3,3,3))
circlize_dendrogram(dend, labels=FALSE)
#legend("topleft", legend = deng_celltypes, 
#       fill = cols)
dev.off()


##############   GTEx  Brain 2013  ####################

library(GTExV6Brain)
library(ggplot2)
library(CountClust)
counts <- exprs(GTExV6Brain)
meta_data <- pData(GTExV6Brain)
gene_names <- rownames(counts)
gene_names <- sapply(1:length(gene_names), function(x) substring(gene_names[x],1,15))

voom_counts <- limma::voom(counts)$E;

celltype_labels <- factor(meta_data$tissue_subtype)
numtypes <- length(unique(celltype_labels))
library(colorspace) # get nice colors

library(RColorBrewer)
# make color scheme
col_vector <- c(RColorBrewer::brewer.pal(8, "Accent"),
                RColorBrewer::brewer.pal(8, "Dark2")[1:4],
                RColorBrewer::brewer.pal(12, "Paired")[1])

species_col <- col_vector[as.numeric(celltype_labels)]

#d_gtex_brain <- dist(t(voom_counts)) # method="man" # is a bit better
#hc_gtex_brain <- hclust(d_gtex_brain, method = "complete")
#save(hc_gtex_brain, file="../rdas/hc_gtex_brain.rda")

hc_gtex_brain <- get(load("../rdas/hc_gtex_brain.rda"))
species_col_adjusted <- species_col[hc_gtex_brain$order]

gtex_brain_celltypes <- rev(levels(celltype_labels))

library(dendextend)
dend <- as.dendrogram(hc_gtex_brain)
# order it the closest we can to the order of the observations:
#dend <- rotate(dend, 1:1259)

# Color the branches based on the clusters:
#dend <- color_branches(dend, k=6) #, groupLabels=iris_species)

# Manually match the labels, as much as possible, to the real classification of the flowers:
#labels_colors(dend) <-
#  rainbow_hcl(numtypes)[sort_levels_values(
#    as.numeric(celltype_labels)[order.dendrogram(dend)]
#  )]

#labels(dend) <- paste(as.character(celltype_labels)[order.dendrogram(dend)],
#                      "(",labels(dend),")",
#                      sep = "")

dend <- dend %>% set("leaves_pch", c(15)) %>% set("leaves_cex", c(2)) %>% set("leaves_col", species_col_adjusted)

#labels(dend) <- rep("",length(species_col))


#dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)
# And plot:
pdf("../plots/dendextend_gtex_brain.pdf")
par(mar = c(3,3,3,3))
plot(dend,
     main = "Clustered GTEx Brain data set
     (Labels give true tissue states)",
     horiz =  TRUE,  nodePar = list(cex = .007))
legend("topleft", legend = gtex_brain_celltypes, fill = rev(col_vector))
dev.off()

pdf("../plots/dendextend_gtex_brain_circle.pdf")
par(mar = c(3,3,3,3))
circlize_dendrogram(dend, labels=FALSE)
#legend("topleft", legend = gtex_brain_celltypes,
#       fill = rev(col_vector), cex=0.6)
dev.off()



###############   dendextend  GTEx whole tissues  #####################

cis_data <- data.frame(data.table::fread("../external_data/GTEX_V6/cis_gene_expression.txt"))
samples_id <- read.table("../external_data/GTEX_V6/samples_id.txt")
tissue_labels <- factor(samples_id[,3])
numtypes <- length(unique(tissue_labels))

#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#color1 <- sample(color, 53)

library(RColorBrewer)
n <- 53
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, 
                            qual_col_pals$maxcolors, rownames(qual_col_pals)))


species_col <- rev(col_vector)[as.numeric(tissue_labels)]


hc_gtex <- get(load("../rdas/hclust_gtex.rda"))
order <- hc_gtex[[3]]
species_col_adjusted <- species_col[order];
gtex_celltypes <- rev(levels(tissue_labels))

library(dendextend)
dend <- as.dendrogram(hc_gtex)
# order it the closest we can to the order of the observations:
dend <- rotate(dend, 1:8555)

# Color the branches based on the clusters:
#dend <- color_branches(dend, k=6) #, groupLabels=iris_species)

# Manually match the labels, as much as possible, to the real classification of the flowers:
#labels_colors(dend) <-
#  rainbow_hcl(numtypes)[sort_levels_values(
#    as.numeric(celltype_labels)[order.dendrogram(dend)]
#  )]

#labels(dend) <- paste(as.character(celltype_labels)[order.dendrogram(dend)],
#                      "(",labels(dend),")",
#                      sep = "")

dend <- dend %>% set("leaves_pch", c(15)) %>% set("leaves_cex", c(2)) %>%
    set("leaves_col", species_col_adjusted)

#labels(dend) <- rep("",length(species_col))


#dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)
# And plot:
pdf("../plots/dendextend_gtex.pdf")
par(mar = c(3,3,3,3))
plot(dend,
     main = "Clustered GTEx all tissues data
     (Labels give true tissue states)",
     horiz =  TRUE,  nodePar = list(cex = .007))
legend("topleft", legend = gtex_celltypes, 
       fill = color1[1:numtypes], ncol=2, cex=0.5)
dev.off()

pdf("../plots/dendextend_gtex_circle.pdf")
par(mar = c(3,3,3,3))
circlize_dendrogram(dend, labels=FALSE)
#legend("topleft", legend = rev(gtex_celltypes), 
#       fill = col_vector[1:numtypes], ncol=2, cex=0.5)
dev.off()