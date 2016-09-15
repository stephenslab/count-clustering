

########  Apply dendextend to the Deng et al and Iris data ################


install.packages("dendextend")
install.packages('dendextendRcpp')


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
# order it the closest we can to the order of the observations:
dend <- rotate(dend, 1:150)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=3) #, groupLabels=iris_species)

# Manually match the labels, as much as possible, to the real classification of the flowers:
labels_colors(dend) <-
  rainbow_hcl(3)[sort_levels_values(
    as.numeric(iris[,5])[order.dendrogram(dend)]
  )]

labels(dend) <- paste(as.character(iris[,5])[order.dendrogram(dend)],
                      "(",labels(dend),")", 
                      sep = "")
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
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


#############  Deng et al   ########################

library(singleCellRNASeqMouseDeng2014)
library(CountClust)
library(ggplot2)
counts <- exprs(Deng2014MouseESC)
meta_data <- pData(Deng2014MouseESC)
gene_names <- rownames(counts)

voom_counts <- limma::voom(counts)$E;

celltype_labels <- meta_data$cell_type
numtypes <- length(unique(celltype_labels))
library(colorspace) # get nice colors
species_col <- rev(rainbow_hcl(numtypes))[as.numeric(celltype_labels)]

d_deng <- dist(t(voom_counts)) # method="man" # is a bit better
hc_deng <- hclust(d_deng, method = "complete")
deng_celltypes <- rev(levels(celltype_labels))

library(dendextend)
dend <- as.dendrogram(hc_deng)
# order it the closest we can to the order of the observations:
dend <- rotate(dend, 1:259)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=6) #, groupLabels=iris_species)

# Manually match the labels, as much as possible, to the real classification of the flowers:
labels_colors(dend) <-
  rainbow_hcl(numtypes)[sort_levels_values(
    as.numeric(celltype_labels)[order.dendrogram(dend)]
  )]

labels(dend) <- paste(as.character(celltype_labels)[order.dendrogram(dend)],
                      "(",labels(dend),")", 
                      sep = "")


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
legend("topleft", legend = deng_celltypes, fill = rainbow_hcl(numtypes))
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
species_col <- rev(rainbow_hcl(numtypes))[as.numeric(celltype_labels)]

d_gtex_brain <- dist(t(voom_counts)) # method="man" # is a bit better
hc_gtex_brain <- hclust(d_gtex_brain, method = "complete")
gtex_brain_celltypes <- rev(levels(celltype_labels))

library(dendextend)
dend <- as.dendrogram(hc_gtex_brain)
# order it the closest we can to the order of the observations:
dend <- rotate(dend, 1:1259)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=6) #, groupLabels=iris_species)

# Manually match the labels, as much as possible, to the real classification of the flowers:
labels_colors(dend) <-
  rainbow_hcl(numtypes)[sort_levels_values(
    as.numeric(celltype_labels)[order.dendrogram(dend)]
  )]

labels(dend) <- paste(as.character(celltype_labels)[order.dendrogram(dend)],
                      "(",labels(dend),")", 
                      sep = "")


dend <- hang.dendrogram(dend,hang_height=0.1)
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
legend("topleft", legend = gtex_brain_celltypes, fill = rainbow_hcl(numtypes))
dev.off()
