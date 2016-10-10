

##############  Deng et al final sub-analysis #####################


library(devtools)

read.data1 = function() {
  x = tempfile()
  download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
  z = get(load((x)))
  return(z)
}

Deng2014MouseESC <- read.data1()


deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

blast_indices <- grep("blast", deng.meta_data$cell_type)

blast.deng.counts <- deng.counts[,blast_indices]

#topic_clus <- maptpx::topics(t(blast.deng.counts), K=3, tol=100);
#save(topic_clus, file="../rdas/blast_deng_topic_fit_k_3_all_genes.rda")

topic_clus <- get(load("../rdas/blast_deng_topic_fit_k_3_all_genes.rda"))
blast <- rep("blastocyst", length(blast_indices))

omega <- topic_clus$omega
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(deng.meta_data$cell_type[grep("blast", deng.meta_data$cell_type)],
                        levels = rev( c("earlyblast","midblast",
                                        "lateblast") ) ) )
rownames(omega) <- annotation$sample_id;


library(CountClust)

StructureGGplot(omega = omega,
                annotation = annotation,
                figure_title = "Deng et al Structure Plot (on 46 genes due to Guo et al), k=4",
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


library(readxl)
guo_genes <- read_excel("../external_data/Deng_Data/guo48genes.xls")

guo_gene_names <- guo_genes$`Gene Symbol`
guo_gene_names <- setdiff(guo_gene_names, c("Actb", "Gapdh"))

matched_indices <- match(guo_gene_names, deng.gene_names)
matched_indices <- matched_indices[!is.na(matched_indices)]

deng_counts_guo <- deng.counts[matched_indices,]

blast.deng.counts.guo <- deng_counts_guo[, blast_indices];

topic_clus <- maptpx::topics(t(blast.deng.counts.guo), K=3, tol=0.1);
save(topic_clus, file="../rdas/blast_deng_topic_fit_k_3_48_genes.rda")

topic_clus <- get(load("../rdas/blast_deng_topic_fit_k_3_48_genes.rda"))
blast <- rep("blastocyst", length(blast_indices))

omega <- topic_clus$omega
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(deng.meta_data$cell_type[grep("blast", deng.meta_data$cell_type)],
                        levels = rev( c("earlyblast","midblast",
                                        "lateblast") ) ) )
rownames(omega) <- annotation$sample_id;


library(CountClust)

StructureGGplot(omega = omega,
                annotation = annotation,
                figure_title = "Deng et al Structure Plot (on 46 genes due to Guo et al), k=4",
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

indices <- ExtractTopFeatures(topic_clus$theta, top_features = 20, method="poisson",
                   options="min")

apply(indices, c(1,2), function(x) return(guo_gene_names[x]))


