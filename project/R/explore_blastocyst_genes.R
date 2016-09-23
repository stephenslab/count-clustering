


##########   exploring the purple and orange cluster genes ####################



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

deng_topics <- get(load("../rdas/deng_topic_fit.rda"))

deng_topics_clust6 <- deng_topics$clust_6

omega <- deng_topics_clust6$omega
theta <- deng_topics_clust6$theta

blast_indices <- grep("blast", deng.meta_data$cell_type)

omega[blast_indices,]

library(CountClust)
theta_mat <- deng_topics_clust6$theta;
top_features <- ExtractTopFeatures(theta_mat, top_features=100,
                                   method="poisson", options="min");
gene_list <- do.call(rbind, lapply(1:dim(top_features)[1],
                                   function(x) deng.gene_names[top_features[x,]]))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(rownames(omega),
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
                                        "lateblast") ) ) )

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

purple_gene_indices <- match(gene_list[2,], deng.gene_names)

deng.purple_gene_counts <- deng.counts[purple_gene_indices,];

topic_clus <- maptpx::topics(t(deng.purple_gene_counts), K=2, tol=0.1)

omega <- topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(deng.meta_data$cell_type,
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
                                        "lateblast") ) ) )

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"),
                figure_title = "Purple cluster genes (100) merged Structure")



orange_gene_indices <- match(gene_list[3,], deng.gene_names)

deng.orange_gene_counts <- deng.counts[orange_gene_indices,];

topic_clus <- maptpx::topics(t(deng.orange_gene_counts), K=2, tol=0.1)

omega <- topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(deng.meta_data$cell_type,
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
                                        "lateblast") ) ) )

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"),
                figure_title = "Orange cluster genes (100) merged Structure")


deng.purple_orange_gene_counts <- deng.counts[c(purple_gene_indices,orange_gene_indices),];

topic_clus <- maptpx::topics(t(deng.purple_orange_gene_counts), K=2, tol=0.1)

omega <- topic_clus$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(deng.meta_data$cell_type,
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
                                        "lateblast") ) ) )

rownames(omega) <- annotation$sample_id;

StructureGGplot(omega = omega,
                annotation = annotation,
                palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"),
                figure_title = "Orange + Purple cluster genes (200) merged Structure")

top_features <- ExtractTopFeatures(topic_clus$theta, top_features=100,
                                   method="poisson", options="min");


#####################  scatter plot  ###################################

theta_mat <- deng_topics_clust6$theta;
top_features <- ExtractTopFeatures(theta_mat, top_features=10,
                                   method="poisson", options="min");
gene_list <- do.call(rbind, lapply(1:dim(top_features)[1],
                                   function(x) deng.gene_names[top_features[x,]]))

purple_prop <- deng_topics_clust6$omega[,2]/(deng_topics_clust6$omega[,2] + deng_topics_clust6$omega[,3])
purple_genes_counts <- deng.counts[top_features[2,],];
purple_genes_voom_counts <- limma::voom(purple_genes_counts)$E;

dat <- data.frame("topic_prop"=purple_prop[blast_indices],
                  "log_cpm_expression"=t(purple_genes_voom_counts[,blast_indices]),
                  "labels"=deng.meta_data$cell_type[blast_indices])

graphList <- vector(mode="list");
library(ggplot2)

graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Upp1, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Tdgf1, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Aqp8, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Fabp5, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Tat, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Pdgfra, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Pyy, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Prdx1, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Col4a1, data=dat, colour=labels)
graphList[[10]] <- qplot(topic_prop, log_cpm_expression.Spp1, data=dat, colour=labels)


library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 2))

orange_prop <- deng_topics_clust6$omega[,3]/(deng_topics_clust6$omega[,3] + deng_topics_clust6$omega[,2])
orange_genes_counts <- deng.counts[top_features[3,],];
orange_genes_voom_counts <- limma::voom(orange_genes_counts)$E;

dat <- data.frame("topic_prop"=orange_prop[blast_indices],
                  "log_cpm_expression"=t(orange_genes_voom_counts[,blast_indices]),
                  "labels"=deng.meta_data$cell_type[blast_indices])

graphList <- vector(mode="list");
library(ggplot2)

graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Actb, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Krt18, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Fabp3, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Id2, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Tspan8, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Gm2a, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Lgals1, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Adh1, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Lrp2, data=dat, colour=labels)
graphList[[10]] <- qplot(topic_prop, log_cpm_expression.BC051665, data=dat, colour=labels)


library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 2))




plot(purple_genes_voom_counts[1,blast_indices], 
     purple_prop[blast_indices], col=deng.meta_data$cell_type,
     )
