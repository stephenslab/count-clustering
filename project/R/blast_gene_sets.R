

###################   Genes Deng   ############################


geneset1 <- c("Creb3l2", "Tcf23", "Snai1","Pdgfra",
              "Gata4")

geneset2 <- c("Klf4", "Fgf4", "Klf2", "Esrrb",
              "Nanog", "Bmp4", "Sox2", "Fn1", "Pecam1")

geneset3 <- c("Tspan8", "Dppa1", "Lcp1", "Aqp3",
              "Id2", "Krt8", "Cebpa", "Eomes",
              "Gata3", "Mbnl3", "Tcfap2a", "Grhl1",
              "Atp12a", "Grhl2", "Cdx2")


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

match(geneset1, deng.gene_names)
match(geneset2, deng.gene_names)
match(geneset3, deng.gene_names)

deng_topics <- get(load("../rdas/deng_topic_fit.rda"))

deng_topics_clust6 <- deng_topics$clust_6

omega <- deng_topics_clust6$omega
theta <- deng_topics_clust6$theta

purple_prop <- deng_topics_clust6$omega[,2]/(deng_topics_clust6$omega[,2] + deng_topics_clust6$omega[,3])

voom_expr <- limma::voom(deng.counts)$E;

blast_indices <- grep("blast", deng.meta_data$cell_type)

dat <- data.frame("topic_prop"=purple_prop[blast_indices],
                  "log_cpm_expression"=t(voom_expr[,blast_indices]),
                  "labels"=droplevels(deng.meta_data$cell_type[blast_indices]))

graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Creb3l2, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Tcf23, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Snai1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Pdgfra, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Gata4, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 3,
                         nrow = 2))



graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Klf4, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Fgf4, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Klf2, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Esrrb, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Nanog, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Bmp4, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Sox2, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Fn1, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Pecam1, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 2))


graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Tspan8, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Dppa1, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Lcp1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Aqp3, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Id2, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Krt8, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Cebpa, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Eomes, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Gata3, data=dat, colour=labels)
graphList[[10]] <- qplot(topic_prop, log_cpm_expression.Mbnl3, data=dat, colour=labels)
graphList[[11]] <- qplot(topic_prop, log_cpm_expression.Tcfap2a, data=dat, colour=labels)
graphList[[12]] <- qplot(topic_prop, log_cpm_expression.Grhl1, data=dat, colour=labels)
graphList[[13]] <- qplot(topic_prop, log_cpm_expression.Atp12a, data=dat, colour=labels)
graphList[[14]] <- qplot(topic_prop, log_cpm_expression.Grhl2, data=dat, colour=labels)
graphList[[15]] <- qplot(topic_prop, log_cpm_expression.Cdx2, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 3))


##################  on 47  genes  #############################

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

gene_names <- get(load(file="../external_data/Deng_Data/TF_gene_names.rda"))


library(readxl)
guo_genes <- read_excel("../external_data/Deng_Data/guo48genes.xls")

guo_gene_names <- guo_genes$`Gene Symbol`
guo_gene_names <- setdiff(guo_gene_names, c("Actb"))

matched_indices <- match(guo_gene_names, deng.gene_names)
matched_indices <- matched_indices[!is.na(matched_indices)]

deng_counts_guo <- deng.counts[matched_indices,]

topic_clus <- maptpx::topics(t(deng_counts_guo), K=6, tol=0.1)
save(topic_clus, file="../rdas/topic_clus_deng_TF_47_genes.rda")

topic_clus <- get(load("../rdas/topic_clus_deng_TF_47_genes.rda"))

omega <- topic_clus$omega;
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(deng.meta_data$cell_type,
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
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


purple_prop <- omega[,2]

voom_expr <- limma::voom(deng.counts)$E;

blast_indices <- grep("blast", deng.meta_data$cell_type)
type <- droplevels(deng.meta_data$cell_type[blast_indices])

dat <- data.frame("topic_prop"=purple_prop[blast_indices],
                  "log_cpm_expression"=t(voom_expr[,blast_indices]),
                  "labels"=factor(type, levels=c("earlyblast", "midblast", "lateblast")))

graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Creb3l2, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Tcf23, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Snai1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Pdgfra, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Gata4, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 3,
                         nrow = 2))



graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Klf4, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Fgf4, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Klf2, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Esrrb, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Nanog, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Bmp4, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Sox2, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Fn1, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Pecam1, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 2))


graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Tspan8, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Dppa1, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Lcp1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Aqp3, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Id2, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Krt8, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Cebpa, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Eomes, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Gata3, data=dat, colour=labels)
graphList[[10]] <- qplot(topic_prop, log_cpm_expression.Mbnl3, data=dat, colour=labels)
graphList[[11]] <- qplot(topic_prop, log_cpm_expression.Tcfap2a, data=dat, colour=labels)
graphList[[12]] <- qplot(topic_prop, log_cpm_expression.Grhl1, data=dat, colour=labels)
graphList[[13]] <- qplot(topic_prop, log_cpm_expression.Atp12a, data=dat, colour=labels)
graphList[[14]] <- qplot(topic_prop, log_cpm_expression.Grhl2, data=dat, colour=labels)
graphList[[15]] <- qplot(topic_prop, log_cpm_expression.Cdx2, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 3))

#################   orange cluster   ##########################

orange_prop <- omega[,3]

voom_expr <- limma::voom(deng.counts)$E;

blast_indices <- grep("blast", deng.meta_data$cell_type)
type <- droplevels(deng.meta_data$cell_type[blast_indices])

dat <- data.frame("topic_prop"=orange_prop[blast_indices],
                  "log_cpm_expression"=t(voom_expr[,blast_indices]),
                  "labels"=factor(type, levels=c("earlyblast", "midblast", "lateblast")))

graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Creb3l2, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Tcf23, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Snai1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Pdgfra, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Gata4, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 3,
                         nrow = 2))



graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Klf4, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Fgf4, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Klf2, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Esrrb, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Nanog, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Bmp4, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Sox2, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Fn1, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Pecam1, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 2))


graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Tspan8, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Dppa1, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Lcp1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Aqp3, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Id2, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Krt8, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Cebpa, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Eomes, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Gata3, data=dat, colour=labels)
graphList[[10]] <- qplot(topic_prop, log_cpm_expression.Mbnl3, data=dat, colour=labels)
graphList[[11]] <- qplot(topic_prop, log_cpm_expression.Tcfap2a, data=dat, colour=labels)
graphList[[12]] <- qplot(topic_prop, log_cpm_expression.Grhl1, data=dat, colour=labels)
graphList[[13]] <- qplot(topic_prop, log_cpm_expression.Atp12a, data=dat, colour=labels)
graphList[[14]] <- qplot(topic_prop, log_cpm_expression.Grhl2, data=dat, colour=labels)
graphList[[15]] <- qplot(topic_prop, log_cpm_expression.Cdx2, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 3))


#################   yellow cluster   ##########################

yellow_prop <- omega[,4]

voom_expr <- limma::voom(deng.counts)$E;

blast_indices <- grep("blast", deng.meta_data$cell_type)
type <- droplevels(deng.meta_data$cell_type[blast_indices])

dat <- data.frame("topic_prop"=yellow_prop[blast_indices],
                  "log_cpm_expression"=t(voom_expr[,blast_indices]),
                  "labels"=factor(type, levels=c("earlyblast", "midblast", "lateblast")))

graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Creb3l2, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Tcf23, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Snai1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Pdgfra, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Gata4, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 3,
                         nrow = 2))



graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Klf4, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Fgf4, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Klf2, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Esrrb, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Nanog, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Bmp4, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Sox2, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Fn1, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Pecam1, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 2))


graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Tspan8, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Dppa1, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Lcp1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Aqp3, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Id2, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Krt8, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Cebpa, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Eomes, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Gata3, data=dat, colour=labels)
graphList[[10]] <- qplot(topic_prop, log_cpm_expression.Mbnl3, data=dat, colour=labels)
graphList[[11]] <- qplot(topic_prop, log_cpm_expression.Tcfap2a, data=dat, colour=labels)
graphList[[12]] <- qplot(topic_prop, log_cpm_expression.Grhl1, data=dat, colour=labels)
graphList[[13]] <- qplot(topic_prop, log_cpm_expression.Atp12a, data=dat, colour=labels)
graphList[[14]] <- qplot(topic_prop, log_cpm_expression.Grhl2, data=dat, colour=labels)
graphList[[15]] <- qplot(topic_prop, log_cpm_expression.Cdx2, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 3))

###############   Deng et al full (without Actb) ###################

deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

deng_gene_names_1 <- setdiff(deng.gene_names, c("Actb"))

matched_indices <- match(deng_gene_names_1, deng.gene_names)
matched_indices <- matched_indices[!is.na(matched_indices)]

deng_counts_1 <- deng.counts[matched_indices,]

topic_clus <- maptpx::topics(t(deng_counts_1), K=6, tol=100);
save(topic_clus, file="../rdas/topic_clus_no_Actb_deng.rda")

topic_clus <- get(load("../rdas/topic_clus_no_Actb_deng.rda"))

omega <- topic_clus$omega;
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(deng.meta_data$cell_type,
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
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

#################  2nd cluster (purple)  ##########################

purple_prop <- omega[,2]

voom_expr <- limma::voom(deng.counts)$E;

blast_indices <- grep("blast", deng.meta_data$cell_type)
type <- droplevels(deng.meta_data$cell_type[blast_indices])

dat <- data.frame("topic_prop"=purple_prop[blast_indices],
                  "log_cpm_expression"=t(voom_expr[,blast_indices]),
                  "labels"=factor(type, levels=c("earlyblast", "midblast", "lateblast")))

graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Creb3l2, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Tcf23, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Snai1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Pdgfra, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Gata4, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 3,
                         nrow = 2))



graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Klf4, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Fgf4, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Klf2, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Esrrb, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Nanog, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Bmp4, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Sox2, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Fn1, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Pecam1, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 2))


graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Tspan8, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Dppa1, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Lcp1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Aqp3, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Id2, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Krt8, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Cebpa, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Eomes, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Gata3, data=dat, colour=labels)
graphList[[10]] <- qplot(topic_prop, log_cpm_expression.Mbnl3, data=dat, colour=labels)
graphList[[11]] <- qplot(topic_prop, log_cpm_expression.Tcfap2a, data=dat, colour=labels)
graphList[[12]] <- qplot(topic_prop, log_cpm_expression.Grhl1, data=dat, colour=labels)
graphList[[13]] <- qplot(topic_prop, log_cpm_expression.Atp12a, data=dat, colour=labels)
graphList[[14]] <- qplot(topic_prop, log_cpm_expression.Grhl2, data=dat, colour=labels)
graphList[[15]] <- qplot(topic_prop, log_cpm_expression.Cdx2, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 3))

#################  3rd cluster (orange)  ##########################

orange_prop <- omega[,3]

voom_expr <- limma::voom(deng.counts)$E;

blast_indices <- grep("blast", deng.meta_data$cell_type)
type <- droplevels(deng.meta_data$cell_type[blast_indices])

dat <- data.frame("topic_prop"=orange_prop[blast_indices],
                  "log_cpm_expression"=t(voom_expr[,blast_indices]),
                  "labels"=factor(type, levels=c("earlyblast", "midblast", "lateblast")))

graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Creb3l2, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Tcf23, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Snai1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Pdgfra, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Gata4, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 3,
                         nrow = 2))



graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Klf4, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Fgf4, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Klf2, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Esrrb, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Nanog, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Bmp4, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Sox2, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Fn1, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Pecam1, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 2))


graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Tspan8, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Dppa1, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Lcp1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Aqp3, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Id2, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Krt8, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Cebpa, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Eomes, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Gata3, data=dat, colour=labels)
graphList[[10]] <- qplot(topic_prop, log_cpm_expression.Mbnl3, data=dat, colour=labels)
graphList[[11]] <- qplot(topic_prop, log_cpm_expression.Tcfap2a, data=dat, colour=labels)
graphList[[12]] <- qplot(topic_prop, log_cpm_expression.Grhl1, data=dat, colour=labels)
graphList[[13]] <- qplot(topic_prop, log_cpm_expression.Atp12a, data=dat, colour=labels)
graphList[[14]] <- qplot(topic_prop, log_cpm_expression.Grhl2, data=dat, colour=labels)
graphList[[15]] <- qplot(topic_prop, log_cpm_expression.Cdx2, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 3))

#################  4th cluster (yellow)  ##########################

yellow_prop <- omega[,4]

voom_expr <- limma::voom(deng.counts)$E;

blast_indices <- grep("blast", deng.meta_data$cell_type)
type <- droplevels(deng.meta_data$cell_type[blast_indices])

dat <- data.frame("topic_prop"=yellow_prop[blast_indices],
                  "log_cpm_expression"=t(voom_expr[,blast_indices]),
                  "labels"=factor(type, levels=c("earlyblast", "midblast", "lateblast")))

graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Creb3l2, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Tcf23, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Snai1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Pdgfra, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Gata4, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 3,
                         nrow = 2))



graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Klf4, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Fgf4, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Klf2, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Esrrb, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Nanog, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Bmp4, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Sox2, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Fn1, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Pecam1, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 2))


graphList <- vector(mode="list");
library(ggplot2)
graphList[[1]] <- qplot(topic_prop, log_cpm_expression.Tspan8, data=dat, colour=labels)
graphList[[2]] <- qplot(topic_prop, log_cpm_expression.Dppa1, data=dat, colour=labels)
graphList[[3]] <- qplot(topic_prop, log_cpm_expression.Lcp1, data=dat, colour=labels)
graphList[[4]] <- qplot(topic_prop, log_cpm_expression.Aqp3, data=dat, colour=labels)
graphList[[5]] <- qplot(topic_prop, log_cpm_expression.Id2, data=dat, colour=labels)
graphList[[6]] <- qplot(topic_prop, log_cpm_expression.Krt8, data=dat, colour=labels)
graphList[[7]] <- qplot(topic_prop, log_cpm_expression.Cebpa, data=dat, colour=labels)
graphList[[8]] <- qplot(topic_prop, log_cpm_expression.Eomes, data=dat, colour=labels)
graphList[[9]] <- qplot(topic_prop, log_cpm_expression.Gata3, data=dat, colour=labels)
graphList[[10]] <- qplot(topic_prop, log_cpm_expression.Mbnl3, data=dat, colour=labels)
graphList[[11]] <- qplot(topic_prop, log_cpm_expression.Tcfap2a, data=dat, colour=labels)
graphList[[12]] <- qplot(topic_prop, log_cpm_expression.Grhl1, data=dat, colour=labels)
graphList[[13]] <- qplot(topic_prop, log_cpm_expression.Atp12a, data=dat, colour=labels)
graphList[[14]] <- qplot(topic_prop, log_cpm_expression.Grhl2, data=dat, colour=labels)
graphList[[15]] <- qplot(topic_prop, log_cpm_expression.Cdx2, data=dat, colour=labels)

library(grid)
library(gridExtra)
a <- do.call("grid.arrange", 
             args = list(grobs=graphList,
                         ncol = 5,
                         nrow = 3))

out <- ExtractTopFeatures(topic_clus$theta, top_features = 100, method="poisson", options="min")

genes_extracted <- apply(out, c(1,2), function(x) return(deng_gene_names_1[x]))

which(!is.na(match(genes_extracted[4,],geneset1)))
