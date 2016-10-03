
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

#gene_names <- get(load(file="../external_data/Deng_Data/TF_gene_names.rda"))


library(readxl)
guo_genes <- read_excel("../external_data/Deng_Data/guo48genes.xls")

guo_gene_names <- guo_genes$`Gene Symbol`

matched_indices <- match(guo_gene_names, deng.gene_names)
matched_indices <- matched_indices[!is.na(matched_indices)]

deng_counts_guo <- deng.counts[matched_indices,]

counts <- t(deng_counts_guo)
K <- 6
shape=NULL
initopics=NULL
tol=0.1
bf=FALSE
kill=2
ord=TRUE
verb=1
admix=TRUE
nbundles=1
use_squarem=FALSE
init.adapt=FALSE
light=1
method_admix=1
sample_init=TRUE
tmax=10000


library(maptpx)
library(slam)
topic_clus <- topics(t(deng_counts_guo), K=6, tol=0.1, start_init=TRUE)

omega <- topic_clus$omega

topic_clus$D$dispersion

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
                figure_title = "Deng et al blastocyst Structure Plot(on 48 genes due to Guo et al)",
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


guo_gene_names <- setdiff(guo_gene_names, c("Actb"))

matched_indices <- match(guo_gene_names, deng.gene_names)
matched_indices <- matched_indices[!is.na(matched_indices)]

deng_counts_guo <- deng.counts[matched_indices,]

library(maptpx)
library(slam)
topic_clus <- topics(t(deng_counts_guo), K=6, tol=0.1, start_init=TRUE)

omega <- topic_clus$omega

topic_clus$D$dispersion

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
                figure_title = "Deng et al blastocyst Structure Plot(on 47 genes due to Guo et al)",
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))




