


#################   Guo 48 genes analysis    ##############################


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

gene_names <- get(load(file="../external_data/Deng_Data/TF_gene_names.rda"))


library(readxl)
guo_genes <- read_excel("../external_data/Deng_Data/guo48genes.xls")

guo_gene_names <- guo_genes$`Gene Symbol`
guo_gene_names <- setdiff(guo_gene_names, c("Actb"))

matched_indices <- match(guo_gene_names, deng.gene_names)
matched_indices <- matched_indices[!is.na(matched_indices)]

deng_counts_guo <- deng.counts[matched_indices,]

###  topic model

topic_clus <- maptpx::topics(t(deng_counts_guo), K=3, tol=0.1)
save(topic_clus, file="topic_fit_deng_counts_guo_k_4.rda")

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


deng_counts_guo_blast <- deng.counts[matched_indices, grep("blast", deng.meta_data$cell_type)]


topic_clus_blast <- CountClust::FitGoM(t(deng_counts_guo_blast), K=2, tol=0.1)
save(topic_clus_blast, file="../rdas/deng_guo_blast_47_genes_k_2.rda")

topic_clus_blast <- CountClust::FitGoM(t(deng_counts_guo_blast), K=3, tol=0.1)
save(topic_clus_blast, file="../rdas/deng_guo_blast_47_genes_k_3.rda")

topic_clus_blast <- maptpx::topics(t(deng_counts_guo_blast), K=4, tol=0.1)
save(topic_clus_blast, file="../rdas/deng_guo_blast_47_genes_k_4.rda")

CountClust::compGoM(t(deng_counts_guo_blast), topic_clus_blast)

omega <- topic_clus_blast$omega

annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(deng.meta_data$cell_type[grep("blast", deng.meta_data$cell_type)],
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


save(deng_counts_guo, file="../rdas/deng_counts_guo_48_genes.rda")


topic_clus <- get(load("topic_fit_deng_counts_guo_k_6.rda"))


library(CountClust)

topic_clus <- maptpx::topics(t(deng_counts_guo_blast), K=5,tol=0.1)
save(topic_clus, file="../rdas/topics_deng_guo_blast-k-5.rda")

out <- ExtractHighCorFeatures(topic_clus$omega,
                              deng_counts_guo_blast,
                              num_genes=20)

apply(out$feature_indices, c(1,2), function(x) return(guo_gene_names[x]))

omega <- topic_clus$omega

annotation <- data.frame(
    sample_id = paste0("X", c(1:NROW(omega))),
    tissue_label = factor(deng.meta_data$cell_type[grep("blast", deng.meta_data$cell_type)],
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





num_genes <- 20;
K <- 6
indices <- matrix(0,K,num_genes)
cor_values <- matrix(0,K,num_genes)

for(k in 1:K){
    cor_genes <- array(0, dim(topic_clus$theta)[1])
    for(l in 1:dim(topic_clus$theta)[1]){
        cor_genes[l] <- abs(cor(deng_counts_guo[l,], topic_clus$omega[,k]))
    }
    indices[k,] <- order(cor_genes, decreasing=TRUE)[1:num_genes]
    cor_values[k,] <- cor_genes[indices[k,]]
}

