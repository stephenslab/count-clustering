


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

library(readxl)
guo_genes <- read_excel("../external_data/Deng_Data/guo48genes.xls")

guo_gene_names <- guo_genes$`Gene Symbol`

matched_indices <- match(guo_gene_names, deng.gene_names)
matched_indices <- matched_indices[!is.na(matched_indices)]

deng_counts_guo <- deng.counts[matched_indices,]

###  topic model

topic_clus <- maptpx::topics(t(deng_counts_guo), K=6, tol=0.1)
save(topic_clus, file="topic_fit_deng_counts_guo_k_6.rda")

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
                figure_title = "Deng et al Structure Plot (on 48 genes due to Guo et al), k=2",
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))


deng_counts_guo_blast <- deng.counts[matched_indices, grep("blast", deng.meta_data$cell_type)]



topic_clus_blast <- maptpx::topics(t(deng_counts_guo_blast), K=2, tol=0.1)

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
                figure_title = "Deng et al blastocyst Structure Plot(on 48 genes due to Guo et al)",
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

