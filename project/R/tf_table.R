
#######################  TF genes  ################################

tf_table <- read.delim("../utilities/TF_genes.txt", sep="\t", header=FALSE)

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

thepage = readLines('http://www.bioguo.org/AnimalTFDB/BrowseGeneral.php?fam=AF-4&spe=Mus_musculus#Introduction')

out <- as.matrix(read.table("../utilities/TF_mouse.txt"))
tf_genes <- as.vector(out)

library(biomaRt)
listMarts(host="www.ensembl.org")
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="www.ensembl.org")

out <- getBM(
  attributes= c("ensembl_gene_id",
                "external_gene_name"),
  filters="ensembl_gene_id",
  values= tf_genes,
  mart= mart)

gene_names_tf <- out$external_gene_name;

matched_indices <- match(gene_names_tf, deng.gene_names)
matched_indices <- matched_indices[!is.na(matched_indices)]

deng_counts_tf <- deng.counts[matched_indices,]

###  topic model

topic_clus <- maptpx::topics(t(deng_counts_tf), K=6, tol=0.01)

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
                figure_title = "Deng et al Structure Plot(on TF genes)",
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

deng_counts_non_tf <- deng.counts[-matched_indices,]

topic_clus <- maptpx::topics(t(deng_counts_non_tf), K=6, tol=100)

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
                figure_title = "Deng et al Structure Plot(non TF genes)",
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))
