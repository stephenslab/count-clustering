


##################   Deng et al (odd chromosomes vs even chromosome genes)  ####################


library(singleCellRNASeqMouseDeng2014)
library(CountClust)
library(ggplot2)
counts <- exprs(Deng2014MouseESC)
meta_data <- pData(Deng2014MouseESC)
gene_names <- rownames(counts)

seq(1,24,2)
library(biomaRt)
mart<- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
chr_odd_genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), 
                    filters = 'chromosome_name', values = seq(1,24,2), mart = mart)

chr_even_genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), 
                       filters = 'chromosome_name', values = seq(2,24,2), mart = mart)

matched_indices_odd <- match(gene_names, chr_odd_genes$external_gene_name)

gene_names_odd <- gene_names[which(!is.na(matched_indices_odd))]
counts_odd <- counts[which(!is.na(matched_indices_odd)),];

matched_indices_even <- match(gene_names, chr_even_genes$external_gene_name)

gene_names_even <- gene_names[which(!is.na(matched_indices_even))]
counts_even <- counts[which(!is.na(matched_indices_even)),];

deng_fit_even <- maptpx::topics(t(counts_even), K=6, tol=100)
save(deng_fit_even, file="../rdas/deng_topic_fit_k_6_even_chromosome.rda")
deng_fit_odd <- maptpx::topics(t(counts_odd), K=6, tol=100)
save(deng_fit_odd, file="../rdas/deng_topic_fit_k_6_odd_chromosome.rda")

omega <- deng_fit_odd$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(meta_data$cell_type,
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
                                        "lateblast") ) ) )
rownames(omega) <- annotation$sample_id; 

StructureGGplot(omega = omega,
                annotation = annotation,
                figure_title = "Deng et al Structure Plot(odd chromosome genes)",
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))

omega <- deng_fit_even$omega

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(meta_data$cell_type,
                        levels = rev( c("zy", "early2cell",
                                        "mid2cell", "late2cell",
                                        "4cell", "8cell", "16cell",
                                        "earlyblast","midblast",
                                        "lateblast") ) ) )
rownames(omega) <- annotation$sample_id; 

StructureGGplot(omega = omega,
                annotation = annotation,
                figure_title = "Deng et al Structure Plot(even chromosome genes)",
                palette = RColorBrewer::brewer.pal(8, "Accent"),
                yaxis_label = "Development Phase",
                order_sample = TRUE,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))
