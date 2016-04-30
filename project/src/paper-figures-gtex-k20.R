## GTEx not thinned

#gom_model_fit <- get(load("../external_data/GTEX_V6/gtexv6fit.k.20.master.rda"))
gom_model_fit <- get(load("../../../../Downloads/gtexv6fit.k.20.master.rda"))
omega <- gom_model_fit$omega
colnames(omega) <- c(1:NCOL(omega))

sample_labels <- read.table("../external_data/GTEX_V6/samples_id.txt",
                            header = TRUE, sep = " ",
                            stringsAsFactors = FALSE)
tissue_labels <- vector("numeric", NROW(sample_labels))
tissue_labels <- sample_labels[ ,3]

tissue_levels <- unique(tissue_labels);

tissue_labels[grep("Gastroe", tissue_labels)] = "Esophagus -Gastroesophageal Jn."
tissue_labels[grep("cingulate", tissue_labels)] = "Brain - Anterior cortex (BA24)"
tissue_labels[grep("EBV", tissue_labels)] = "Cells -EBV-lymphocytes"
tissue_labels[grep("Suprapubic", tissue_labels)] = "Skin - Unexposed (Suprapubic)"
tissue_labels[grep("Lower Leg", tissue_labels)] = "Skin - Exposed (Lower Leg)"
tissue_labels[grep("accumbens", tissue_labels)] = "Brain - N. accumbens"
tissue_labels[grep("Caudate", tissue_labels)] = "Brain - Caudate"
tissue_labels[grep("Putamen", tissue_labels)] = "Brain - Putamen"


factor_levels <- c(
    "Brain - Spinal cord (cervical c-1)",
    "Brain - Substantia nigra",
    "Brain - Amygdala",
    "Brain - Hypothalamus",
    "Brain - Hippocampus",
    "Brain - Putamen",
    "Brain - Caudate",
    "Brain - N. accumbens",
    "Brain - Frontal Cortex (BA9)",
    "Brain - Anterior cortex (BA24)",
    "Brain - Cortex",
    "Brain - Cerebellum",
    "Brain - Cerebellar Hemisphere",
    "Cells - Transformed fibroblasts",
    "Cells -EBV-lymphocytes",
    "Spleen",
    "Whole Blood",
    "Muscle - Skeletal",
    "Liver",
    "Pancreas",
    "Stomach",
    "Kidney - Cortex",
    "Adrenal Gland",
    "Colon - Transverse",
    "Small Intestine - Terminal Ileum",
    "Heart - Atrial Appendage",
    "Heart - Left Ventricle",
    "Minor Salivary Gland",
    "Skin - Sun Exposed (Lower leg)",
    "Skin - Unexposed (Suprapubic)",
    "Lung",
    "Ovary",
    "Thyroid",
    "Pituitary",
    "Testis",
    "Nerve - Tibial",
    "Breast - Mammary Tissue",
    "Adipose - Visceral (Omentum)",
    "Adipose - Subcutaneous",
    "Artery - Coronary",
    "Artery - Tibial",
    "Artery - Aorta",
    "Esophagus - Mucosa",
    "Vagina",
    "Cervix - Endocervix",
    "Esophagus -Gastroesophageal Jn.",
    "Colon - Sigmoid",
    "Esophagus - Muscularis",
    "Cervix - Ectocervix",
    "Fallopian Tube",
    "Prostate",
    "Uterus",
    "Bladder")


annotation <- data.frame(
    sample_id = paste0("X", 1:length(tissue_labels)),
    tissue_label = factor(tissue_labels,
                          levels = rev(factor_levels) ) )

cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Paired"))[c(3,4,7,8,11,12,5,6,9,10)],
           RColorBrewer::brewer.pal(12, "Set3")[c(1,2,5,8,9)],
           RColorBrewer::brewer.pal(9, "Set1")[c(9,7)],
           RColorBrewer::brewer.pal(8, "Dark2")[c(3,4,8)])


# pdf("plots/gtex-figures/gtex-k20-thin01.pdf",
#     height = 9, width=4)
CountClust::StructureGGplot(omega = omega,
                            annotation= annotation,
                            palette = cols1,
                            yaxis_label = "",
                            order_sample = TRUE,
                            split_line = list(split_lwd = .1,
                                              split_col = "white"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 5,
                                             axis_label_face="bold"))
#dev.off()



### Thinned data GTEX V6

gom_model_fit <- get(load("../external_data/GTEX_V6/gtexv6fit.k20.thin.0001.rda"))
#gom_model_fit <- get(load("../../../../Downloads/gtexv6fit.k20.thin.0001.rda"))
omega_thin <- gom_model_fit$omega


cols2 <- array(0, length(cols1));
cols2[3] <- cols1[2]
cols2[9] <- cols1[9]
cols2[7] <- cols1[8]
cols2[6] <- cols1[12]
cols2[16] <- cols1[19]
cols2[15] <- cols1[17]
cols2[8] <- cols1[7]
cols2[18] <- cols1[20]
cols2[17] <- cols1[18]
cols2[10] <- cols1[16]
cols2[11] <- cols1[14]
cols2[5] <- cols1[6]
cols2[12] <- cols1[11]
cols2[4] <- cols1[3]
cols2[14] <- cols1[15]
cols2[2] <- cols1[4]
cols2[which(cols2=="0")] <- setdiff(cols1, cols2)

# pdf("plots/gtex-figures/gtex-k20-thin0001.pdf",
#     height = 9, width=4)
CountClust::StructureGGplot(omega = omega_thin,
                            annotation= annotation,
                            palette = cols2,
                            yaxis_label = "",
                            order_sample = TRUE,
                            split_line = list(split_lwd = .1,
                                              split_col = "white"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 5,
                                             axis_label_face="bold"))
#dev.off()



gom_model_fit <- get(load("../external_data/GTEX_V6/gtexv6fit.k20.thin.01.rda"))
#gom_model_fit <- get(load("../../../../Downloads/gtexv6fit.k20.thin.01.rda"))
omega_thin <- gom_model_fit$omega;

# CountClust::StructureGGplot(omega = omega_thin,
#                             annotation= annotation,
#                             palette = cols1,
#                             yaxis_label = "",
#                             order_sample = TRUE,
#                             split_line = list(split_lwd = .1,
#                                               split_col = "white"),
#                             axis_tick = list(axis_ticks_length = .1,
#                                              axis_ticks_lwd_y = .1,
#                                              axis_ticks_lwd_x = .1,
#                                              axis_label_size = 5,
#                                              axis_label_face="bold"))

cols3 <- array(0, length(cols1));
cols3[2] <- cols1[9]
cols3[6] <- cols1[2]
cols3[8] <- cols1[8]
cols3[10] <- cols1[12]
cols3[16] <- cols1[19]
cols3[15] <- cols1[17]
cols3[9] <- cols1[7]
cols3[20] <- cols1[20]
cols3[18] <- cols1[18]
cols3[14] <- cols1[10]
cols3[12] <- cols1[14]
cols3[7] <- cols1[6]
cols3[11] <- cols1[16]
cols3[3] <- cols1[3]
cols3[5] <- cols1[5]
cols3[13] <- cols1[15]
cols3[4] <- cols1[4]
cols3[which(cols3=="0")] <- setdiff(cols1, cols3)

# pdf("plots/gtex-figures/gtex-k20-thin01.pdf",
#     height = 9, width=4)
CountClust::StructureGGplot(omega = omega_thin,
                            annotation= annotation,
                            palette = cols3,
                            yaxis_label = "",
                            order_sample = TRUE,
                            split_line = list(split_lwd = .1,
                                              split_col = "white"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 5,
                                             axis_label_face="bold"))
#dev.off()





######<---- brain tissue plots

omega <- read.table("../rdas/omega_cis_genes_brain.txt",
#omega <- read.table("../../../../Downloads/omega_cis_genes_brain.txt",
                    header = TRUE, sep = " ",
                    stringsAsFactors = FALSE)
dim(omega)
colnames(omega) <- c(1:NCOL(omega))
head(omega)

# make cell sample labels
# want a version consistent with majority of the literature
sample_labels <- read.table("../rdas/samples_id.txt",
                            header = TRUE, sep = " ",
                            stringsAsFactors = FALSE)

tissue_labels <- vector("numeric", NROW(sample_labels))
tissue_labels <- sample_labels[ ,3]

tissue_levels <- unique(tissue_labels);

tissue_labels[grep("accumbens", tissue_labels)] = "Brain - N. accumbens"
tissue_labels[grep("Caudate", tissue_labels)] = "Brain - Caudate"
tissue_labels[grep("Putamen", tissue_labels)] = "Brain - Putamen"
tissue_labels[grep("cingulate", tissue_labels)] = "Brain - Anterior cortex (BA24)"

brain_labels <- tissue_labels[grep("Brain", tissue_labels)]

# assign tissue labels
rownames(omega) <- paste0("X", 1:length(brain_labels))
annotation <- data.frame(
    sample_id = paste0("X", 1:length(brain_labels)),
    tissue_label = factor(brain_labels,
                          levels = rev(c("Brain - Spinal cord (cervical c-1)",
                                         "Brain - Substantia nigra",
                                         "Brain - Amygdala",
                                         "Brain - Hypothalamus",
                                         "Brain - Hippocampus",
                                         "Brain - Putamen",
                                         "Brain - Caudate",
                                         "Brain - N. accumbens",
                                         "Brain - Frontal Cortex (BA9)",
                                         "Brain - Anterior cortex (BA24)",
                                         "Brain - Cortex",
                                         "Brain - Cerebellum",
                                         "Brain - Cerebellar Hemisphere") ) ) )

# define colors of the clusers
cols <- c("blue", "darkgoldenrod1", "cyan", "red")

# pdf("plots/gtex-figures/brain-barplot.pdf",
#     height = 4, width = 4)
CountClust::StructureGGplot(omega = omega,
                            annotation= annotation,
                            palette = cols,
                            yaxis_label = "",
                            order_sample = TRUE,
                            split_line = list(split_lwd = .4,
                                              split_col = "white"),
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))
#dev.off()