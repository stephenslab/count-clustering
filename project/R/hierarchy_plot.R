

#####################  Plot hierarchical maps  ###########################


hierarchy_prop_1 <- get(load("../rdas/hierarchy_prop_1.rda"))
hierarchy_prop_2 <- get(load("../rdas/hierarchy_prop_2.rda"))
hierarchy_prop_3 <- get(load("../rdas/hierarchy_prop_3.rda"))
hierarchy_prop_4 <- get(load("../rdas/hierarchy_prop_4.rda"))

samples_id <- read.table("../external_data/GTEX_V6/samples_id.txt")
tissue_labels <- samples_id[,3]

tt <- table(tissue_labels)
tissue_names <- names(tt)

tissues_to_consider <- tissue_names[which(as.numeric(table(tissue_labels))>60)]

colnames(hierarchy_prop_1) <- tissues_to_consider
rownames(hierarchy_prop_1) <- tissues_to_consider

hierarchy_F1 <- hierarchy_prop_1 + t(hierarchy_prop_1)

hierarchy_F1[hierarchy_F1 < 0.2] = 0
hierarchy_F1[hierarchy_F1 >= 0.2] = 1

cape::myImagePlot(1 - hierarchy_F1)


colnames(hierarchy_prop_2) <- tissues_to_consider
rownames(hierarchy_prop_2) <- tissues_to_consider

hierarchy_F2 <- hierarchy_prop_2 + t(hierarchy_prop_2)

hierarchy_F2[hierarchy_F2 < 0.2] = 0
hierarchy_F2[hierarchy_F2 >= 0.2] = 1

cape::myImagePlot(1 - hierarchy_F2)


colnames(hierarchy_prop_3) <- tissues_to_consider
rownames(hierarchy_prop_3) <- tissues_to_consider

hierarchy_F3 <- hierarchy_prop_3 + t(hierarchy_prop_3)

hierarchy_F3[hierarchy_F3 < 0.2] = 0
hierarchy_F3[hierarchy_F3 >= 0.2] = 1

cape::myImagePlot(1 - hierarchy_F3)

colnames(hierarchy_prop_4) <- tissues_to_consider
rownames(hierarchy_prop_4) <- tissues_to_consider

hierarchy_F4 <- hierarchy_prop_4 + t(hierarchy_prop_4)

hierarchy_F4[hierarchy_F4 < 0.2] = 0
hierarchy_F4[hierarchy_F4 >= 0.2] = 1

cape::myImagePlot(1 - hierarchy_F4)


