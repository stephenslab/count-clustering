
pc_data <- get(load("../rdas/prcomp_gtex_v6.rda"));
pc_data_x <- pc_data$x;

dim(pc_data_x)
pc_data_x_filtered  <- pc_data_x[,1:10]
save(pc_data_x_filtered, file="../rdas/pc_gtex_v6_x.rda")
