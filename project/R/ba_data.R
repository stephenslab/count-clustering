

############   BA  GTEx data  ########################

data <- data.frame(data.table::fread("../external_data/GTEX_V6/cis_gene_expression.txt"))
matdata <- data[,-(1:2)]

sample_names <- read.table("../external_data/GTEX_V6/samples_id.txt")

ba_indices <- grep("BA", sample_names[,3])

ba_data <- matdata[, ba_indices]

colnames(ba_data) <- sample_names[ba_indices,3]

out <- t(ba_data)

save(out, file="../external_data/GTEX_V6/ba_data_gtex.rda")
