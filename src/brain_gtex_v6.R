###  Brain GTEx V6 data

data <- data.frame(fread('../data/GTEX_V6/cis_gene_expression.txt'));
matdata <- data[,-(1:2)];
samples_id=read.table("/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/samples_id.txt")[,3];
brain_indices <- grep("Brain", samples_id);

brain_data <- matdata[,brain_indices];
colnames(brain_data) <- samples_id[brain_indices];
brain_data_frame <- cbind.data.frame(data[,2],brain_data);

write.table(brain_data_frame, '../data/GTEX_V6/cis_gene_expression_brain.txt');

gene_names <- cbind.data.frame(data[,2]);
colnames(gene_names) <- "cis_gene_names"
write.table(gene_names, "../data/GTEX_V6/gene_names_GTEX_V6.txt");
