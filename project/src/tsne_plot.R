
## gtex t-SNe plot

library(qtlcharts)
library(tsne)
tsne_samples <- read.table("../internal_data/tsne_samples.txt");

sample_id_sub <- read.table("../internal_data/samples_id.txt")[,3];

options(warn=-1)
tsne_plot <- suppressWarnings(suppressMessages(iplot(tsne_samples[,1],tsne_samples[,2],as.numeric(sample_id_sub),sample_id_sub)))

htmlwidgets::saveWidget(tsne_plot, file="../plots/tsne_tissues.html",selfcontained=TRUE);


topicmodel_tsne_samples <- read.table("../internal_data/topicmodel_samples_tsne_15.txt");

sample_id_sub <- read.table("../internal_data/samples_id.txt")[,3];

options(warn=-1)
topicmodel_tsne_plot <- suppressWarnings(suppressMessages(iplot(topicmodel_tsne_samples[,1],topicmodel_tsne_samples[,2],as.numeric(sample_id_sub),sample_id_sub)))

htmlwidgets::saveWidget(topicmodel_tsne_plot, file="../plots/topicmodel_tsne_tissues.html",selfcontained=FALSE);
