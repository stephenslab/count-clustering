##  A simulation mechanism to compare hierarchical and admixture methods of clustering

rm(list=ls())

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(maptpx)))
suppressMessages(suppressWarnings(library(gplots)))
suppressMessages(suppressWarnings(library(philentropy)))
suppressMessages(suppressWarnings(library(plyr)))

data= data.frame(fread('../data/gtex_data/gtex_thinned_version_1.txt'))[,-1];
samples_id=read.table("../data/gtex_data/samples_id.txt");

test_indices <- which(samples_id[,3]=='Heart - Left Ventricle' | samples_id[,3]=='Muscle - Skeletal');

N_run <- 200;

misclass_admix <- array(0, N_run);
misclass_hierarchy <- array(0, N_run);

for(run in 1:N_run)
{
  test_indices_sampled <- sample(test_indices, 200, replace=FALSE);

  samples <- data[,test_indices_sampled];

  tissue_ids_of_samples <- samples_id[test_indices_sampled,3];
  tissue_lab_hierarchy <- tissue_ids_of_samples[heatmap(distance(t(samples),method="euclidean"),hclustfun = function(x) hclust(x,method="average"))$rowInd];
  tissue_lab_hierarchy <- as.numeric(revalue(droplevels(as.factor(tissue_lab_hierarchy)), c('Heart - Left Ventricle'='0', 'Muscle - Skeletal'='1')))-1;


  Topic_Clus = topics(t(samples), K=2, tol=0.005);
  docweights_samples=Topic_Clus$omega;

  tissue_lab_admix <- tissue_ids_of_samples[heatmap(distance(docweights_samples,method="euclidean"),hclustfun = function(x) hclust(x,method="average"))$rowInd];
  tissue_lab_admix <- as.numeric(revalue(droplevels(as.factor(tissue_lab_admix)), c('Heart - Left Ventricle'='0', 'Muscle - Skeletal'='1')))-1;

  truth1 <- c(rep(0, length(which(tissue_lab_admix==0))), rep(1,length(which(tissue_lab_admix==1))))
  truth2 <- c(rep(1, length(which(tissue_lab_admix==1))), rep(0,length(which(tissue_lab_admix==0))))

  misclass_admix[run] <- min(length(which(truth1 !=tissue_lab_admix)), length(which(truth2 !=tissue_lab_admix)))/ length(tissue_lab_admix);
  misclass_hierarchy[run] <- min(length(which(truth1 !=tissue_lab_hierarchy)), length(which(truth2 !=tissue_lab_hierarchy)))/ length(tissue_lab_hierarchy);
  cat("Run number:", run);
}

write.table(misclass_admix,"../internal_data/misclass_admix_sim.txt");
write.table(misclass_hierarchy,"../internal_data/misclass_hierarchy_sim.txt");

misclass_admix <- as.vector(as.matrix(read.table('../internal_data/misclass_admix_sim.txt')));
misclass_hierarchy <- as.vector(as.matrix(read.table('../internal_data/misclass_hierarchy_sim.txt')));

plot(misclass_admix, misclass_hierarchy, lwd=2,pch=20, lty=1, ylab="",xlab="")
title(main="Scatter plot")
title(xlab="misclassification proportion (admixture)")
title(ylab="misclassification proportion (hierarchical)")

plot(density(misclass_hierarchy),col="red",ylim=c(0,3));
lines(density(misclass_admix),col="green",ylim=c(0,3))
