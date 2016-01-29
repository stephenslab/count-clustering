
rm(list=ls())
setwd('/Users/kushal/Documents/count-clustering/src/')
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(maptpx)))
suppressMessages(suppressWarnings(library(gplots)))
suppressMessages(suppressWarnings(library(philentropy)))
suppressMessages(suppressWarnings(library(plyr)))

data= data.frame(fread('../data/gtex_data/gtex_thinned_version_1.txt'))[,-1];
samples_id=read.table("../data/gtex_data/samples_id.txt");


unique_tissues <- unique(samples_id[,3]);

hierarchy_mat <- matrix(0,length(unique_tissues), length(unique_tissues));
admixture_mat <- matrix(0,length(unique_tissues), length(unique_tissues));

thinning_coef <- 1
clus <- 2

unique_tissues <- unique_tissues[-which(as.numeric(table(samples_id[,3])) < 50)];

#data <- apply(data,c(1,2), function(x) rbinom(1,x,thinning_coef));

for(i in 2:length(unique_tissues))
{
  for(j in 1:(i-1))
  {
    test_indices <- which(samples_id[,3]==toString(unique_tissues[i]) | samples_id[,3]==toString(unique_tissues[j]) );
    test_indices_sampled <- sample(test_indices, 50, replace=FALSE);
    sample_data <- t(data[,test_indices_sampled]);
    #sample_data <- apply(sample_data,c(1,2),function(x) rbinom(1,x,thinning_coef));
    samples_id_sub <- droplevels(samples_id[test_indices_sampled,3]);
   # samples_id_sub <- revalue(samples_id_sub,c(paste0(levels(samples_id_sub)[1])='0', paste0(levels(samples_id_sub)[2])='1')
  #  samples_id_sub <- revalue(samples_id_sub,c('Adipose - Subcutaneous'='0', 'Adipose - Visceral (Omentum)'='1')
    
    hc <- hclust(dist(sample_data), "complete");
    clus_lab <- cutree(hc, k = clus)
    class_lab <- samples_id_sub;
    freq_tabs <- xtabs(~class_lab + clus_lab);
    recall_mat <- t(apply(freq_tabs, 1, function(x) x/sum(x)));
    prec_mat <- apply(freq_tabs, 2, function(x) x/sum(x));
    F_mat <- 2/(1/prec_mat + 1/recall_mat);
    F_score <- (rowSums(freq_tabs)%*%apply(F_mat,1, function(x) max(x)))/sum(freq_tabs);
    
    Topic_Clus = topics(sample_data, K=clus, tol=0.0005);
    docweights_samples=Topic_Clus$omega;
    hc_topic <- hclust(dist(docweights_samples), "complete")
    clus_lab <- cutree(hc_topic, k = clus)
    class_lab <- samples_id_sub;
    freq_tabs <- xtabs(~class_lab + clus_lab);
    recall_mat <- t(apply(freq_tabs, 1, function(x) x/sum(x)));
    prec_mat <- apply(freq_tabs, 2, function(x) x/sum(x));
    F_mat_topic <- 2/(1/prec_mat + 1/recall_mat);
    F_score_topic <- (rowSums(freq_tabs)%*%apply(F_mat_topic,1, function(x) max(x)))/sum(freq_tabs);
    
    hierarchy_mat[i,j] <- F_score;
    admixture_mat[i,j] <- F_score_topic;
  }
  
  cat("We are at tissue: ", i)
}


    
    
    

