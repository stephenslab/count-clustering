##  A simulation mechanism to compare hierarchical and admixture methods of clustering

rm(list=ls())
setwd('/Users/kushal/Documents/count-clustering/src/')
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

log2 <- function(x)
{
  y <- matrix(0, nrow(x), ncol(x));
  
  for(m in 1:nrow(x))
  {
    for(n in 1:ncol(x))
    {
      if(x[m,n]<=0) {
        y[m,n]=0;
      }else {
        y[m,n] <- log(x[m,n]);
     }    
    }
  }
  return(y)
}

for(run in 1:N_run)
{
  test_indices_sampled <- sample(test_indices, 50, replace=FALSE);
  samples_id_sub <- droplevels(samples_id[test_indices_sampled,3]);

  ##  assign  numeric class labels
  
  samples_id_sub <- revalue(samples_id_sub,c('Heart - Left Ventricle'='0', 'Muscle - Skeletal'='1'))
  sample_data <- t(data[,test_indices_sampled]);
  hc <- hclust(dist(sample_data), "ave");
  nclus <- 11
  
  F_score <- array(0,nclus-1);
  E_score <- array(0,nclus-1);
  
  for(clus in 2:nclus)
  {
    clus_lab <- cutree(hc, k = clus)
    class_lab <- samples_id_sub;
    freq_tabs <- xtabs(~class_lab + clus_lab);
    recall_mat <- t(apply(freq_tabs, 1, function(x) x/sum(x)));
    prec_mat <- apply(freq_tabs, 2, function(x) x/sum(x));
    F_mat <- 2/(1/prec_mat + 1/recall_mat);
    F_score[clus-1] <- (rowSums(freq_tabs)%*%apply(F_mat,1, function(x) max(x)))/sum(freq_tabs);
    E_score_clus <- - colSums(prec_mat * log2(prec_mat));
    E_score[clus-1] <- (colSums(freq_tabs)%*%E_score_clus)/sum(freq_tabs);
  }
  
  col = c(rgb(seq(0,1,length=15),1,seq(0,1,length=15)), rgb(1,seq(1,0,length=15),seq(1,0,length=15)));
  
  png(filename = "../plots/heart_muscle_hierarchical_heatmap_average.png")
  heatmap.2(as.matrix(dist(sample_data)),
            labCol=samples_id[test_indices_sampled,3],
            labRow=samples_id[test_indices_sampled,3],
            scale="none", trace="none", 
            distfun=function(x) dist(x,method="euclidean"), 
            col=col, hclustfun = function(x) hclust(x,method="average"));
  dev.off()
  
  num_topics <- c(2,3,4,5);
  
  E_score_topics <- matrix(0,length(num_topics),nclus-1);
  F_score_topics <- matrix(0,length(num_topics),nclus-1);

  for(num in 1:length(num_topics))
  {
    Topic_Clus = topics(sample_data, K=num_topics[num], tol=0.005);
    docweights_samples=Topic_Clus$omega;
    png(filename = "../plots/heart_muscle_admix_heatmap_average.png")
    heatmap.2(as.matrix(dist(docweights_samples)),
              labCol=samples_id[test_indices_sampled,3],
              labRow=samples_id[test_indices_sampled,3],
              scale="none", trace="none", 
              distfun=function(x) dist(x,method="euclidean"), 
              col=col, hclustfun = function(x) hclust(x,method="average"));
    dev.off()
    
    
    hc_topic <- hclust(dist(docweights_samples), "ave")
    for(clus in 2:nclus)
    {
      clus_lab <- cutree(hc_topic, k = clus)
      class_lab <- samples_id_sub;
      freq_tabs <- xtabs(~class_lab + clus_lab);
      recall_mat <- t(apply(freq_tabs, 1, function(x) x/sum(x)));
      prec_mat <- apply(freq_tabs, 2, function(x) x/sum(x));
      F_mat <- 2/(1/prec_mat + 1/recall_mat);
      F_score_topics[num,clus-1] <- (rowSums(freq_tabs)%*%apply(F_mat,1, function(x) max(x)))/sum(freq_tabs);
      E_score_topics_clus <- - colSums(prec_mat * log2(prec_mat));
      E_score_topics[num,clus-1] <- (colSums(freq_tabs)%*%E_score_topics_clus)/sum(freq_tabs);
    }
  }
  
  F_score.frame <- rbind(F_score,F_score_topics);
  
  matplot(t(F_score.frame), type="l", ylim=c(0.5,1),
          ylab="F score", xlab="tree cuts",
          main="F score comparison",lwd=2, lty=1)
  legend("topright",c("hierarchical","admix-2","admix-3","admix-4","admix-5"),
         col=1:5,lty=rep(1,5),lwd=rep(2,5), cex=0.7)
  
  d2 <- cophenetic(hc_topic)
  
  d1 <- cophenetic(hc)
  
  cor(d1,d2)
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
