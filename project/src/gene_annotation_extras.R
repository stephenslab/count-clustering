### Shallow green cluster (Liver)

```{r echo=FALSE, eval=TRUE}

gene_names_mat_poisson[12,]
lapply(1:10, function(n) out[grep(gene_names_mat_poisson[12,n], out$query),])

```

### Yellow cluster (Cells transformed fibroblasts)

```{r echo=FALSE, eval=TRUE}

gene_names_mat_poisson[9,]
lapply(1:10, function(n) out[grep(gene_names_mat_poisson[9,n], out$query),])

```

### Shady Blue cluster (Arteries)

```{r echo=FALSE, eval=TRUE}

gene_names_mat_poisson[3,]
lapply(1:10, function(n) out[grep(gene_names_mat_poisson[3,n], out$query),])

```

### Blue cluster (Brain)

```{r echo=FALSE, eval=TRUE}

gene_names_mat_poisson[2,]
lapply(1:10, function(n) out[grep(gene_names_mat_poisson[2,n], out$query),])

```

One thing we find is that we are getting a gene ENSG00000244734 in multiple tissue clusters. It is a hemoglobin related gene. We look at the topic expression for that gene.

```{r echo=FALSE, eval=TRUE}

colnames(topics_theta)=paste("V_",1:12);
topics_theta[grep("ENSG00000244734", gene_names),]

```

Note that this has very high expression in the 10th cluster (grey cluster) which corresponds to blood, but since it is divergent from other tissues, other tissues are also reporting it. Is that what we want??

We also look at the Brain sub tissues and cluster them separately. Here Brain cerebellum and cerebellur hemisphere are pretty different from the basal ganglia and hippocampus etc, that is because cerebellum and cerebellur hemisphere contain very high percentage of neurons compared to other parts. Does that egt reflected in the clusters we get? What does gene annotations say?


## Brain samples clustering (Structure)

```{r brain_structure, echo=FALSE, eval=TRUE}

samples_id=read.table("/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/samples_id.txt")

brain_indices=which(samples_id[,2]=='Brain');
brain_subgroups_id =samples_id[brain_indices,3];

docweights= read.table('/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/topics_brain_thinned_omega_clus_4.txt');

K=4
barplot(t(docweights),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
brain_subgroups_id[394]=brain_subgroups_id[395];
labels=match(unique(brain_subgroups_id), brain_subgroups_id);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(docweights)[1]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(brain_subgroups_id),las=2);


```


## Brian topics expression data 

```{r brain_expr, echo=FALSE, eval=TRUE}


data <- data.frame(fread("/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/gtex_thinned_version_1.txt"));

brain_samples=data[,brain_indices];

#Topic_Clus <- topics(t(brain_samples),K=4, tol=0.001);

#write.table(Topic_Clus$omega, '/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/topics_brain_thinned_omega_clus_4.txt');


#write.table(Topic_Clus$theta, '/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/topics_brain_thinned_theta_clus_4.txt');

gene_snp_names_gtex <- data.frame(fread('/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/gene_snp_names.txt',header=FALSE));
gene_names <- substring(gene_snp_names_gtex[,2],1,15);


topics_theta <- read.table('/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/topics_brain_thinned_theta_clus_4.txt');

features <- ExtractTopFeatures(topics_theta,top_features = 20, method="poisson")
features_vec <- unique(as.vector(features));

class <- as.numeric(apply(topics_theta[features_vec,], 1, which.max))

imp_gene_names <- gene_names[features_vec];

imp_genes_per_class <- lapply(1:dim(topics_theta)[2], function(x) imp_gene_names[which(class==x)]);

```

We use the KL Poisson model to get a set of genes that were highly expressed in one cluster compared to another.

## Green cluster (Cerebellum and Cerebellur Hemisphere)

```{r echo=FALSE, eval=TRUE}

out <- queryMany(imp_genes_per_class[[3]], scopes="ensembl.gene", fields=c("name", "summary"), species="human");

out

```

## Blue cluster (Spinal cord + Substantia nigra)

```{r echo=FALSE, eval=TRUE}

out <- queryMany(imp_genes_per_class[[4]], scopes="ensembl.gene", fields=c("name", "summary"), species="human");

out
```


## Red cluster (Ganglia?)

```{r echo=FALSE, eval=TRUE}

out <- queryMany(imp_genes_per_class[[1]], scopes="ensembl.gene", fields=c("name", "summary"), species="human");

out
```


## green cluster (cerebellum)

```{r echo=FALSE, eval=TRUE}

out <- queryMany(imp_genes_per_class[[2]], scopes="ensembl.gene", fields=c("name", "summary"), species="human");

out
```

<!-- Two genes that seem to crop up a numebr of times in the above analysis are **ENSG00000197971** and **ENSG00000132639**.

Let us see how the topic expressions for these two genes look like

**ENSG00000197971** -->
  
  
  
  ```{r, echo=FALSE, eval=FALSE}

topics_theta[grep("ENSG00000197971", gene_names),]
```

```{r, echo=FALSE, eval=FALSE}

topics_theta[grep("ENSG00000132639", gene_names),]

```




```{r, echo=FALSE, eval=FALSE}

## Bernoulli K-L divergence comparison

indices_mat_bernoulli=matrix(0,dim(topics_theta)[2],20);

for(k in 1:dim(topics_theta)[2])
{
  temp_mat <- KL_score_bernoulli[[k]][,-k];
  vec <- apply(temp_mat, 1, function(x) max(abs(x)) )
  indices_mat_bernoulli[k,] = order(vec, decreasing = TRUE)[1:20]
}

significant_indices_bernoulli = unique(as.vector(indices_mat_bernoulli));

topics_theta[significant_indices_bernoulli,];

#gene_names_eberwine=data.frame(fread("/Users/kushal/Documents/Matthew Stephens Project/counts_clustering/counts_data.txt", header=FALSE))[,2];

gene_snp_names_gtex <- data.frame(fread('/Users/kushal/Documents/gtex-viz/gtex.Kushal/data/gene_snp_names.txt',header=FALSE));
gene_names <- substring(gene_snp_names_gtex[,2],1,15);

gene_names_sig_bernoulli = gene_names[significant_indices_bernoulli];

gene_names_mat_bernoulli = matrix(gene_names[as.vector(indices_mat_bernoulli)],nrow=dim(topics_theta)[2]);

color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");

```





```{r, echo=FALSE, eval=FALSE}

## sqrt z score comparison

temp_mat <- sqrt_zscore;

vec <- apply(temp_mat, 1, function(x) max(abs(x)) )
indices_mat_sqrt_zscore = order(vec, decreasing = TRUE)[1:500];

gene_names_sig_zscore = gene_names_eberwine[indices_mat_sqrt_zscore];

apply(topics_theta[indices_mat_sqrt_zscore,],1,which.max)


```

