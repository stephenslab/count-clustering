---
layout: page
title: "Model based clustering of RNA-seq data"
tagline:
---

### GTEX V6 analysis

  The Genotype Tissue Expression (GTEx) Project is a large scale project collecting tissue samples from actual human tissues. Check the [paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010069/).

  We apply model based clustering of the bulk-RNA GTEx V6 data (check the official release: [GTEx portal](http://www.gtexportal.org/home/)). The data contained 8555 tissue samples coming from 53 different tissues and we focussed on 16069 genes chosen using a filtering criterion. Below, we present the cluster analysis and Structure plot visualizations for number of clusters K=15 on all the tissues as well as clustering of brain samples using K=4.

  *[GTEX V6 analysis](project/src/gtex_v6_structure_genes.html)
  
  An alternative representation of the clustering besides the Structure plot is t-SNE representation. We present the results from applying t-SNE on the original GTEx read counts data and also on the topic proportions from the topic model fit.
  
  *[GTEX V6: t-SNE representation](project/src/tissues_tSNE_2.html)
  
  We annotate the genes that drive the clusters. For each cluster, we find a set of few top genes that distinguish that cluster from the rest (we term these as cluster annotating genes).
  
  *[GTEX V6: Gene Annotations](project/src/gene_annotation_2.html)

### Jaitin et al (2014) single cell analysis

  We apply the model based clustering on the dataset due to Jaitin *et al* 2014, [paper](http://science.sciencemag.org/content/343/6172/776).

  Jaitin *et al* sequenced over 4000 single cells from mouse spleen. Following the original authors protocol, we also filtered out 16 genes that they found to show significant batch-specific expression. Here we analyze 1041 of these cells that were categorized as CD11c+ in the *sorting markers* column of their data ([link](http://compgenomics.weizmann.ac.il/tanay/?page_id=519)), and which had total number of reads mapping to non-ERCC genes greater than 600. (We believe these cells correspond roughly to the 1040 cells in their Figure S7.)

  Below we fit our model on the above data for K=7.

  *[Jaitin et al 2014 sc-RNA analysis](project/src/jaitin_structure_genes.html)

### Deng et al (2014) single cell analysis

  We apply the model based clustering on the dataset due to Deng *et al* 2014, [paper](http://science.sciencemag.org/content/343/6167/193).

  Deng *et al* collected expression data from individual cells from zygote to blastocyst stages of mouse preimplantation development. Deng *et al*'s analysis focussed particularly on allele-specific expression from the two contributing mouse strains (CAST/EiJ and C57BL/6J).

  Here we present the topic model fit and Structure plot on this data for a range of values of K (number of clusters) to see how patterns change with increasing number of topics.

  [Deng et al 2014 sc-RNA analysis](project/src/deng_structure_all_genes.html)
  
  The cluster annotations of the clusters for the data are given in the script below.
  
  [Deng et al 2014 cluster annotations](project/src/deng_cluster_annotations.html)

### Other applications of CountClust

  We apply a batch correction procedure `BatchCorrectedCounts()` to remove known technical effects using a voom type framework in the package `CountClust`. We present 3 simulation scenarios to present the effectiveness of the batch correction mechanism. 
  
  [Batch Correction scenarios using CountClust](project/src/batch_correction_scenarios.html)
  
  We also validated the clusters we obtained in GTEx V6 and the Deng 2014 data by considering some of the genes with interesting biological properties and checking the trends of log of expression across the different tissue samples for that gene.
  
  [Gene expression study of genes from cluster annotation of GTEx + Deng data](project/src/extracted_genes_expr_study.html)
  
<!---   It seems Testis and LCL cluster together which seems counter-intuitive. We show that the cluster corresponding to the two has representativeness from genes specific to both Testis and LCLs.
   
   [testis-lcl cluster analysis](project/src/lcl_testis_cluster_analysis.html)
-->
  
