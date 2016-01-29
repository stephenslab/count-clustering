---
layout: page
title: "Model based clustering of RNA-seq data"
tagline:
---

### GTEX V6 analysis

  The Genotype Tissue Expression (GTEx) Project is a large scale project collecting tissue samples from actual human tissues. Check the [paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010069/).

  We apply model based clustering of the bulk-RNA GTEx V6 data (check the official release: [GTEx portal](http://www.gtexportal.org/home/)). The data contained $8555$ tissue samples coming from $53$ different tissues and we focussed on $16069$ genes chosen using a filtering criterion. Below, we present the cluster analysis for $K=15$ on all the tissues as well as clustering of brain samples using $K=4$.

  *[GTEX V6 analysis](project/src/gtex_v6_structure_genes.html)

### Jaitin et al (2014) single cell analysis

  We apply the model based clustering on the dataset due to Jaitin \emph{et al} 2014, [paper] (http://science.sciencemag.org/content/343/6172/776).

  Jaitin \textit{et al} sequenced over $4000$ single cells from mouse spleen. Following the original authors protocol, we also filtered out 16 genes that they found to show significant batch-specific expression. Here we analyze 1041 of these cells that were categorized as $CD11c+$ in the \textit{sorting markers} column of their data ([link](http://compgenomics.weizmann.ac.il/tanay/?page_id=519)), and which had total number of reads mapping to non-ERCC genes greater than $600$. (We believe these cells correspond roughly to the 1040 cells in their Figure S7.)

  Below we fit our model on the above data for $K=7$.

  *[Jaitin et al 2014 sc-RNA analysis](project/src/jaitin_structure_genes.html)

### Deng et al (2014) single cell analysis

  We apply the model based clustering on the dataset due to Deng \emph{et al} 2014, [paper](http://science.sciencemag.org/content/343/6167/193).

  Deng \textit{et al} collected expression data from individual cells from zygote to blastocyst stages of mouse preimplantation development \cite{Deng2014}. Deng \textit{et al}'s analysis focussed particularly on allele-specific expression from the two contributing mouse strains (CAST/EiJ and C57BL/6J).

  Here we present the topic model fit and Structure plot on this data for a range of values of $K$ to see how patterns change with increasing number of topics.

  *[Deng et al 2014 sc-RNA analysis](project/src/deng_structure_all_genes.html)


