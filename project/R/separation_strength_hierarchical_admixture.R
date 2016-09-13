

#######  Hierarchical vs Admixture: Separation Strength  ###############


counts <- data.frame(data.table::fread("../external_data/GTEX_V6/cis_gene_expression.txt"))
matdata <- counts[,-(1:2)]
voom_out  <- limma::voom(matdata)
voom_weights <- voom_out$weights
voom_data <- t(voom_out$E)

samples_id <- read.table("../external_data/GTEX_V6/samples_id.txt")
tissue_labels <- samples_id[,3]
library(DESeq2)

library(edgeR)

cpm_data <- cpm(matdata, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)

tt <- table(tissue_labels)
tissue_names <- names(tt)

tissues_to_consider <- tissue_names[which(as.numeric(table(tissue_labels))>60)]

admix_prop <- matrix(0, length(tissues_to_consider), length(tissues_to_consider))
hierarchy_prop <- matrix(0, length(tissues_to_consider), length(tissues_to_consider))

library(maptpx)
library(slam)
source("../../../maptpx/R/count.R")
source("../../../maptpx/R/topics.R")
source("../../../maptpx/R/tpx.R")

for(m in 2:length(tissues_to_consider)){
  for(n in 1:(m-1)){
    tissue1 <- tissues_to_consider[m];
    tissue2 <- tissues_to_consider[n];
    index1 <- which(tissue_labels==tissue1)
    index2 <- which(tissue_labels==tissue2)
    
    matdata1 <- matdata[,index1[1:50]]
    matdata2 <- matdata[,index2[1:50]]
    pooled_data <- cbind(matdata1, matdata2)
    
    cpm_data1 <- cpm_data[,index1[1:50]]
    cpm_data2 <- cpm_data[,index2[1:50]]
    cpm_data_pooled <- cbind(cpm_data1, cpm_data2)
    
    factor1 <- c(rep(1,50), rep(2,50))
    dd <- dist(t(cpm_data_pooled))
    
    hclusters <- cutree(hclust(dd),2)
    tab <- xtabs(~hclusters + factor1)
    misclass1 <- (min(tab[1,1]+tab[2,2], tab[1,2]+ tab[2,1]))/100;
    cat("Hierarchical clustering: misclass", misclass1, "\n")
    
    suppressWarnings(topics_fit <- maptpx::topics(t(pooled_data), K=3, tol=100))
    tclusters <- apply(topics_fit$omega, 1, function(x) return(which.max(x)))
    levels(tclusters) <- c("1","2")
    tab2 <- xtabs(~tclusters + factor1)
    if(dim(tab2)[1]==2){
     misclass2 <- (min(tab2[1,1]+tab2[2,2], tab2[1,2]+ tab2[2,1]))/100;
    }else{
      misclass2 <- 0.5
    }
    cat("Admixture clustering: misclass", misclass2, "\n")
    
    admix_prop[m,n] <- misclass2;
    hierarchy_prop[m,n] <- misclass1;
    admix_prop[n,m] <- admix_prop[m,n]
    hierarchy_prop[n,m] <- hierarchy_prop[m,n]
  }
  cat("We are at tissue", m, "\n")
}

save(admix_prop, file="../rdas/admix_prop_Sep2016.rda")
save(hierarchy_prop, file="../rdas/hierarchy_prop_Sep2016.rda")

admix_prop <- get(load("../rdas/admix_prop_Sep2016.rda"))
hierarchy_prop <- get(load("../rdas/hierarchy_prop_Sep2016.rda"))

colnames(hierarchy_prop) <- tissues_to_consider
rownames(hierarchy_prop) <- tissues_to_consider
colnames(admix_prop) <- tissues_to_consider
rownames(admix_prop) <- tissues_to_consider

admix_prop[admix_prop==0.5] <- hierarchy_prop[admix_prop==0.5]


myImagePlot(as.matrix(hierarchy_prop))
myImagePlot(as.matrix(admix_prop))

hierarchy_F <- hierarchy_prop + t(hierarchy_prop)
admixture_F <- admix_prop + t(admix_prop)

hierarchy_F[hierarchy_F < 0.2] = 0
hierarchy_F[hierarchy_F >= 0.2] = 1

admixture_F[admixture_F < 0.2] = 0
admixture_F[admixture_F >= 0.2] = 1

cape::myImagePlot(1-as.matrix(hierarchy_F))
cape::myImagePlot(1-as.matrix(admixture_F))


admix_prop <- get(load("../rdas/admix_prop_Sep2016.rda"))
hierarchy_prop <- get(load("../rdas/hierarchy_prop_Sep2016.rda"))

admix_prop[admix_prop==0.5] <- hierarchy_prop[admix_prop==0.5]

admix_prop[upper.tri(admix_prop)] <- 0
hierarchy_prop[upper.tri(hierarchy_prop)] <- 0

admixture_F_thinned <- apply(admix_prop, c(1,2), function(x) (return(rbinom(1,1,prob=min(5*x,1)))))
hierarchy_F_thinned <- apply(hierarchy_prop, c(1,2), function(x) (return(rbinom(1,1,prob=min(5*x,1)))))

a1 <- data.frame(as.matrix(hierarchy_F_thinned + t(hierarchy_F_thinned)))
a2 <- data.frame(as.matrix(admixture_F_thinned + t(admixture_F_thinned)))

rownames(a1) <- tissues_to_consider
colnames(a1) <- tissues_to_consider
rownames(a2) <- tissues_to_consider
colnames(a2) <- tissues_to_consider

cape::myImagePlot(1-a1)
cape::myImagePlot(1-a2)




myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}