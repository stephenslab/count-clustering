

###  SFA on Deng et al (2014) data

library(devtools)

read.data1 = function() {
  x = tempfile()
  download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
  z = get(load((x)))
  return(z)
}

Deng2014MouseESC <- read.data1()

deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

deng.counts[1:3,1:3]

voom_deng_counts <- limma::voom(t(deng.counts));
sqrt_deng_counts <- sqrt(t(deng.counts));

write.table(t(deng.counts), file="../sfa_inputs/deng_counts", 
            row.names=FALSE,col.names=FALSE, quote=FALSE, sep="\t")

write.table(voom_deng_counts, file="../sfa_inputs/voom_deng_counts", 
            row.names=FALSE,col.names=FALSE, quote=FALSE, sep="\t")

write.table(sqrt_deng_counts, file="../sfa_inputs/sqrt_deng_counts", 
            row.names=FALSE,col.names=FALSE, quote=FALSE, sep="\t")

### Go to the sfa directory

#./sfa_mac -gen ../translog_data.txt -g 259 -n 22431 -k 10 -iter 50 -r 800 -mn -o output/sfa_alex

#./sfa_mac -gen ../../sfa_inputs/sqrt_deng_counts -g 259 -n 22431 -k 10 -iter 50 -r 800 -mn -mg -o ../../sfa_outputs/Deng2014/sqrt_deng_counts
#./sfa_mac -gen ../../sfa_inputs/sqrt_deng_counts -g 259 -n 22431 -k 10 -iter 50 -r 800 -mn -mg -o ../../sfa_outputs/Deng2014/sqrt_deng_counts
#./sfa_mac -gen ../../sfa_inputs/sqrt_deng_counts -g 22431 -n 259 -k 10 -iter 50 -r 800 -mn -mg -t  -o ../../sfa_outputs/Deng2014/sqrt_deng_counts_transpose
# ./sfa_linux -gen ../../sfa_inputs/sqrt_deng_counts -g 22431 -n 259 -k 10 -iter 100 -r 500 -mn -mg -t -o ../../sfa_outputs/Deng2014/sqrt_deng_counts_transpose

#./sfa_mac -gen ../../sfa_inputs/sqrt_deng_counts -g 22431 -n 259 -k 10 -iter 50 -r 800 -mn -mg -t  -o ../../sfa_outputs/Deng2014/sqrt_deng_counts_transpose

# ./SFAmix --nf 10 --y ../sfa_inputs/sqrt_deng_counts  --out ../sfa_outputs/Deng2014SFAmix --sep tab


lambda_out <- read.table("../sfa_outputs/Deng2014/sqrt_deng_counts_transpose_lambda.out")
f_out <- read.table("../sfa_outputs/Deng2014/sqrt_deng_counts_transpose_F.out")

par(mar=c(10,4,2,2))

for(k in 1:10){
barplot(t(lambda_out[,k]), axisnames=F, space=0, border=NA,
        main="",las=1,ylim=c(min(lambda_out[,k]),max(lambda_out[,k])),
        cex.axis=1,cex.main=0.5)
title(main=paste("Taxonomic Factor Loadings, loading:",k), cex.main=0.7)
match_labs=match(unique(deng.meta_data[,1]),deng.meta_data[,1]);
match_labs_suffix=c(match_labs[2:length(match_labs)]-1,dim(lambda_out)[1]);
match_labs_prefix=match_labs[1:(length(match_labs))];
labs=match_labs_prefix + 0.5*(match_labs_suffix - match_labs_prefix);

axis(1,at=labs,unique(deng.meta_data[,1]),las=2, cex.axis=0.5);
abline(v=match_labs_suffix)
}

par(mar=c(10,4,2,2))

for(k in 1:10){
  barplot(t(f_out)[,k], axisnames=F, space=0, border=NA,
          main="",las=1,ylim=c(min(f_out[k,]),max(f_out[k,])),
          cex.axis=1,cex.main=0.5, xaxt="n")
  title(main=paste("Taxonomic Factor Loadings, loading:",k), cex.main=0.7)
  match_labs=match(unique(deng.meta_data[,1]),deng.meta_data[,1]);
  match_labs_suffix=c(match_labs[2:length(match_labs)]-1,dim(f_out)[2]);
  match_labs_prefix=match_labs[1:(length(match_labs))];
  labs=match_labs_prefix + 0.5*(match_labs_suffix - match_labs_prefix);
  
  axis(1,at=labs,unique(deng.meta_data[,1]),las=2, cex.axis=0.5);
  abline(v=match_labs_suffix)
}

