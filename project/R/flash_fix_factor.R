

###########   fixed factor FLASH:  on Lindo data  #########################

signature_counts <- get(load("../../../ancient-damage/summary_data/signature-counts-clubbed-Lindo2016.rda"))

signature_set <- colnames(signature_counts)

topic_fit <- maptpx::topics(signature_counts, K=2);

fit_factors <- topic_fit$theta;
fit_omega <- topic_fit$omega;

library(flashr)


res <- signature_counts

K <- 2
theta_mat <- as.numeric()

for(k in 1:K){
    flash_fit <- flash(t(res),
                       tol=1e-5, maxiter_r1 = 100,
                       partype="known",
                       sigmae2_true = t(res+0.000001),
                       factor_value = fit_omega[,k], fix_factor = TRUE,
                       nonnegative=FALSE,
                       ash_para = list(control=list(maxiter=1000)))

     lfit <- flash_fit$l;
     theta_fit <- lfit/sum(lfit);
     theta_mat <- cbind(theta_mat, theta_fit)

     res <- res - fit_omega[,k]%*%t(theta_fit)
}

pdf("../../../ancient-damage/figures/lindo_flash_fixfactor_known_type.pdf")
par(mfrow=c(1,2))
for(k in 1:K){
plot(fit_factors[,k], col="red", type="l", xlab="signature index", ylab="factor")
lines(theta_mat[,k], col="blue")
legend("topright", fill=c("red", "blue"),
       legend=c("topicmodel", "flash_topicmodel"),
       box.lwd = 0.4, cex=0.5, lty=2,
       lwd=0.5)
title(paste("Factor:",k))
}
dev.off()



##########  fixed factor FLASH:  Deng et al  ##############################

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

deng_fit <- get(load("../rdas/deng_topic_fit.rda"))
deng_fit_clus6 <- deng_fit$clust_5

loadings <- deng_fit_clus6$omega
factors <- deng_fit_clus6$theta

res <- t(deng.counts)

K <- 6
theta_mat <- as.numeric()

for(k in 1:K){
  flash_fit <- flash(t(res),
                     tol=1e-5, maxiter_r1 = 30,
                     partype="known",
                     sigmae2_true = t(res+0.000001),
                     factor_value = loadings[,k], fix_factor = TRUE,
                     nonnegative=FALSE,
                     ash_para = list(control=list(maxiter=500)))

  lfit <- flash_fit$l;
  theta_fit <- lfit/sum(lfit);
  theta_mat <- cbind(theta_mat, theta_fit)

  res <- res - loadings[,k]%*%t(theta_mat[,k])
  cat("Finished step: ", k, "\n")
}

topic_fit <- list("omega"=loadings,
                  "theta"=theta_mat)
save(topic_fit, file="../rdas/topic_fit_flash_fixfactor_known.rda")

topic_fit <- get(load("../rdas/topic_fit_flash_fixfactor_known.rda"))

pdf("../plots/deng-figures/flash_deng_fixfactor_known.pdf")
par(mfrow=c(3,2))
for(k in 1:6){
plot(factors[,k], col="red", type="l", xlab="signature index", ylab="factor")
lines(topic_fit$theta[,k], col="blue")
legend("topright", fill=c("red", "blue"),
       legend=c("topicmodel", "flash_topicmodel"),
       box.lwd = 0.4, cex=0.5, lty=2,
       lwd=0.5)
title(paste("Factor:",k))
}
dev.off()


########################################################################






