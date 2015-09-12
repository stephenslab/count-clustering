
## create a chunk that performs one cell operation for a phase diagram

phase_plot_cell <- function(n.out, G, alpha, gamma)
{
  if(gamma > 2/G) stop('The value of gamma should be less than 2/(number of genes)');
  library(maptpx)
  omega_sim <- cbind(seq(alpha,1-alpha,length.out=n.out), 1- seq(alpha,1-alpha,length.out=n.out));
  K <- dim(omega_sim)[2];
  #barplot(t(omega_sim),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of   clusters=",K),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
  
  freq <- rbind(c(gamma,(2/G)-gamma,rep(1/G,G-2)),c((2/G)-gamma,gamma,rep(1/G,G-2)));
  counts <- t(do.call(cbind,lapply(1:dim(omega_sim)[1], function(x)  rmultinom(1,1000,prob=omega_sim[x,]%*%freq))));
  
  Topic_clus <- topics(counts, K=2,tol=0.001);
  K=2
  #barplot(t(Topic_clus$omega),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of   clusters=",K),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
  topics_theta <- Topic_clus$theta;
  K= dim(Topic_clus$omega)[2];
  
  KL_score_poisson <- lapply(1:K, function(n) 
  {
    out <- t(apply(topics_theta, 1, function(x)
    {
      y=x[n] *log(x[n]/x) + x - x[n];
      return(y)
    }));
    return(out)
  })
  
  
  indices_mat_poisson=matrix(0,dim(topics_theta)[2],10);
  
  for(k in 1:dim(topics_theta)[2])
  {
    temp_mat <- KL_score_poisson[[k]][,-k];
    vec <- order(temp_mat, decreasing=TRUE)[1:2]
    #  print(vec)
  }
  
  deviance <- min(mean((Topic_clus$omega -omega_sim)^2), mean((1-Topic_clus$omega -omega_sim)^2));
  flag=0;
  if(all(sort(vec)==c(1,2))) flag=flag+1;
  if(deviance < 0.1) flag=flag+2;
  return(flag)
}

n.out <- 200;
G <- 100;
phase_mat <- matrix(0,nrow=100,ncol=100);
flag1=1; flag2=1
for (alpha in seq(0,1,length.out=100)) 
{
  for(gamma in seq(0,2/G,length.out=100))
  {
    phase_mat[flag1,flag2] = phase_plot_cell(n.out,G,alpha,gamma)
    flag2=flag2+1;
  }
  flag1=flag1+1;
  flag2=1;
}

require(geoR)
phase_mat <- as.matrix(read.table('/Users/kushal/Documents/Matthew Stephens Project/counts_clustering/mat.phase.200.txt'));
#phase_mat1 <- as.matrix(read.table('/Users/kushal/Documents/Matthew Stephens Project/counts_clustering/mat.phase.200.txt'));

G=200;
out <- expand.grid(x=seq(0,1,length.out=100),y=seq(0,2/G,length.out=100))
phase_mat_vec <- matrix(phase_mat,nrow=100*100,byrow=TRUE);
plot(out[,1],out[,2],col=phase_mat_vec+1, 
     main=paste0("Phase diagram: no. of genes=",G),
     xlab="alpha (topic prop fluc)",
     ylab="gamma (freq fluc)")
legend("topleft",col=1:4,legend=c("no match","omega match","theta match","omega/theta match"),fill=1:4)
