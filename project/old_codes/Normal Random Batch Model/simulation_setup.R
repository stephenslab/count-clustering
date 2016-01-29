

############   Normal  Topic Model with Random Effect  ####################

###  We have solved the problem of normal topic model, now we introduce a random

###  effect in the normal topic model and see how the methods work.  We shall 

###  focus on the case with a single rando effect in the model 


##################   Set  up 1   for Counts  Table  ########################
K=4;
G=100;
N=500;

alpha_true=matrix(rnorm((K)*G,0.5,1),nrow=(K)); ### the matrix of fixed effects

Label.Batch=c(rep(1,N/4),rep(2,N/4),rep(3,N/4),rep(4,N/4));

B=max(Label.Batch);

sigmab_true=2;

beta_true=matrix(0,B,G);       ###  the matrix of the random effect

for(g in 1:G)
{
	beta_true[,g]=rnorm(B,mean=0,sd=sigmab_true);
}

over_dis=0.3;

noise_true=matrix(0,N,G);

for(n in 1:N)
{
	noise_true[n,]=over_dis*rnorm(G,0,1);
}

library(gtools)
T=10;
omega_true=matrix(rbind(rdirichlet(T*10,c(3,4,2,6)),rdirichlet(T*10,c(1,4,6,3)),
			rdirichlet(T*10,c(4,1,2,2)),rdirichlet(T*10,c(2,6,3,2)),
			rdirichlet(T*10,c(3,3,5,4))), nrow=N);


H_true=alpha_true;

###  generating the table 

read_counts=matrix(0,N,G);

for(n in 1:N)
{
	read_counts[n,]=omega_true[n,]%*%H_true +beta_true[Label.Batch[n],]+noise_true[n,];
}

windows()
k=K
barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)








