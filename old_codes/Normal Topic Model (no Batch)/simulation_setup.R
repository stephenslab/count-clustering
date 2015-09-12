

############   Normal  Topic Model with Random Effect  ####################

###  We have solved the problem of normal topic model, now we introduce a random

###  effect in the normal topic model and see how the methods work.  We shall 

###  focus on the case with a single rando effect in the model 


##################   Set  up 1   for Counts  Table  ########################
K=4;
G=100;
N=100;

alpha_true=1000+matrix(rnorm((K)*G,0.5),nrow=(K)); ### the matrix of fixed effects

Label.Batch=c(rep(1,N/5),rep(2,N/5),rep(3,N/5),rep(4,N/5),rep(5,N/5));

B=max(Label.Batch);

sigmab_true=2;

beta_true=matrix(0,B,G);       ###  the matrix of the random effect

for(g in 1:G)
{
	beta_true[,g]=rnorm(B,mean=0,sd=sigmab_true);
}

library(gtools)
T=2;
omega_true=matrix(rbind(rdirichlet(T*10,c(3,4,2,6)),rdirichlet(T*10,c(1,4,6,3)),
			rdirichlet(T*10,c(4,1,2,2)),rdirichlet(T*10,c(2,6,3,2)),
			rdirichlet(T*10,c(3,3,5,4))), nrow=N);


H_true=alpha_true

###  generating the table 

read_counts=matrix(0,N,G);

for(n in 1:N)
{
	read_counts[n,]=omega_true[n,]%*%H_true +0*rnorm(G,0,1);
}



#############################################################################





