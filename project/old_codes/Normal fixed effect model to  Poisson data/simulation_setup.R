
#############   Simulation Set up  Poisson distribution ###################

############  Poisson Topic Model with Random Effect  #####################


 

K=4;
G=100;
N=100;

alpha_true=matrix(rnorm((K)*G,0.5),nrow=(K)); ### the matrix of fixed effects

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


H_true=alpha_true+1;
alpha_true=H_true

###  generating the table 


read_counts=matrix(0,N,G);

for(n in 1:N)
{
	for(g in 1:G)
	{

		mean=exp(omega_true[n,]%*%H_true[,g] +beta_true[Label.Batch[n],g]+0*rnorm(1,0,1));
		read_counts[n,g]=rpois(1,mean);
	}
}


	




