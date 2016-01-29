
#

counts=read_counts

scale=1; K=4; N=500; G=100;

####  The starting value of omega (topic prop weights)  ####################

# Use a preset seed so the example is reproducable.
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
seed=54321))

omega_preprocess=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);

log_counts=log(counts+1e-07);
y=matrix(log_counts,1,N*G);


omega0=omega_preprocess
alpha0=matrix(rnorm((K)*G,1,1),nrow=(K));
param0=c(matrix(omega0,1,N*K),matrix(alpha0,1,K*G));

library(SQUAREM)

Poisson_preprocess <- function(param_vec,y)
{
	log_counts=matrix(y,N,G);
	#N=dim(log_counts)[1]; G=dim(log_counts)[2];
	omega0=matrix(param_vec[1:(N*K)],N,K);
	alpha0=matrix(param_vec[-(1:(N*K))],K,G);
	
	svd_omega=svd(omega0);
	temp1=t(svd_omega$v)%*%diag(1/svd_omega$d^2,dim(omega0)[2])%*%svd_omega$v;
	temp2=t(omega0)%*%log_counts;
	# temp2=t(omega0)%*%pnorm_counts;
	temp1=solve(t(omega0)%*%omega0);
	H = temp1%*%temp2;

	###  Estimation of the matrix W (or omega) 

	omega=matrix(0,N,K);
	for(n in 1:N)
	{
		omega_vec=omega0[n,];
		counts_vec=log_counts[n,];
		res=optim(reverse_transform(omega_vec), function(v) loglik_norm(transform(v),counts_vec,t(H)) );
		omega[n,]=transform(res$par);
	}

	param_vec_omega=matrix(omega,1,N*K);
	param_vec_alpha=matrix(H,1,K*G);
	param_new_vec=c(param_vec_omega,param_vec_alpha);
	param_vec=param_new_vec
	return(param_new_vec)
}

options(warn=-1)
system.time(res <- squarem(p=param0,y=y, fixptfn=Poisson_preprocess, control=list(maxiter = 100, trace = FALSE)));


################################################################################################

omega_initial=matrix(res$par[1:(N*K)],N,K);
alpha_initial=matrix(res$par[-(1:(N*K))],K,G)

docweights=omega_initial;
library(permute);
library("BioPhysConnectoR");
perm_set=rbind(1:K,allPerms(1:K));
diff=array(0,dim(perm_set)[1]);
for (p in 1:dim(perm_set)[1])
{
	temp=docweights[,perm_set[p,]];
	diff[p]=fnorm(temp,omega_true);
}

p_star=which(diff==min(diff));
docweights=docweights[,perm_set[p_star,]];

omega_initial=docweights;
alpha_initial=alpha_initial[perm_set[p_star,],];

windows()
k=K
barplot(t(omega_initial),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


