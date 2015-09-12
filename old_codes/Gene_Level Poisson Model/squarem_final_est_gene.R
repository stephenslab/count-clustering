

##########  SQUAREM  Final Estimation Mechanism  ########################


counts=read_counts
y=matrix(counts,1,N*G);

scale=1; K=4; N=500; G=100;

####  The starting value of omega (topic prop weights)  ####################

# Use a preset seed so the example is reproducable.
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
seed=54321))

omega0=omega_initial;
alpha0=alpha_initial;

sigmab_0=0.3;

beta0=matrix(0,B,G);       ###  the matrix of the random effect

for(g in 1:G)
{
	beta0[,g]=rnorm(B,mean=0,sd=sigmab_0);
}


param0=c(matrix(omega0,1,N*K),matrix(alpha0,1,K*G),matrix(beta0,1,B*G));

library(lme4)

noise=noise_true;

counts_final.em <- function(param_vec,y)
{
	counts=matrix(y,N,G);
	omega0=matrix(param_vec[1:(N*K)],N,K);
	param_vec_hyp=param_vec[-(1:(N*K))];
	beta0=matrix(tail(param_vec,B*G),B,G);
	alpha0=matrix(param_vec_hyp[(1:(K*G))],K,G);

	################   Estimation of alpha  ###########################

	alpha=matrix(0,K,G);
	for(k in 1:K)
	{
		for(g in 1:G)
		{	
			temp1=0; temp2=0;
			for(n in 1:N)
			{
				temp1=temp1+omega0[n,k]*counts[n,g];
				temp2=temp2+omega0[n,k]*exp(beta0[Label.Batch[n],g]+noise[n,g]);
			}
			alpha[k,g]=log(temp1/temp2);
		}	
	}

	###############   Estimation of omega  ############################

	omega=matrix(0,N,K);
	Z=array(0,c(N,K,G));

	for(n in 1:N)
	{
		for(g in 1:G)
		{
			for(k in 1:K)
			{
				lambda=beta0[Label.Batch[n],g]+alpha[k,g]+noise[n,g];
				Z[n,k,g]=omega0[n,k]*dpois(counts[n,g],exp(lambda));
			}
			Z[n,,g]=(Z[n,,g]+1e-07)/sum(Z[n,,g]+1e-07)
		}
	}

	for(n in 1:N)
	{
		for(k in 1:K)
		{
			omega[n,k]=sum(Z[n,k,])/G;
		}
	}

	docweights=omega;
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
	omega=docweights[,perm_set[p_star,]];


	###############  Estimation of beta  #################################

	beta=matrix(0,B,G);

	for(b in 1:B)
	{
		for(g in 1:G)
		{
			temp1=0; temp2=0;

			for(n in 1:N)
			{
				if(Label.Batch[n]==b)
				{
					temp1=temp1+sum(omega[n,])*counts[n,g];
					temp2=temp2+(omega[n,]%*%exp(alpha[,g]+noise[n,g]));
				}
			}
			beta[b,g]=log(temp1/temp2);
		}
	}

	param_vec_omega=matrix(omega,1,N*K);
	param_vec_alpha=matrix(alpha,1,K*G);
	param_vec_beta=matrix(beta,1,B*G);
	param_new_vec=c(param_vec_omega,param_vec_alpha,param_vec_beta);
	param_vec=param_new_vec
	return(param_new_vec)
}


options(warn=-1)
system.time(res <- squarem(p=param0,y=y, fixptfn=counts_final.em, control=list(maxiter = 20, trace = FALSE)));

omega_final=matrix(res$par[1:(N*K)],N,K);
alpha_final=matrix(res$par[-(1:(N*K))],K,G)

docweights=omega_final;
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

omega_final=docweights;
alpha_final=alpha_final[perm_set[p_star,],];

windows()
k=K
barplot(t(omega_final),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


windows()
k=K
barplot(t(omega_initial),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)






