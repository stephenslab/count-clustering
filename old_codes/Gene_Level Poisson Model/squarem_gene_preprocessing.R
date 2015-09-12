counts=read_counts;

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

counts_batch_process.em <- function(param_vec_batch,y_batch)
{
	Nb=length(y_batch)/G;
	log_counts_batch=matrix(y_batch,Nb,G);
	omega0_batch=matrix(param_vec_batch[1:(Nb*K)],Nb,K);
	alpha0_batch=matrix(param_vec_batch[-(1:(Nb*K))],K,G);

	temp2=t(omega0_batch)%*%log_counts_batch;
	# temp2=t(omega0_batch)%*%pnorm_meth;
	temp1=solve(t(omega0_batch)%*%omega0_batch);
	H = temp1%*%temp2;

###  Estimation of the matrix W (or omega) 

	omega_batch=matrix(0,Nb,K);
	for(n in 1:Nb)
	{
		omega_vec=omega0_batch[n,];
		counts_vec=log_counts_batch[n,];
		res=optim(reverse_transform(omega_vec), function(v) loglik_norm(transform(v),counts_vec,t(H)) );
		omega_batch[n,]=transform(res$par);
	}
	param_vec_omega=matrix(omega_batch,1,Nb*K);
	param_vec_alpha=matrix(H,1,K*G);
	param_new_vec=c(param_vec_omega,param_vec_alpha);
	param_vec=param_new_vec
	return(param_new_vec)
}



omega0=omega_preprocess
alpha0=matrix(rnorm((K)*G,1,1),nrow=(K));
param0=c(matrix(omega0,1,N*K),matrix(alpha0,1,K*G));



counts_process.em <- function(param_vec,y)
{
	log_counts=matrix(y,N,G);
	omega0=matrix(param_vec[1:(N*K)],N,K);
	alpha0=matrix(param_vec[-(1:(N*K))],K,G);
	mean_effect=array(0,c(B,K,G));

	omega=matrix(0,N,K);
	alpha=matrix(0,K,G);

	for(b in 1:B)
	{
		log_counts_batch=log_counts[which(Label.Batch==b),];
		omega0_batch=omega0[which(Label.Batch==b),];
		Nb=length(which(Label.Batch==b));
		y_batch=matrix(log_counts_batch,1,Nb*G)
		alpha0_batch=matrix(rnorm((K)*G,1,1),nrow=(K));
		param0_batch=c(matrix(omega0_batch,1,Nb*K),matrix(alpha0_batch,1,K*G));

		options(warn=-1)
		system.time(res <- squarem(p=param0_batch,y_batch=y_batch, 
		fixptfn=counts_batch_process.em, control=list(maxiter = 20, trace = FALSE)));

		omega_batch=matrix(res$par[1:(Nb*K)],Nb,K);
		alpha_batch=matrix(res$par[-(1:(Nb*K))],K,G);
		
		docweights=omega_batch;
		library(permute);
		library("BioPhysConnectoR");
		perm_set=rbind(1:K,allPerms(1:K));
		diff=array(0,dim(perm_set)[1]);
		for (p in 1:dim(perm_set)[1])
		{
			temp=docweights[,perm_set[p,]];
			diff[p]=fnorm(temp,omega_true[which(Label.Batch==b),]);
		}

		p_star=which(diff==min(diff));
		docweights=docweights[,perm_set[p_star,]];

		omega_batch=docweights;

		alpha_batch=alpha_batch[perm_set[p_star,],];

		omega[which(Label.Batch==b),]=omega_batch;
		mean_effect[b,,]=alpha_batch;
	}

	for(k in 1:K)
	{
		for(g in 1:G)
		{
			alpha[k,g]=mean(mean_effect[,k,g]);
		}
	}

	param_vec_omega=matrix(omega,1,N*K);
	param_vec_alpha=matrix(alpha,1,K*G);
	param_new_vec=c(param_vec_omega,param_vec_alpha);
	param_vec=param_new_vec
	return(param_new_vec)
}

options(warn=-1)
system.time(res <- squarem(p=param0,y=y, fixptfn=counts_process.em, 
				control=list(maxiter = 100, trace = FALSE)));


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

######################  Using the alpha log likelihood to modify alpha est  #########################

alpha0=alpha_initial;
omega0=omega_initial; 

	
omega_preprocess=omega_initial;
alpha=alpha_initial;
beta=beta_true
noise=noise_true



	





