

#############  SQUAREM  for Normal Topic Model (No Random effect) #############

counts=read_counts; 
N=dim(counts)[1];  G=dim(counts)[2]; 
y=matrix(counts,1,N*G);

scale=1; K=4;

require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
seed=54321))


omega0=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);
alpha0=matrix(rnorm((K)*G,1,1),nrow=(K)); ### the matrix of fixed effects

windows()
k=K
barplot(t(omega0),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

param0=c(matrix(omega0,1,N*K),matrix(alpha0,1,K*G));

library(SQUAREM)

normal_topic <-function(param_vec,y)
{
	counts=matrix(y,N,G);
	omega0=matrix(param_vec[1:(N*K)],N,K);
	alpha0=matrix(param_vec[-(1:(N*K))],K,G);
	
###########   Estimating the effect size alpha  ########################

	svd_omega=svd(omega0);
	temp1=t(svd_omega$v)%*%diag(1/svd_omega$d^2,dim(omega0)[2])%*%svd_omega$v;
	temp2=t(omega0)%*%counts;
	temp1=solve(t(omega0)%*%omega0);
	alpha = temp1%*%temp2;

###########  Estimating the topic proportions ########################

	
	omega=matrix(0,dim(counts)[1],K);
	for(n in 1:dim(counts)[1])
	{
		omega_vec=omega0[n,];
		counts_vec=counts[n,];
		res=optim(reverse_transform(omega_vec), function(v) loglik_norm(transform(v),counts_vec,t(alpha)) );
		omega[n,]=transform(res$par);
	}

	param_vec_omega=matrix(omega,1,N*K);
	param_vec_alpha=matrix(alpha,1,K*G);
	param_new_vec=c(param_vec_omega,param_vec_alpha);
	param_vec=param_new_vec
	return(param_new_vec)
}

options(warn=-1)
system.time(res <- squarem(p=param0,y=y, fixptfn=normal_topic, control=list(maxiter = 100, trace = FALSE)));

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
par(mar=c(8,5.25,2.5,2.5))

# - get rid of space space between leftmost bar and y axis
par(xaxs="i")

k=K
# Make plot 
# - plot data
barplot(t(omega_final),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

	


	
