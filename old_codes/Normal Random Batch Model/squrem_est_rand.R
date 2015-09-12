

#######  SQUAREM implementaion Normal Random Batch Topic Model ################

omega_loglik_rand = function(omega_vec,counts_vec,alpha,beta,lab)
{

	G=length(counts_vec);
	sum=0;
	for(g in 1:G)
	{
		sum=sum+0.5*(counts_vec[g] - omega_vec%*%alpha[,g] -beta[lab,g])^2;
		#cat(sum,"\n");
	}

	return(sum);
}

counts=read_counts;
y=matrix(counts,1,N*G);
scale=3; K=4;
omega0=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);
alpha0=matrix(rnorm((K)*G,0.5,1),nrow=(K)); ### the matrix of fixed effects


windows()
k=K
barplot(t(omega0),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

param0=c(matrix(omega0,1,N*K),matrix(alpha0,1,K*G));

library(SQUAREM)
library(lme4)

normal_topic.randeff <- function(param_vec,y)
{
	counts=matrix(y,N,G);
	omega0=matrix(param_vec[1:(N*K)],N,K);
	alpha0=matrix(param_vec[-(1:(N*K))],K,G);
	
##############  Estimating the alpha, beta and overdis ####################
	alpha=matrix(0,K,G);
	beta=matrix(0,B,G);
	

	for(g in 1:G)
	{
		counts_col=counts[,g];
		rand_col=as.factor(Label.Batch);
		#noise_col=1:N;
		fit=lmer(counts_col~omega0[,1]+omega0[,2]+omega0[,3]+omega0[,4]+(1|rand_col)-1,
						REML=TRUE);
		#vec_intercept=c(0,rep(as.numeric(fixef(fit))[1],(length(as.numeric(fixef(fit)))-1)));
		alpha[,g]=as.numeric(fixef(fit));
		beta[,g]=as.numeric(as.matrix(ranef(fit)$rand_col));
	}

###################  Topic proportion omega esimate ########################

	omega=matrix(0,N,K);
	for(n in 1:dim(counts)[1])
	{
		omega_vec=omega0[n,];
		counts_vec=counts[n,];
		lab=Label.Batch[n];
		#noise_vec=noise[n,]
		res=optim(reverse_transform(omega_vec), function(v) omega_loglik_rand(transform(v),counts_vec,alpha,beta,lab) );
		omega[n,]=transform(res$par);
	}

	param_vec_omega=matrix(omega,1,N*K);
	param_vec_alpha=matrix(alpha,1,K*G);
	param_new_vec=c(param_vec_omega,param_vec_alpha);
	param_vec=param_new_vec
	return(param_new_vec)
}

options(warn=-1)
system.time(res <- squarem(p=param0,y=y, fixptfn=normal_topic.randeff, control=list(maxiter = 20, trace = FALSE)));

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

	


