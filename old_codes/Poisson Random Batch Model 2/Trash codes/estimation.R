
##########  Estimation of the Poisson Topic Model Version 2 ################

omega_loglik_Poisson_2 = function(omega_vec,counts_vec,alpha,beta,lab)
{

	G=length(counts_vec);
	sum=0;
	for(g in 1:G)
	{
		lambda=omega_vec%*%exp(alpha[,g] +beta[lab,g]);
		# cat(lambda,g,"\n")
		sum=sum +lambda - counts_vec[g]*log(lambda) ;
		# cat(sum,"\n");
	}

	return(sum);
}

###############  Counts set up Parameter Initialization ##################


counts=read_counts;
scale=3; K=4;
Label.Batch=c(rep(1,N/5),rep(2,N/5),rep(3,N/5),rep(4,N/5),rep(5,N/5));

omega0=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);
sigmab_0=2;

beta0=matrix(0,B,G);       ###  the matrix of the random effect

for(g in 1:G)
{
	beta0[,g]=rnorm(B,mean=0,sd=sigmab_0);
}

alpha0=matrix(rnorm((K)*G,1.5),nrow=(K)); ### the matrix of fixed effects


library(nlme)
library(lme4)
beta=matrix(0,B,G);
alpha=matrix(0,K,G);


##############  True log likelihood calculation ##########################

true_loglik=0;
for(n in 1:dim(counts)[1])
{
	omega_vec=omega_true[n,];
	counts_vec=counts[n,];
	lab=Label.Batch[n];
	res=optim(reverse_transform(omega_vec), function(v) omega_loglik_Poisson_2(transform(v),counts_vec,alpha_true,beta_true,lab) );
	true_loglik=true_loglik+res$value+(scale/K)*sum(log(omega_true[n,]));
}

###############  Estimation of the Parameters  ############################

iteration=1;
MaxIter=500;

while(iteration <=MaxIter)
{
	mean_intensity=array(0,c(N,K,G));

	for(n in 1:N)
	{
		for(g in 1:G)
		{
			for(k in 1:K)
			{
				mean_intensity[n,k,g]=omega0[n,k]*dpois(counts[n,g],exp(alpha0[k,g] +beta0[Label.Batch[n],g]));
			}
		}
	}

	c_bkg_hat=array(0,c(B,K,G));



	for(g in 1:G)
	{
		for(k in 1:K)
		{
			for(b in 1:B)
			{
				temp=0; 
				for(n in 1:N)
				{
					if(Label.Batch[n]==b)
					{
						temp=temp+counts[n,g]*(mean_intensity[n,k,g]+10^(-5))/sum(mean_intensity[n,,g]+10^(-5));
					}
				}
				c_bkg_hat[b,k,g]=temp;
			}
		}
	}

	for(g in 1:G)	
	{
		response=matrix(c_bkg_hat[,,g],nrow=B*K);
		fixed_eff = as.factor(rep(1:K,each=B));
		rand_eff = as.factor(rep(1:B,K));
		fit=glmer(response~fixed_eff+(1|rand_eff),family=poisson(),REML=TRUE);
		vec_intercept=c(0,rep(as.numeric(fixef(fit))[1],(length(as.numeric(fixef(fit)))-1)));
		alpha[,g]=as.numeric(fixef(fit))+vec_intercept;
		beta[,g]=as.numeric(as.matrix(ranef(fit)$rand_eff));
	}

	omega=matrix(0,N,K);
	curr_loglik=0;
	for(n in 1:dim(counts)[1])
	{
		omega_vec=omega0[n,];
		counts_vec=counts[n,];
		lab=Label.Batch[n];
		res=optim(reverse_transform(omega_vec), function(v) omega_loglik_Poisson_2(transform(v),counts_vec,alpha,beta,lab) );
		omega[n,]=transform(res$par);
		curr_loglik=curr_loglik+res$value+(scale/K)*sum(log(omega[n,]));
	}

	omega0=omega; iteration=iteration+1;
	diff=(curr_loglik-true_loglik)/true_loglik;
	cat("The relative diff is",diff,"\n");
}



docweights=omega0;
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


windows()
par(mar=c(8,5.25,2.5,2.5))

# - get rid of space space between leftmost bar and y axis
par(xaxs="i")

k=K
# Make plot 
# - plot data
barplot(t(docweights),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

windows()
barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


##################################################


