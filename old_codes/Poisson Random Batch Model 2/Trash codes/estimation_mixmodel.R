
##########  Estimation of the Poisson Topic Model Version 2 ################

omega_loglik_Poisson_2 = function(omega_vec,counts_vec,alpha,beta,lab)
{

	G=length(counts_vec);
	sum=0;
	for(g in 1:G)
	{
		index=sample(1:K,1,omega_vec,replace=T);
		lambda=exp(alpha[index,g] +beta[lab,g]);
		# cat(lambda,g,"\n")
		sum=sum +lambda -counts_vec[g]*log(lambda) ;
		# cat(sum,"\n");
	}

	return(sum);
}

omega_marginal_loglik =  function(omega_vec, counts_vec, alpha,sigma_b,lab)
{
	G=length(counts_vec);
	sum=0;
	for(g in 1:G)
	{
		sum=sum+(omega_vec%*%exp(alpha[,g]))*exp(0.5*sigma_b^2)
					-(omega_vec%*%alpha[,g])*counts_vec[g]
					+sum(omega_vec)*sum(log(1:max(counts[n,g],1)));
	}
	return(sum)
}

		

###############  Counts set up Parameter Initialization ##################


counts=read_counts;
scale=3; K=4;
Label.Batch=c(rep(1,N/2),rep(2,N/2));

omega0=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);
omega0=matrix(rep(rep(1/K,K),N),nrow=N);
sigmab_0=1;
alpha0=matrix(rnorm((K)*G,1,5),nrow=(K)); ### initial matrix of fixed effects
alpha0=alpha_true+matrix(rnorm(K*G,0,1),nrow=K);

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
sigmab_0=array(sigmab_true,G);

while(iteration <=MaxIter)
{
	
#	for(g in 1:G)
#	{
#		counts_col=counts[,g];
#		rand_col=as.factor(Label.Batch);
#		fixed_eff=log(omega0%*%exp(alpha0[,g]));
#		
#		fit=glmer(counts_col~(1|rand_col)+offset(fixed_eff)-1,
#			    family=poisson(),control=glmerControl(
#					optimizer="bobyqa",optCtrl=list(maxfun=100000)));
#		beta[,g]=as.numeric(as.matrix(ranef(fit)$rand_col));
#		predict.val[,g]=predict(fit,re.form=NULL);
#		sigmab[g]=sqrt(as.numeric(VarCorr(fit)));
#	}
#
#

	mean_intensity=array(0,c(B,K,G));

	for(b in 1:B)
	{
		for(k in 1:K)
		{
			for(g in 1:G)
			{
				mean_intensity[b,k,g]=exp(alpha0[k,g]+beta0[b,g]);
			}
		}
	}

	omega_new=array(0,c(N,K,G))

	for(n in 1:N)
	{
		for(k in 1:K)
		{
			
			for(g in 1:G)
			{
				omega_new[n,k,g]=omega0[n,k]*dpois(counts[n,g],mean_intensity[Label.Batch[n],k,g]);
			}
			
		}
	}

	indicator_new=matrix(0,N,G);
	

	for(n in 1:N)
	{
		for(g in 1:G)
		{
			class_prob_vec=(omega_new[n,,g]+1e-07)/sum(omega_new[n,,g]+1e-07);
			indicator_new[n,g]=sample(1:K,1,class_prob_vec,replace=T);
		}
	}

	prop=length(which((indicator-indicator_new)!=0))/(N*G);

	alpha=matrix(0,K,G);
	beta=matrix(0,B,G);
	predict.val=matrix(0,N,G);

	for(g in 1:G)
	{
		Z=array(0,N);
		for(n in 1:N)
		{
			Z[n]=sample(1:K,1,omega_true[n,],replace=T);
		}
		
		fixed_eff=factor(indicator_new[,g]);
		counts_col=counts[,g];
		rand_eff=as.factor(Label.Batch);
		fit=glmer(counts_col~fixed_eff+(1|rand_eff)-1,
					family=poisson());
		alpha[,g]=as.numeric(fixef(fit));
		beta[,g]=as.numeric(as.matrix(ranef(fit)$rand_eff));
		
	}

	omega=matrix(0,N,K);
	curr_loglik=0;
	for(n in 1:dim(counts)[1])
	{
		omega_vec=omega0[n,];
		counts_vec=counts[n,];
		lab=Label.Batch[n];
		### res=optim(reverse_transform(omega_vec), function(v) omega_marginal_loglik(transform(v),counts_vec,alpha_true,sigmab_true,lab) );
		res=optim(reverse_transform(omega_vec), 
			function(v) omega_loglik_Poisson_2(transform(v),
					counts_vec,alpha,beta,lab) );

		omega[n,]=transform(res$par);
		curr_loglik=curr_loglik+res$value+(scale/K)*sum(log(omega[n,]));
	}

	omega0=omega; iteration=iteration+1;
	alpha0=alpha; 
	beta0=beta;
	diff=(curr_loglik-true_loglik)/true_loglik;
	cat("The relative diff is",diff,"\n");
	cat("proportion of misclassification",prop,"\n");
	windows()
	barplot(t(omega0),col=2:(k+1),axisnames=F,space=0,
			border=NA,main=paste("No. of clusters=",k),
				las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

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

	barplot(t(docweights),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

	windows()
	barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)



