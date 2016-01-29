

##########  Estimation of the Poisson Topic Model Version 2 ################

omega_loglik_Poisson_2 = function(omega_vec,counts_vec,alpha,beta,lab,noise_vec)
{

	G=length(counts_vec);
	sum=0;
	for(g in 1:G)
	{
		lambda=omega_vec%*%exp(alpha[,g] +beta[lab,g]+noise_vec[g]);
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
sigmab_0=1;
alpha0=matrix(rnorm((K)*G,1,5),nrow=(K)); ### initial matrix of fixed effects
alpha0=alpha_true+matrix(rnorm(K*G,0,0.9),nrow=K);

library(nlme)
library(lme4)
beta=matrix(0,B,G);
alpha=matrix(0,K,G);
noise=matrix(0,N,G);


##############  True log likelihood calculation ##########################

true_loglik=0;
for(n in 1:dim(counts)[1])
{
	omega_vec=omega_true[n,];
	counts_vec=counts[n,];
	lab=Label.Batch[n];
	noise_vec=noise_true[n,];
	res=optim(reverse_transform(omega_vec), function(v) omega_loglik_Poisson_2(transform(v),counts_vec,alpha_true,beta_true,lab,noise_vec) );
	true_loglik=true_loglik+res$value+(scale/K)*sum(log(omega_true[n,]));
}

###############  Estimation of the Parameters  ############################

iteration=1;
MaxIter=500;
sigmab_0=array(sigmab_true,G);

while(iteration <=MaxIter)
{
	beta=matrix(0,B,G);
	sigmab=array(0,G);
	predict.val=matrix(0,N,G);

	for(g in 1:G)
	{
		counts_col=counts[,g];
		rand_col=as.factor(Label.Batch);
		fixed_eff=log(omega0%*%exp(alpha0[,g]));
		seq=as.factor(1:N);
		
		fit=glmer(counts_col~(1|rand_col)+(1|seq)
				+offset(fixed_eff)-1,family=poisson(),
					control=glmerControl(optimizer="Nelder_Mead",
						optCtrl=list(maxfun=200000)));
		#fit=lmer(log(counts_col+1e-07)~(1|rand_col)+(1|seq)+offset(fixed_eff)-1);
				
		beta[,g]=as.numeric(as.matrix(ranef(fit)$rand_col));
		noise[,g]=as.numeric(as.matrix(ranef(fit)$seq));
		predict.val[,g]=predict(fit,re.form=NULL);
		#sigmab[g]=sqrt(as.numeric(VarCorr(fit)));
	}

	

	omega=matrix(0,N,K);
	curr_loglik=0;
	for(n in 1:dim(counts)[1])
	{
		omega_vec=omega0[n,];
		counts_vec=counts[n,];
		lab=Label.Batch[n];
		noise_vec=noise[n,];
		### res=optim(reverse_transform(omega_vec), function(v) omega_marginal_loglik(transform(v),counts_vec,alpha_true,sigmab_true,lab,noise_vec) );
		res=optim(reverse_transform(omega_vec), function(v) 
			omega_loglik_Poisson_2(transform(v),counts_vec,
							alpha0,beta,lab,noise_vec) );

		omega[n,]=transform(res$par);
		curr_loglik=curr_loglik+res$value+(scale/K)*sum(log(omega[n,]));
	}
	
	residuals.fitted=predict.val-beta[Label.Batch,];
	alpha=matrix(0,K,G);

	for(g in 1:G)
	{
		beta_col=beta[,g];
		counts_col=counts[,g];
		noise_col=noise[,g];
		#res=optim(alpha0[,g],function(v) alpha_normal_distance(v,residuals.fitted[,g],omega));
		res=optim(alpha0[,g],function(v) 
				alpha_likelihood(v,counts_col,omega,beta_col,noise_col),
					method="L-BFGS-B");
		alpha[,g]=res$par;
	}

	omega0=omega; iteration=iteration+1;
	alpha0=alpha; 
	beta0=beta;sigmab_0=sigmab;
	diff=(curr_loglik-true_loglik)/true_loglik;
	cat("The relative diff is",diff,"\n");
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
	
	#docweights=docweights[,c(1,2,3,4)]


	
	windows()
	par(mar=c(8,5.25,2.5,2.5))

	# - get rid of space space between leftmost bar and y axis
	par(xaxs="i")

	barplot(t(docweights),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


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
	
	docweights=docweights[,c(1,3,2,4)]


	
	windows()
	par(mar=c(8,5.25,2.5,2.5))

	# - get rid of space space between leftmost bar and y axis
	par(xaxs="i")

	barplot(t(docweights),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

	windows()
	barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)





