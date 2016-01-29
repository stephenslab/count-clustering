
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

############################################################################

#######################  Initialization  #################################
counts=read_counts;
scale=3; K=4;
omega0=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);
library(nlme)
library(lme4)
beta=matrix(0,B,G);
alpha=matrix(0,K,G);


###########################  True log likelihood values ####################

true_loglik=0;
for(n in 1:dim(counts)[1])
{
	omega_vec=omega_true[n,];
	counts_vec=counts[n,];
	lab=Label.Batch[n];
	res=optim(reverse_transform(omega_vec), function(v) omega_loglik(transform(v),counts_vec,alpha_true,beta_true,lab) );
	true_loglik=true_loglik+res$value;
}

#######################  Estimation of parameters #######################

iteration=1;
MaxIter=500;


while(iteration <= MaxIter)
{
	for(g in 1:G)
	{
		counts_col=counts[,g];
		rand_col=as.factor(Label.Batch);
		
		fit=lmer(counts_col~omega0[,1]+omega0[,2]+omega0[,3]+omega0[,4]+(1|rand_col)-1,
						REML=TRUE);
		alpha[,g]=as.numeric(fixef(fit));
		beta[,g]=as.numeric(as.matrix(ranef(fit)$rand_col));
	}

	omega=matrix(0,N,K);
	#curr_loglik=0;
	for(n in 1:dim(counts)[1])
	{
		omega_vec=omega0[n,];
		counts_vec=counts[n,];
		lab=Label.Batch[n];
		res=optim(reverse_transform(omega_vec), function(v) omega_loglik_rand(transform(v),counts_vec,alpha,beta,lab) );
		omega[n,]=transform(res$par);
		#curr_loglik=curr_loglik+res$value;
	}

	omega0=omega; iteration=iteration+1;
	#diff=(curr_loglik-true_loglik)/true_loglik;
	#cat("The relative diff is",diff,"\n");

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
}


windows()
barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)



#############################  The  End  ###################################


	

		