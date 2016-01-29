
omega_loglik_Poisson = function(omega_vec,counts_vec,alpha,beta,lab,noise_vec)
{

	G=length(counts_vec);
	sum=0;
	for(g in 1:G)
	{
		lambda=exp(omega_vec%*%alpha[,g] +beta[lab,g]+noise_vec[g]);
		#cat(lambda,"\n")
		sum=sum +lambda - counts_vec[g]*log(lambda) ;
		#cat(sum,"\n");
	}

	return(sum);
}

#####################  Counts set up initials ##########################

counts=read_counts;
scale=3; K=4;
# omega0=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);

omega0=omega_initial
library(nlme)
library(lme4)
beta=matrix(0,B,G);
alpha=matrix(0,K,G);
noise=matrix(0,N,G);


###########################  True log likelihood values ####################

true_loglik=0;
for(n in 1:dim(counts)[1])
{
	omega_vec=omega_true[n,];
	counts_vec=counts[n,];
	lab=Label.Batch[n];
	noise_vec=noise_true[n,];
	res=optim(reverse_transform(omega_vec), function(v) omega_loglik_Poisson(transform(v),counts_vec,alpha_true,beta_true,lab,noise_vec) );
	true_loglik=true_loglik+res$value+(scale/K)*sum(log(omega_true[n,]));
}


############################  Estimation  ###############################

iteration=1;
MaxIter=500;


while(iteration <= MaxIter)
{
	for(g in 1:G)
	{
		counts_col=counts[,g];
		rand_col=as.factor(Label.Batch);
		seq=as.factor(1:N);
		
		fit=glmer(counts_col~omega0[,2]+omega0[,3]+omega0[,4]+(1|rand_col)
					+(1|seq),family=poisson());
		vec_intercept=c(0,rep(as.numeric(fixef(fit))[1],(length(as.numeric(fixef(fit)))-1)));
		alpha[,g]=as.numeric(fixef(fit))+vec_intercept;
		beta[,g]=as.numeric(as.matrix(ranef(fit)$rand_col));
		noise[,g]=as.numeric(as.matrix(ranef(fit)$seq));
	}

	omega=matrix(0,N,K);
	curr_loglik=0;
	for(n in 1:dim(counts)[1])
	{
		omega_vec=omega0[n,];
		counts_vec=counts[n,];
		lab=Label.Batch[n];
		noise_vec=noise[n,];
		res=optim(reverse_transform(omega_vec), function(v) omega_loglik_Poisson(transform(v),counts_vec,alpha,beta,lab,noise_vec) );
		omega[n,]=transform(res$par);
		curr_loglik=curr_loglik+res$value+(scale/K)*sum(log(omega[n,]));
	}
	
	omega0=omega; iteration=iteration+1;
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
	

	windows()
	par(mar=c(8,5.25,2.5,2.5))

	# - get rid of space space between leftmost bar and y axis
	par(xaxs="i")

	k=K
	# Make plot 
	# - plot data
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

	

