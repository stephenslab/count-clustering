

#####  omega_initial  for Poisson Random Batch Model distribution ##########
counts=read_counts;
mean_temp=array(0,c(B,K,G));
omega_preprocess=matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);

 
for(b in 1:B)
{
	
	counts_batch=counts[which(Label.Batch==b),];
	log_counts_batch=log(counts_batch+0.5)
	scale=3; K=4;

	#omega0=matrix(rdirichlet(dim(counts_batch)[1],
	#		c(scale/K,scale/K,scale/K,scale/K)), nrow=dim(counts_batch)[1]);

	# omega0=omega_true[Label.Batch==b,];

	omega0=omega_preprocess[which(Label.Batch==b),];


	iteration=1;
	MaxIter=500
	diff2=1000; diff1=100;

	while(iteration < MaxIter)
	{

		####   Estimation of the matrix H

		svd_omega=svd(omega0);
		temp1=t(svd_omega$v)%*%diag(1/svd_omega$d^2,dim(omega0)[2])%*%svd_omega$v;
		temp2=t(omega0)%*%log_counts_batch;
		temp1=solve(t(omega0)%*%omega0);
		H = temp1%*%temp2;

		###  Estimation of the matrix W (or omega) 

		omega=matrix(0,dim(log_counts_batch)[1],K);
		for(n in 1:dim(log_counts_batch)[1])
		{
			omega_vec=omega0[n,];
			counts_vec=log_counts_batch[n,];
			res=optim(reverse_transform(omega_vec), function(v) loglik_norm(transform(v),counts_vec,t(H)) );
			omega[n,]=transform(res$par);
		}
		diff2=diff1;
		diff1=fnorm(log_counts_batch,omega%*%H);
		cat("The difference is",diff1,"\n");
		omega0=omega;
		iteration=iteration+1;
	}

	mean_temp[b,,]=H;
	omega_preprocess[which(Label.Batch==b),]=omega0;
}

docweights=omega_preprocess;
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

k=K
# Make plot 
# - plot data
windows()
barplot(t(docweights),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

omega_initial=docweights;

	

