

K=4;
G=100;
N=500;

alpha_true=matrix(rnorm((K)*G,1,3),nrow=(K)); ### the matrix of fixed effects

alpha_true=matrix(rep(c(-3,-1,1,3),G),ncol=G);
Label.Batch=c(rep(1,N/2),rep(2,N/2));

B=max(Label.Batch);

sigmab_true=1;

beta_true=matrix(0,B,G);       ###  the matrix of the random effect

for(g in 1:G)
{
	beta_true[,g]=rnorm(B,mean=0,sd=sigmab_true);
}

library(gtools)
T=10;
omega_true=matrix(rbind(rdirichlet(T*10,c(3,4,2,6)),rdirichlet(T*10,c(1,4,6,3)),
			rdirichlet(T*10,c(4,1,2,2)),rdirichlet(T*10,c(2,6,3,2)),
			rdirichlet(T*10,c(3,3,5,4))), nrow=N);


H_true=alpha_true;
alpha_true=H_true

###  generating the table 


read_counts=matrix(0,N,G);
indicator=matrix(0,N,G);

for(n in 1:N)
{
	for(g in 1:G)
	{
		index=sample(1:K,1,omega_true[n,],replace=T);
		mean=exp(alpha_true[index,g] +beta_true[Label.Batch[n],g]);
		read_counts[n,g]=rpois(1,mean);
		indicator[n,g]=index;
	}
}


i#########  Topic model fit   ####################

############  usual topic model  fit  ##########################


	Topic_Clus=topics(read_counts,K,kill=0,tol=0.01);
	docweights_topics=Topic_Clus$omega;
	library(permute);
	library("BioPhysConnectoR");
	perm_set=rbind(1:K,allPerms(1:K));
	diff=array(0,dim(perm_set)[1]);
	for (p in 1:dim(perm_set)[1])
	{
		temp=docweights_topics[,perm_set[p,]];
		diff[p]=fnorm(temp,omega_true);
	}

	p_star=which(diff==min(diff));
	docweights_topics=docweights_topics[,perm_set[p_star,]];

