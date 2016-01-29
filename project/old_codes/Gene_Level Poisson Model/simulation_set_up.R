
##################   Simulation set up (Gene Level Poisson) #################

K=4;
G=100;
N=500;

alpha_true=matrix(rnorm((K)*G,1,1),nrow=(K)); ### the matrix of fixed effects

Label.Batch=c(rep(1,N/2),rep(2,N/2)); ##  the batch labels

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

###  generating the table 



over_dis=0.3;

noise_true=matrix(0,N,G);

for(n in 1:N)
{
	noise_true[n,]=over_dis*rnorm(G,0,1);
}


read_counts=matrix(0,N,G);
indicator=matrix(0,N,G);

for(n in 1:N)
{
	for(g in 1:G)
	{
		index=sample(1:K,1,omega_true[n,],replace=T);
		mean=exp(alpha_true[index,g] +beta_true[Label.Batch[n],g]+noise_true[n,g]);
		read_counts[n,g]=rpois(1,mean);
		indicator[n,g]=index;
	}
}

k=K
windows()
barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


i#########  Topic model fit   ####################

############  usual topic model  fit  ##########################

	library(maptpx)

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

	windows()
	par(mar=c(8,5.25,2.5,2.5))

	# - get rid of space space between leftmost bar and y axis
	par(xaxs="i")
	k=K;

	
	windows()
	barplot(t(docweights_topics),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


	




