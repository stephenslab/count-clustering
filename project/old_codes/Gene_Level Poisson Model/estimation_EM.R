
omega_loglik_Poisson_2 = function(omega_vec,counts_vec,alpha,beta,lab,noise_vec)
{

	G=length(counts_vec);
	sum=0;
	for(g in 1:G)
	{
		index=sample(1:K,1,omega_true[n,],replace=T);
		lambda=exp(alpha[index,g] +beta[Label.Batch[n],g]+noise_vec[g]);

		sum=sum +lambda - counts_vec[g]*log(lambda) ;
		#cat(sum,"\n");
	}

	return(sum);
}


##########################   Estimation Process Iterative Algorithm #############################

omega0=omega_initial;
alpha0=alpha_initial;
noise=noise_true

sigmab_0=0.3;

beta0=matrix(0,B,G);       ###  the matrix of the random effect

for(g in 1:G)
{
	beta0[,g]=rnorm(B,mean=0,sd=sigmab_0);
}

iteration=1; MaxIter=500;

while(iteration<=MaxIter)
{
	################   Estimation of alpha  ###########################

	alpha=matrix(0,K,G);
	for(k in 1:K)
	{
		for(g in 1:G)
		{	
			temp1=0; temp2=0;
			for(n in 1:N)
			{
				temp1=temp1+omega0[n,k]*counts[n,g];
				temp2=temp2+omega0[n,k]*exp(beta0[Label.Batch[n],g]+noise[n,g]);
			}
			alpha[k,g]=log(temp1/temp2);
		}	
	}

	###############   Estimation of omega  ############################

	omega=matrix(0,N,K);
	Z=array(0,c(N,K,G));

	for(n in 1:N)
	{
		for(g in 1:G)
		{
			for(k in 1:K)
			{
				lambda=beta0[Label.Batch[n],g]+alpha[k,g]+noise[n,g];
				Z[n,k,g]=omega0[n,k]*dpois(counts[n,g],exp(lambda));
			}
			Z[n,,g]=(Z[n,,g]+1e-07)/sum(Z[n,,g]+1e-07)
		}
	}

	for(n in 1:N)
	{
		for(k in 1:K)
		{
			omega[n,k]=sum(Z[n,k,])/G;
		}
	}

	windows()
	barplot(t(omega),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


	###############  Estimation of beta  #################################

	beta=matrix(0,B,G);

	for(b in 1:B)
	{
		for(g in 1:G)
		{
			temp1=0; temp2=0;

			for(n in 1:N)
			{
				if(Label.Batch[n]==b)
				{
					temp1=temp1+sum(omega[n,])*counts[n,g];
					temp2=temp2+(omega[n,]%*%exp(alpha[,g]+noise[n,g]));
				}
			}
			beta[b,g]=log(temp1/temp2);
		}
	}

	##############  Estimation of the term epsilon #########################

#	noise=matrix(0,N,G)

#	for(n in 1:N)
#	{
#		for(g in 1:G)
#		{
#			temp1=sum(omega[n,])*(counts[n,g])+1e-07;
#			temp2=(omega[n,]%*%exp(alpha[,g]+beta[Label.Batch[n],g]))+1e-07;
#			noise[n,g]=log(temp1/temp2);
#		}
#	}

	##############  Rearranging the rows of omega and plotting it  ##############

	
	docweights=omega;
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
	barplot(t(docweights),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


	omega0=docweights; alpha0=alpha; beta0=beta; noise=noise_true;
	iteration=iteration+1;
}





		
