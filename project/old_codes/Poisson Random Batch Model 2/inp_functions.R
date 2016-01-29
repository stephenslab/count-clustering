

alpha_normal_distance <- function(alpha_vec,res_counts_col,omega)
{
	N=dim(counts)[1];
	G=dim(counts)[2];
	K=dim(omega)[2];
	temp=0
	for(n in 1:N)
	{
		temp=temp+(res_counts_col[n] - log(omega[n,]%*%exp(alpha_vec)))^2;
	}
	return (temp)
}


alpha_likelihood <- function(alpha_vec,counts_col,omega,beta_col,noise_col)
{
	N=dim(counts)[1];
	G=dim(counts)[2];
	K=dim(omega)[2];
	temp=0
	for(n in 1:N)
	{
		temp=temp-dpois(counts_col[n],omega[n,]%*%exp(alpha_vec+beta_col[Label.Batch[n]]+noise_col[n]));
	}
	return(temp)
}


#########################################################################

	newcounts_topics=array(0,c(N,K,G));
	for(n in 1:N)
	{
		for(g in 1:G)
		{
			for(k in 1:K)
			{
				newcounts_topics[n,k,g]=rpois(1,exp(alpha0[k,g]+beta0[Label.Batch[n],g]));
	
			}
		}
	}

###  this is the code for generating pseudocounts from alpha0 and beta0- from

### one mixture component

############################################################################

mean_est=array(0,c(B,K,G));
	mean_est_2=array(0,c(B,K,G));

	alpha=matrix(0,K,G);



	for(g in 1:G)
	{
		sigmab_temp=mean(sigmab_0);
		for(k in 1:K)
		{
			temp1=0; temp2=0;
			for(b in 1:B)
			{
				
				for(n in 1:N)
				{
					if(Label.Batch[n]==b)
					{
						temp1=temp1+counts[n,g]*omega0[n,k];
						temp2=temp2+omega0[n,k]*exp(0.5*sigmab_temp^2);
					}
				}
				mean_est_2[b,k,g]=exp(alpha0[k,g]+beta0[b,g])
			}
			alpha[k,g]=log(temp1/temp2);
		}
	}

###  estimation of alpha using the marginal likelihood method (integrating over 

###  beta) a lower bound to the actual log likelihood 


############################################################################



#	pseudo_counts=array(0,c(N,K,G));
#
#	for(n in 1:N)
#	{
#		for(k in 1:K)
#		{
#			for(g in 1:G)
#			{
#				pseudo_counts[n,k,g]=rpois(1,mean_est[Label.Batch[n],k,g]);
#			}
#		}
#	}

##  Another pseduo counts generation algorithm 


	#b=1;
	#newcounts_topics[which(Label.Batch==b),k,g]

#	beta=matrix(0,B,G);
#     alpha=matrix(0,K,G);

#	for(g in 1:G)
#	{
#		pseudo_count_response=matrix(pseudo_counts[,,g],nrow=N*K);
#		rand_eff=as.factor(rep(Label.Batch,K));
#		fixed_eff=as.factor(rep(1:K,each=N));
#		fit=glmer(pseudo_count_response~fixed_eff+(1|rand_eff),family=poisson());
#		vec_intercept=c(0,rep(as.numeric(fixef(fit))[1],(length(as.numeric(fixef(fit)))-1)));
#		alpha[,g]=as.numeric(fixef(fit))+vec_intercept;
#		beta[,g]=as.numeric(as.matrix(ranef(fit)$rand_eff));
#	}

##  estimation of alpha and beta based on pseudo counts technique  ######


###########################################################################

	alpha=matrix(0,K,G);
	for(g in 1:G)
	{
		for(k in 1:K)
		{
			temp1=0; temp2=0;
			for(n in 1:N)
			{
				temp1=temp1+counts[n,g]*omega0[n,k];
				temp2=temp2+omega0[n,k]*exp(beta[Label.Batch[n],g]);
					
			}
			alpha[k,g]=log(temp1/temp2);
		}
	}

##  estimaing alpha when omega and beta are known  ######################

###################### Using the Z indicator variables  ####################

	alpha=matrix(0,K,G);

	for(g in 1:G)
	{
		Z=array(0,N);
		for(n in 1:N)
		{
			Z[n]=sample(1:K,1,omega[n,],replace=T);
		}
		
		fixed_eff=factor(Z);
		counts_col=counts[,g];
		rand_eff=beta[Label.Batch,g];
		fit=glm(counts_col~fixed_eff+offset(rand_eff)-1,
					family=poisson());
		alpha[,g]=as.numeric(coef(fit));
	}


