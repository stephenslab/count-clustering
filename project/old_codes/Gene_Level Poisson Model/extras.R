
#########################  Estimation of the alpha  ##########################

	alpha=matrix(0,K,G);

	for(g in 1:G)
	{
		for(k in 1:K)
		{
			temp1=0; temp2=0;
			for(n in 1:N)
			{
				temp1=temp1+counts[n,g]*omega0[n,k];
				temp2=temp2+omega0[n,k]*exp(beta0[Label.Batch[n],g]);
			}
			alpha[k,g]=log(temp1/temp2);
		}
	}



##########################  Estimation of the omega ##########################

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




############ Determining the initial values  of alpha, beta ##################

counts=read_counts  ###  we start with the read counts table

alpha_initial=matrix(0,K,G);

for(g in 1:G)
{
	for(b in 1:B)
	{
		counts_vec=counts[which(Label.Batch==b),g];
		log_counts_vec=log(counts_vec+0.5);
		#hist(log_counts_vec,nclass=30)
		kmeans_func=kmeans(log_counts_vec,centers=4)
		means=kmeans_func$centers
		alpha_initial[,g]=means-beta_true[b,g];
	}
}





library(HTSCluster)

emInit(matrix(counts[,g]),4,conds=1, lib.size = TRUE, lib.type = "TC", alg.type = "EM")


simulate <- PoisMixSim(n = 500, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions

PoisMixMin(counts[,g],4,conds=1)


library(cluster)
counts1=counts[Label.Batch==1,]
cla=clara(log(counts1+1),4)
centers=cla$medoids

mean_true=matrix(0,K,G);

for(g in 1:G)
{
	mean_true[,g]=alpha_true[,g]+beta_true[1,g];
}


factor_analysis=factanal(counts1,K)




	

