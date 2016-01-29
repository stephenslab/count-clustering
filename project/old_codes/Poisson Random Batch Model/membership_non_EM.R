

############  Gene membership Probabilities (Random Batch Model) ###########

counts=read_counts;


###  Take alpha, beta and omega as finally obtained from the EM algorithm

membership=matrix(0,K,G);

for(g in 1:G)
{
	for(n in 1:N)
	{
		mean_intensity=array(0,c(N,K));
		for(k in 1:K)
		{
			mean_intensity[n,k]=omega[n,k]*exp(alpha[k,g]+beta[Label.Batch[n],g]);
		}
		mean_intensity[n,]=mean_intensity[n,]/sum(mean_intensity[n,]);
	}
	for(k in 1:K)
	{
		membership[k,g]=sum(mean_intensity[,k])/sum(mean_intensity);
	}
}


