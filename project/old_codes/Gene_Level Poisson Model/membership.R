

#############  Gene Membership Probabilities (EM model )  #################

counts=read_counts;


###  Take alpha, beta and omega as finally obtained from the EM algorithm

membership=matrix(0,K,G);

for(g in 1:G)
{
	for(k in 1:K)
	{
		temp1=0; temp2=0;
		for(n in 1:N)
		{

			mean=exp(alpha[k,g]+beta[Label.Batch[n],g]);
			temp1=omega[n,k]*dpois(counts[n,g],mean);
		}
		membership[k,g]=temp1;
	}
	membership[,g]=membership[,g]/(sum(membership[,g]));
}


