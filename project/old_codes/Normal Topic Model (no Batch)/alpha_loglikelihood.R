

#########  This additional script is focused on estimating alpha #############

# meth_loglik_EM_alpha = function(alpha_val,meth,omega,phi,gene_label,k_index)
# {
#	mu=pnorm(alpha_val);
#	first_term=digamma(mu*phi[gene_label])*phi[gene_label]*sum(omega[,k_index]);
#	sec_term=(digamma((1-mu)*phi[gene_label])*phi[gene_label])*sum(omega[,k_index]);
#	third_term =phi[gene_label]*(t(omega[,k_index])%*%log(meth[,gene_label]/(1-meth[,gene_label])));
#	out=(first_term-sec_term+third_term)*dnorm(alpha_val);
#	return(out)
# }


meth_loglik_EM_alpha = function(alpha_vec,meth,omega,phi,gene_label)
{
	g=gene_label;
	mu=pnorm(alpha_vec);
	sum=0;
	for(n in 1:N)
	{
		out=0
		for(k_index in 1:K)
		{
			out=out+omega[n,k_index]*dbeta(meth[n,g],
					shape1=mu[k_index]*phi[g],
						shape2=(1-mu[k_index])*phi[g],ncp=0);
		}
		sum=sum-log(out);
	}
	return(abs(sum))
}

optim(alpha_true[,gene_label], function(v) meth_loglik_EM_alpha(v,meth,omega_true,phi_true,gene_label),
		method="Nelder-Mead")

meth_loglik_EM_alpha(alpha_true[,gene_label],meth,omega_true,phi_true,gene_label)




