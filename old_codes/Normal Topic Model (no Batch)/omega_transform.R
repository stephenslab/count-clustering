
	a <- omega_vec
	reverse_transform=function(x) 
	{
		out=array(0,K-1);
		for(i in 1:(K-1))
		{
			out[i]=log((x[i]+10^(-5))/(x[K]+10^(-5)));
		}
		return(out)
	}


	# Data
	

	# Log-Likelihood
	loglik_norm <- function(u,y,x) sum((y - x %*% u)^2)

	# Transform the parameters: we just have
	# to find a bijection between R^3 
	# and {(a,b,c,d) \in [0,1]^4 : a+b+c+d=1}.

	transform <- function(v) 
	{
  	# Ensure they in [0,1]
  	temp =c(exp(v),1);
  	out=temp/(1+sum(exp(v)));
  	return(out)
	}
# Minimize the log-likelihood
r <- optim(reverse_transform(omega_vec), function(v) loglik_norm(transform(v),counts_vec,t(H)), method = "L-BFGS-B" )
val=transform(r$par); 