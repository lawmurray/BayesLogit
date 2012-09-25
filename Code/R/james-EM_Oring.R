
EM_Oring <- function(x,y,n)
{
# x: N by p matrix
# y: N by 1 vector, avg response
# n: N by 1 vector, # of obs at distinct x

  p = ncol(x)
  N = nrow(x)

  Z = colSums(x*(y-1/2)*n)
  w = rep(0,N)
  beta = rep(0,p)

  while (TRUE)
  {

    # w: posterior mean
    psi = drop(x%*%beta)
    for ( i in 1:N )
      if ( psi[i]<0.01 )
      {
        b = psi[i]/2
        w[i] = n[i] * cosh(psi[i]/2)^(-1) / 4 * (1 + b^2/6 + b^4/120 + b^6/5040)
      }
      else
        w[i] = n[i] * tanh(psi[i]/2) / psi[i] / 2

    # beta: posterior mode
    prebeta = beta
    ch = chol(t(x)%*%(x*w))

    print(eigen(t(x)%*%(x*w))$value[1]);

    beta = backsolve(ch, forwardsolve(t(ch),Z))
    if ( max(abs(beta-prebeta)) < 1e-8 ) break
  }

  beta
}


#data = read.table("orings.dat",header=TRUE)
#attach(data)
#failure = 2*failure-1

## x = c(53,56,57,63,66,67,68,69,70,72,73,75,76,78,79,80,81)
## y = c(1,1,1,0,0,0,0,0,3/4,0,0,1/2,0,0,0,0,0)
## n = c(1,1,1,1,1,3,1,1,4,1,1,2,2,1,1,1,1)
## ans = EM_Oring(cbind(1,x),y,n)



