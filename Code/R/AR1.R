
ar1.llh.C <- function(beta, mu, phi, W, m0, C0, alpha=NULL)
{
  T   = ncol(beta) - 1;
  N.b = length(mu);
  N   = length(m0);
  N.a = N - N.b

  ## Check
  not.ok = rep(0, 8);
  if (not.ok[1] <- length(mu)  != N.b)
    { print("length(mu)!=N.b") ; }
  if (not.ok[2] <- length(phi) != N.b)
    { print("length(phi)!=N.b"); }
  if (not.ok[3] <- (ncol(W) != N.b || nrow(W) != N.b))
    { print("W is not N.b x N.b"); }
  if (not.ok[4] <- length(m0) != N)
    { print("length(m0) != N"); }
  if (not.ok[5] <- (nrow(C0) != N || ncol(C0) != N))
    { print("C0 is not N x N"); }
  if (not.ok[6] <- N.b > N)
    { print("N.b > N"); }
  if (not.ok[7] <- (nrow(beta) != N.b && ncol(beta) != (T+1)))
    { print("beta is not N.b x T"); }

  if (!is.null(alpha)) {
    if (not.ok[8] <- (length(alpha) != N.a))
      { print("length(alpha) != N.a"); }
  }
  else {
    alpha = 0;
  }
  
  if (!prod(!not.ok)) {
    cat("ar1.llh.C: problem.  Returning NA.\n");
    return(NA)
  }    

  log.dens  = 0
  
  OUT <- .C("ar1_llh", alpha, beta,
            mu, phi, W,
            m0, C0,
            as.integer(N.b), as.integer(N), as.integer(T),
            log.dens,
            PACKAGE="BayesLogit");
     
  OUT[[11]]

}
