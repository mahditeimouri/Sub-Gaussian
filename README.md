# Sub-Gaussian
$\colr{blue}{An R code computing the observed Fisher information (OFI) matrix for the sub-Gaussian elliptically alpha-stable distribution.
The arguments for the main function OFI(X, alpha, Sigma, Mu) are X, alpha, Sigma, and Mu defined as follows. }$

$\color{red}{X}$     := an n by d matrix of n realizations each following a d-dimensional sub-Gaussian elliptically alpha-stable distribution,

alpha := the ML estimator for alpha, that is, the tail index parameter,

Sigma := the ML estimator for dispersion matrix,

Mu    := the ML estimator for location vector.

stochastic = FALSE := A logical statement for computing the partial derivative of log-likelihood function w.r.t tail index. By default it is FALSE indicating that Gaussian quadrature (GQ) is used for this mean and otherwise the MCMC method is used. Note that, in general, GQ is faster than the MCMC method. 

The output of function OFI(X, alpha, Sigma, Mu) has two parts consisting of: i- the asymptotic standard error for the ML estimators of the tail index, lower triangular elements of the dispersion matrix, and location parameters; and ii- the OFI matrix.
