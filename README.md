# Sub-Gaussian
An R code computing the observed Fisher information (OFI) matrix for the sub-Gaussian elliptically alpha-stable distribution.
The arguments for the main function OFI(X, alpha, Sigma, Mu) are X, alpha, Sigma, and Mu defined as follows. 

X         := an n by d matrix of n realizations each following a d-dimensional sub-Gaussian elliptically alpha-stable distribution,

alpha := the ML estimator for alpha, that is, the tail index parameter,

Sigma := the ML estimator for dispersion matrix,

Mu    := the ML estimator for location vector.

The output of function OFI(X, alpha, Sigma, Mu) has two parts consisting of: i- the asymptotic standard error for the ML estimators of the tail index, lower triangular elements of the dispersion matrix, and location parameters; and ii- the OFI matrix.
