# Sub-Gaussian
An R code computing the observed Fisher information (OFI) matrix for the sub-Gaussian elliptically alpha-stable distribution.
The arguments for the main function OFI( $\color{blue}{X}$, $\color{blue}{\alpha}$ $\color{blue}{\Sigma}$, $\color{blue}{\boldsymbol{\mu}}$, $\color{blue}{\text{MCMC = FALSE}}$) are $\color{blue}{X}$ , $\color{blue}{\alpha}$, $\color{blue}{\Sigma}$, $\color{blue}{\boldsymbol{\mu}}$, and ${\text{MCMC = FALSE}}$ defined as follows.

$\color{blue}{X}$     := an n by d matrix of n realizations each following a d-dimensional sub-Gaussian elliptically alpha-stable distribution,

$\color{blue}{\alpha}$ := the ML estimator for alpha, that is, the tail index parameter,

$\color{blue}{\Sigma}$ := the ML estimator for dispersion matrix,

$\color{blue}{\boldsymbol{\mu}}$    := the ML estimator for location vector.

$\color{blue}{\tet{MCMC = FALSE}}$ := A logical statement for computing the partial derivative of log-likelihood function w.r.t tail index. By default it is FALSE indicating that Gaussian quadrature (GQ) is used for this mean and otherwise the MCMC method is used. Note that, in general, GQ is faster than the MCMC method. 

The output of function OFI(X, alpha, Sigma, Mu, MCMC = FALSE) has two parts consisting of: i- the asymptotic standard error for the ML estimators of the tail index, lower triangular elements of the dispersion matrix, and location parameters; and ii- the OFI matrix.
