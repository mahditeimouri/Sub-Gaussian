arrange_sigma <- function(x)
{
	Dim <- length( x[1, ] )
	y <- rep( 0, Dim*(Dim + 1)/2 )
	index <- matrix(0, nrow = Dim*(Dim + 1)/2, ncol = 2)
	k <- 1
		for(i in 1:Dim)
		{
			for(j in i:Dim)
			{
				y[k] <- x[i, j]
				index[k, ] <- c(i, j)
				k <- k + 1
			}
		}
return( list(y = y, index = index) )
}

OFI <- function(X, alpha, Sigma, Mu, stochastic = FALSE)
{
	X <- as.matrix(X)
	if ( is.null(X) )       stop("data mst be given in a matrix form.")
	if ( any( is.na(X) ) )  stop("NAs values are not allowed for matrix of observations.")
	Dim_x <- dim(X)
	n <- Dim_x[1]
	d <- Dim_x[2]
	if( length( Mu ) != d ) stop( "Length of location vector must be equal to d." )
	if( length( Sigma[ ,1] ) != d & length( Sigma[1, ] ) != d ) stop( "dispersion matrix must be square." )
	if( any( eigen(Sigma)$values < 0 ) ) stop( "dispersion matrix must be positive definite." )
	Dim <- d*(d + 1)/2 + d + 1
	D  <- mahalanobis(X, Mu, Sigma)
	N  <- 3000
	M  <- 20
	K1 <- 60
	K2 <- c(1:K1)
	T1 <- (K1*alpha + d)*(K1*alpha + 2)*(K1 + 1)^(-2/alpha) + 2
	E1 <- rep(NA, n)
	epsilon <- 10e-9
	index_T1 <- which( D >= T1 )
	if( all(index_T1 == FALSE ) ){ index_C <- c(1:n) }else{ index_C <- c(1:n)[-index_T1] }
	D_T1 <- D[ index_T1 ]
	D_C  <- D[ index_C ]
	n_T1 <- length( index_T1 )
	n_C  <- n - n_T1
	derive_alpha <- rep(NA, n)
	derive_Sigma <- matrix(NA, nrow = d*d, ncol = n)
	derive_Mu    <- matrix(NA, nrow = d, ncol = n)
	p <- rstable(N, alpha/2, 1, (cos(pi*alpha/4))^(2/alpha), 0, 1)
	C <- 1/sqrt( (4*pi)^d*det(Sigma) )
	T3 <- 6*( (d + alpha)/(2*alpha) )^(1/alpha)
	Sigma_inv <- solve(Sigma)
	node   <- c(-0.996893,-0.983668,-0.960021,-0.926200,-0.882560,-0.829565,-0.767777,-0.697850,-0.620526,
			      -0.536624,-0.447033,-0.352704,-0.254636,-0.153869,-0.051471,
			       0.051471, 0.153869, 0.254636, 0.352704, 0.447033, 0.536624, 0.620526, 0.697850, 0.767777,
			       0.829565, 0.882560, 0.926200, 0.960021, 0.983668, 0.996893)
	weight <- c( 0.007968, 0.018466, 0.028784, 0.038799, 0.048402, 0.057493, 0.065974, 0.073755, 0.080755,
			       0.086899, 0.092122, 0.096368, 0.099593, 0.101762, 0.102852,
			       0.102852, 0.101762, 0.099593, 0.096368, 0.092122, 0.086899, 0.080755, 0.073755, 0.065974,
			       0.057493, 0.048402, 0.038799, 0.028784, 0.018466, 0.007968)
	fun1 <- function(u, a, d, d0, M) M*(M*u)^d*exp(-d0*u^2*M^2/4)*a*(M*u)^(a-1)*exp(-u^a*M^a)*( 1/a+log(M*u) )
	fun2 <- function(u, a, d, d0, M) M*(M*u)^d*exp(-d0*u^2*M^2/4)*a*(M*u)^(a-1)*exp(-u^a*M^a)
	fun5 <- function(x, a, d, y) (a + d - 1)*log(x) - x^a - y*x^2/4 
	fun6 <- function(x, a, d, y) (a + d - 1)/x - a*x^(a-1) - y*x/2 
	Sum  <- rep(0, M)
	Dim_sigma <- d*(d + 1)/2
	S      <- rep( 0, d + Dim_sigma + 1 )
	S_name <- rep(NA, length(S) )
	seq_Dim_sigma <- rep(NA, Dim_sigma)
	index  <- 1
	pdf_T1 <- C/pi*sapply(1:n_T1, function(i) sum( (-1)^(K2 - 1)*exp( lgamma( K2*alpha/2 + 1 ) - 
				lgamma(K2 + 1) + lgamma( (K2*alpha + d)/2 ) )*sin(K2*pi*alpha/2)*( 4/D[ index_T1[i] ] )^( (K2*
				alpha + d)/2 ) ) )
	E1[ index_T1 ] <- C/pi*sapply(1:n_T1, function(i) sum( (-1)^(K2 - 1)*exp( lgamma( K2*alpha/2 + 1 ) -
						lgamma(K2 + 1) + lgamma((K2*alpha+d + 2)/2))*sin(K2*pi*alpha/2)*(4/D[ index_T1[i] ])^
						( (K2*alpha + d + 2)/2) ) )/pdf_T1
	derive_alpha[ index_T1 ] <- C/pi*sapply(1:n_T1, function(i) sum( (-1)^(K2 - 1)*exp( lgamma( K2*alpha/2 +
								1 ) - lgamma(K2) + lgamma(K2*alpha/2 + d/2) )*( D[ index_T1[i] ]/4 )^(-K2*
								alpha/2 - d/2)/2*( sin(K2*pi*alpha/2 )*digamma(K2*alpha/2 + 1  ) + pi*
								cos( K2*pi*alpha/2 )+ sin(K2*pi*alpha/2)*(digamma(K2*alpha/2 + d/2) - 
								log( D[ index_T1[i] ]/4 ) ) ) ) )/pdf_T1
	pdf_C <- sapply(1:n_C, function(i) sum( p^(-d/2)*exp( -0.25*D[ index_C[i] ]/p ), na.rm = TRUE ) )
	E1[ index_C ] <- sapply(1:n_C, function(i) sum( p^(-d/2 - 1)*exp( -0.25*D[ index_C[i] ]/p ), na.rm = TRUE ) )/ pdf_C
		if(n_C > 0){
			for(i in 1:n_C)
			{
				delta <- mahalanobis( c(X[ index_C[i], ] - Mu) , rep(0, d), Sigma )
				dis   <- norm( solve( t( chol(Sigma) ) )%*%c( X[ index_C[i], ] - Mu ), type = "2" )
				E <- rexp(1)
				xi    <- delta/E
					if( alpha <= 1.2 )
					{
						Sum <- rep(0, M)
						for(r in 1:M)
						{
							if( stochastic == TRUE )
							{
								up <- exp( -d/2 + d/2*log(2*d/xi) )
								j <- 1
								while (j < 2)
								{
									w  <- rweibull(1, shape = alpha, scale = 1)
									ex <- w^d*exp( -xi*w^2/4 )
										if ( runif(1) < ex/up )
										{
											w_y <- w
											j <- j + 1
										}
								}
								Sum[r] <- 1/alpha + log( w_y ) - w_y^alpha*log( w_y )
							}else{
								T2     <- sqrt(-4*log(epsilon))/xi
								Sum[r] <- sum( weight/2*fun1( (node + 1)/2, alpha, d, D[ index_C[i] ], T2 ) )/
											sum( weight/2*fun2( (node + 1)/2, alpha, d, D[ index_C[i] ], T2 ) )
							}
						}
					derive_alpha[ index_C[i] ] <- mean( Sum )
					}
					else 
					{
						derive_alpha[ index_C[i] ] <- -sum( weight*(T3*(node + 1)/2)^(d/2 + alpha)*log(T3*(node + 1)/2)*
														exp( -(T3*(node + 1)/2)^alpha )*besselJ(dis*T3*(node + 1)/2, d/2 - 1) )/
														sum( weight*(T3*(node + 1)/2)^(d/2 )*exp( -(T3*(node + 1)/2)^alpha )*
														besselJ(dis*T3*(node + 1)/2, d/2 - 1) )
					}
				}
			}
				Fisher <- matrix(0, Dim, Dim)
				for(i in 1:n)
				{
					derive_Mu[,i] <- Sigma_inv%*%c(X[i, ] - Mu)%*%E1[i]/2
					derive_Sigma[, i] <- - Sigma_inv/2 + Sigma_inv%*%( c(X[i, ] - Mu)%o%c(X[i, ] - Mu)*
										E1[i] )%*%Sigma_inv/4
					S <- c( derive_alpha[i], c( arrange_sigma( matrix( derive_Sigma[, i], d, d) )$y ),
						c( derive_Mu[, i] ) )
					Fisher <- Fisher + S%o%S
				}
					for(i in 1:d)
					{
						for(j in i:d)
						{
							seq_Dim_sigma[ index ] <- paste(i, sep = "", j)
							index <- index + 1
						}
					}
    range_alpha <- seq( 1 , 1 )
    range_sigma <- seq( 2, 1 + Dim_sigma )
    range_mu    <- seq( 2 + Dim_sigma, 1 + Dim_sigma + d )
    S_name[ range_alpha ] <- paste( expression( alpha )                          )
    S_name[ range_sigma ] <- paste( expression( sigma ), sep = "", seq_Dim_sigma )
    S_name[ range_mu    ] <- paste( expression( mu    ), sep = "", 1:d           )
    colnames( Fisher )  <- S_name 	  
    rownames( Fisher )  <- S_name 
    out_std <- sqrt( diag( solve( Fisher, tol = 10e-100 ) ) )
return( list( std_error = out_std, OFI = Fisher ) )
}
