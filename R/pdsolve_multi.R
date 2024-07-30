orig_pd.solve <- function(M,sym=TRUE,semi=TRUE,...)
{
  NAMES <- rev( dimnames(M) ) # dimnames for inverse matrix
  DIM <- dim(M)
  if(is.null(DIM))
  {
    M <- as.matrix(M)
    DIM <- dim(M)
  }
  if(DIM[1]==0) { return(M) }

  # check for Inf & invert those to 0 (and vice versa)
  INF <- diag(M)==Inf
  ZERO <- diag(M)<=0 & sym
  if(any(INF) || any(ZERO))
  {
    # 1/Inf == 0 # correlations not accounted for
    if(any(INF)) { M[INF,] <- M[,INF] <- 0 }

    # 1/0 == Inf
    if(any(ZERO))
    {
      M[ZERO,] <- M[,ZERO] <- 0
      diag(M)[ZERO] <- Inf
    }

    # regular inverse of remaining dimensions
    REM <- !(INF|ZERO)
    if(any(REM)) { M[REM,REM] <- pd.solve(M[REM,REM,drop=FALSE],sym=sym,semi=semi,...) }

    dimnames(M) <- NAMES
    return(M)
  }

  if(semi)
  {
    if(DIM[1]==1)
    {
      M <- matrix(1/M,c(1,1))

      dimnames(M) <- NAMES
      return(M)
    }
    if(DIM[1]==2)
    {
      DET <- M[1,1]*M[2,2]-M[1,2]*M[2,1]
      if(DET<=0) { return(diag(Inf,2)) } # force positive definite / diagonal
      SWP <- M[1,1] ; M[1,1] <- M[2,2] ; M[2,2] <- SWP
      M[1,2] <- -M[1,2]
      M[2,1] <- -M[2,1]
      M <- M/DET

      dimnames(M) <- NAMES
      return(M)
    }
  }

  # symmetrize
  if(sym) { M <- He(M) }

  # rescale
  W <- abs(diag(M))
  W <- sqrt(W)
  ZERO <- W<=.Machine$double.eps
  if(any(ZERO)) # do not divide by zero or near zero
  {
    if(any(!ZERO))
    { W[ZERO] <- min(W[!ZERO]) } # do some rescaling... assuming axes are similar
    else
    { W[ZERO] <- 1 } # do no rescaling
  }
  W <- W %o% W

  # now a correlation matrix that is easier to invert
  M <- M/W

  # try ordinary inverse
  M.try <- try(qr.solve(M,tol=0),silent=TRUE)
  # fall back on decomposition
  if( class(M.try)[1] == "matrix" && (!sym || all(diag(M.try>=0))) )
  { M <- M.try }
  else
  { M <- PDfunc(M,func=function(m){1/m},sym=sym,semi=semi,...) }

  # back to covariance matrix
  M <- M/W

  # symmetrize
  if(sym) { M <- He(M) }

  dimnames(M) <- NAMES
  return(M)
}

orig_multi.pd.solve <- function(arr,n,d){
  return(vapply(1:n,function(i){orig_pd.solve(arr[i,,])},diag(d))) # (OBS*DIM,OBS*DIM,n)
}

vect_pd.solve <- function(M,sym=TRUE,semi=TRUE,...)
{
  NAMES <- rev( dimnames(M) ) # dimnames for inverse matrix
  DIM <- dim(M)
  if(is.null(DIM))
  {
    n <- length(M)
    DIM <- c(n,min(n, 1))
    M <- array(M, DIM)
  }
  if(DIM[2]==0) { return(M) }

  # check for Inf & invert those to 0 (and vice versa)
  INF <- diag(M)==Inf
  ZERO <- diag(M)<=0 & sym
  if(any(INF) || any(ZERO))
  {
    # 1/Inf == 0 # correlations not accounted for
    if(any(INF)) { M[INF,] <- M[,INF] <- 0 }

    # 1/0 == Inf
    if(any(ZERO))
    {
      M[ZERO,] <- M[,ZERO] <- 0
      diag(M)[ZERO] <- Inf
    }

    # regular inverse of remaining dimensions
    REM <- !(INF|ZERO)
    if(any(REM)) { M[REM,REM] <- pd.solve(M[REM,REM,drop=FALSE],sym=sym,semi=semi,...) }

    dimnames(M) <- NAMES
    return(M)
  }

  if(semi)
  {
    if(DIM[2]==1)
    {
      M <- matrix(1/M,c(1,1))

      dimnames(M) <- NAMES
      return(M)
    }
    if(DIM[2]==2)
    {
      DET <- M[,1,1]*M[,2,2]-M[,1,2]*M[,2,1]
      if(DET<=0) { return(diag(Inf,2)) } # force positive definite / diagonal
      SWP <- M[,1,1] ; M[,1,1] <- M[,2,2] ; M[,2,2] <- SWP
      M[,1,2] <- -M[,1,2]
      M[,2,1] <- -M[,2,1]
      M <- M/DET

      dimnames(M) <- NAMES
      return(M)
    }
  }

  # symmetrize
  if(sym) { M <- He(M) }

  # rescale
  W <- abs(diag(M))
  W <- sqrt(W)
  ZERO <- W<=.Machine$double.eps
  if(any(ZERO)) # do not divide by zero or near zero
  {
    if(any(!ZERO))
    { W[ZERO] <- min(W[!ZERO]) } # do some rescaling... assuming axes are similar
    else
    { W[ZERO] <- 1 } # do no rescaling
  }
  W <- W %o% W

  # now a correlation matrix that is easier to invert
  M <- M/W

  # try ordinary inverse
  M.try <- try(qr.solve(M,tol=0),silent=TRUE)
  # fall back on decomposition
  if( class(M.try)[1] == "matrix" && (!sym || all(diag(M.try>=0))) )
  { M <- M.try }
  else
  { M <- PDfunc(M,func=function(m){1/m},sym=sym,semi=semi,...) }

  # back to covariance matrix
  M <- M/W

  # symmetrize
  if(sym) { M <- He(M) }

  dimnames(M) <- NAMES
  return(M)
}

inv <- function(M,x,y){
  return(as.vector(1/M))
}

load('tests/testthat/testdata/test_multi_pd_solve.rda')

testthat::expect_identical(
    #lapply(args, function(a) do.call(orig_multi.pd.solve, a)),
    lapply(args, function(a) do.call(inv, a)),
    expected
)
