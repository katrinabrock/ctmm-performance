library(data.table)
library(microbenchmark)
#source('../ctmm/R/1.R')
devtools::load_all('../ctmm')
inputs_dt <- fread('data/langevin_in.csv')

orig <- function(dt, K, tau, sigma, Omega2, f, nu, TT, DIM){
  if(K<=1) # IID-BM-OU
  {
    # IID limit
    Green <- array(0,c(1,1))
    Sigma <- array(1,c(1,1))

    if(K)
    {
      if(tau[1]==Inf) # BM
      {
        if(dt<Inf) { Green[1,1] <- 1 }
        # absorbing 1/tau into sigma # VAR -> Diffusion
        Sigma[1,1] <- 2*dt
      }
      else if(dt<Inf) # (BM,OU,IID]
      {
        dtau <- dt/tau
        c0 <- exp(-dtau)
        Green[1,1] <- c0
        Sigma[1,1] <- dexp2(dtau,Exp=c0)
      }
    } # >IID
  } # IID-BM-OU
  else if(K==2) # IOU-OUF-OUO
  {
    #IID limit
    Green <- rbind( c(0,0) , c(0,0) )
    Sigma <- rbind( c(1,0) , c(0,Omega2) )

    fdt <- f*dt

    if(tau[1]==Inf) # IOU
    {
      dtau <- dt/tau[2]
      Exp <- exp(-dtau)
      DExp <- dexp1(dtau,Exp) # 1-exp(dt/tau[2])

      if(dt<Inf)
      {
        Green[1,1] <- 1
        Green[1,2] <- tau[2]*DExp
        Green[2,2] <- Exp
      }

      # remember that sigma is D=sigma/tau[1]
      DExp2 <- DExp^2 # (1-exp(-dt/tau[2]))^2
      Sigma[1,1] <- clamp( 2*dt - tau[2]*(2*DExp+DExp2) ,0,Inf) # does this still underflow?
      Sigma[2,2] <- dexp2(dtau,Exp)/tau[2]
      if(dt<Inf) { Sigma[c(2,3)] <- clamp(DExp2,0, sqrt(Sigma[1,1]*Sigma[2,2]) ) } # how does this get so far off?
      # 0 at dt=Inf
    } # END IOU
    else if(dt<Inf) # (IOU,OUF/OUO,IID]
    {
      # function representation choice
      nudt <- nu*dt
      EXP <- (tau[1]>tau[2] && nudt>0.8813736)
      if(EXP) # exponential functions
      {
        dtau <- dt/tau
        dift <- diff(tau)
        Exp0 <- exp(-dtau)
        Exp <- Exp0/dift
        c0 <- diff(Exp*tau)
        c1 <- -diff(Exp)
        c2 <- diff(Exp/tau)
      }
      else # trigonometric and hyperbolic-trigonometric functions
      {
        Exp <- exp(-fdt)

        if(tau[1]>tau[2]) # hyperbolic-trigonometric
        {
          Sin0 <- sinh(nudt)
          Sinc0 <- sinch(nudt,Sin0)
          Cos0 <- cosh(nudt)
        }
        else # trigonometric
        {
          Sin0 <- sin(nudt)
          Sinc0 <- sinc(nudt,Sin0)
          Cos0 <- cos(nudt)
        }
        SincE <- Sinc0*Exp
        CosE <- Cos0*Exp

        c0 <- CosE + fdt*SincE
        c1 <- -(Omega2*dt)*SincE
        c2 <- -Omega2*(CosE - fdt*SincE)
      } # end function representation

      Green[1,1] <- c0
      Green[2,1] <- c1
      Green[1,2] <- -c1/Omega2
      Green[2,2] <- -c2/Omega2

      # initially canceling terms
      if(EXP)
      {
        dift2 <- dift^2
        T2 <- tau^2
        S1 <- dexp2(dtau[1],Exp0[1])
        S2 <- dexp2(dtau[2],Exp0[2])
        S12 <- 2*tau[1]*tau[2]*dexp1(fdt,Exp0[1]*Exp0[2])
        Sigma[1,1] <- (T2[1]*S1 - S12 + T2[2]*S2)/dift2
        Sigma[2,2] <- (T2[2]*S1 - S12 + T2[1]*S2)/dift2 * Omega2
      }
      else
      {
        CROSS <- fdt*Sinc0*Exp
        OUTER <- Cos0^2*dexp2(fdt,Exp) - CROSS^2
        CROSS <- 2*Cos0*Exp*CROSS
        Sin2 <- Sin0^2

        if(tau[1]>tau[2])
        {
          Sigma[1,1] <- OUTER - Sin2 - CROSS
          Sigma[2,2] <- (OUTER - Sin2 + CROSS) * Omega2
        }
        else
        {
          Sigma[1,1] <- OUTER + Sin2 - CROSS
          Sigma[2,2] <- (OUTER + Sin2 + CROSS) * Omega2
        }
      }

      # initially vanishing terms
      c12 <- c1^2
      Sigma[1,1] <- Sigma[1,1] - c12/Omega2
      Sigma[c(2,3)] <- TT*c12
      Sigma[2,2] <- Sigma[2,2] - c12
    } # end OUF/OUO
  }

  # fix the dimension of the filter
  if(DIM==1) # 1D filter
  { Sigma <- sigma * Sigma }
  else # 2D filter
  {
    if(length(sigma)==1) { sigma <- diag(sigma,DIM) }

    K <- max(1,K)

    Sigma <- outer(Sigma,sigma) # (k,k,d,d)
    Sigma <- aperm(Sigma,c(1,3,2,4)) # (k,k,d,d) -> (k,d,k,d)
    dim(Sigma) <- c(K*DIM,K*DIM)

    # BM/IOU prior fix
    NAN <- is.nan(Sigma)
    if(any(NAN)) { Sigma[NAN] <- 0 }

    Green <- outer(Green,diag(DIM)) # (k,k,d,d)
    Green <- aperm(Green,c(1,3,2,4)) # (k,d,k,d)
    dim(Green) <- c(K*DIM,K*DIM)
  }

  return(list(Green=Green, Sigma=Sigma))

}