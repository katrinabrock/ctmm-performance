
devtools::load_all('../ctmm')
Sys.setenv(RSTUDIO_PANDOC=pandoc::pandoc_locate())
my_profile <- profvis::profvis({
pure_r_optim <- function(dt, K, tau, sigma, Omega2, f, nu, TT, DIM){
  fix_dim <- function(Green, Sigma){
  #browser()
  # fix the dimension of the filter
  if(DIM==1) # 1D filter
  { Sigma <- sigma * Sigma 
  return(list(Green=Green, Sigma=Sigma))
  }
  # 2D filter
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
  return(list(Green=Green, Sigma=Sigma))
  }
  # K = 0, 1, 2
  # K = 0 or 1 # IID-BM-OU
  # K = 2 # IOU-OUF-OUO
  tau_inf <- tau[1] == Inf
  dt_leq_inf <- dt < Inf  
  
  if (K == 1 & tau_inf & dt_leq_inf) {
    # BM
    Green <- array(1,c(1,1))
    # absorbing 1/tau into sigma # VAR -> Diffusion
    Sigma <- array(2*dt,c(1,1))
    return(fix_dim(Green, Sigma))
  } 
  if (K == 1 & tau_inf ) { # ! dt_leq_inf
    Green <- array(0,c(1,1))
    Sigma <- array(2*dt ,c(1,1))
    return(fix_dim(Green, Sigma))
  }
  if (K == 1 & dt_leq_inf){ # ! tau_inf
    # (BM,OU,IID]
    dtau <- dt/tau
    c0 <- exp(-dtau)
    Green <- array(c0,c(1,1))
    Sigma <- array(dexp2(dtau,Exp=c0),c(1,1))
    return(fix_dim(Green, Sigma))
  } 
  if (K == 1 | K == 0) { # ! dt_leq_inf ! tau_inf
    # IID limit
    Green <- array(0,c(1,1))
    Sigma <- array(1,c(1,1))
    return(fix_dim(Green, Sigma))
  } 
  if (K == 2  & tau_inf & dt_leq_inf) {
    #IID limit
    #IOU
    fdt <- f*dt
    dtau <- dt/tau[2]
    Exp <- exp(-dtau)
    DExp <- dexp1(dtau,Exp) # 1-exp(dt/tau[2])

    Green <- array(c(
      1, 0, tau[2] * DExp, Exp
    ), c(2,2))

    # remember that sigma is D=sigma/tau[1]
    DExp2 <- DExp^2 # (1-exp(-dt/tau[2]))^2

    Sigma <- array( c(0, 0, 0, 0) , c(2,2))
    Sigma[1,1] <- clamp( 2*dt - tau[2]*(2*DExp+DExp2) ,0,Inf) # does this still underflow?
    Sigma[2,2] <- dexp2(dtau,Exp)/tau[2]
    Sigma[c(2,3)] <- clamp(DExp2,0, sqrt(Sigma[1,1]*Sigma[2,2]) ) 
    # 0 at dt=Inf
    return(fix_dim(Green, Sigma))
  }
  if (K ==2  & tau_inf) { # ! dt_leq_inf
    #IID limit
    Green <- array( c(0,0 , 0,0) , c(2,2))

    # IOU
    Exp <- exp(-dtau)
    DExp <- dexp1(dtau,Exp) # 1-exp(dt/tau[2])
    DExp2 <- DExp^2 # (1-exp(-dt/tau[2]))^2

    # does this still underflow?
    Sigma <- array(c(
      clamp( 2*dt - tau[2]*(2*DExp+DExp2) ,0,Inf), 0, 0, dexp2(dtau,Exp)/tau[2]
    ) , c(2,2))
    # END IOU
    return(fix_dim(Green, Sigma))
  }
  if (K ==2  & !tau_inf &!dt_leq_inf) {
    #IID limit
    Green <- rbind( c(0,0) , c(0,0) )
    Sigma <- rbind( c(1,0) , c(0,Omega2) )
    return(fix_dim(Green, Sigma))
  }

  tau1_g_tau2 <-  tau[1] > tau[2] 
  nudt <- nu*dt
  nudt_over_thresh <- nudt>0.8813736

  if (K == 2 & dt_leq_inf & tau1_g_tau2 & nudt_over_thresh){ # ! tau_eq_inf
    # function representation choice
    dtau <- dt/tau
    dift <- diff(tau)
    Exp0 <- exp(-dtau)
    Exp <- Exp0/dift
    c0 <- diff(Exp*tau)
    c1 <- -diff(Exp)
    c2 <- diff(Exp/tau)

    Green <- array(c(
      c0, c1, -c1/Omega2, -c2/Omega2
    ) , c(2,2))

    # initially canceling terms
    dift2 <- dift^2
    T2 <- tau^2
    S1 <- dexp2(dtau[1],Exp0[1])
    S2 <- dexp2(dtau[2],Exp0[2])
    S12 <- 2*tau[1]*tau[2]*dexp1(f*dt,Exp0[1]*Exp0[2])

    # initially vanishing terms
    c12 <- c1^2

    Sigma <- array(c(
      (T2[1]*S1 - S12 + T2[2]*S2)/dift2 - c12/Omega2,
      TT*c12,
      TT*c12,
      (T2[2]*S1 - S12 + T2[1]*S2)/dift2 * Omega2 - c12
    ), c(2,2))
    return(fix_dim(Green, Sigma))
  }

  if (K == 2 & dt_leq_inf & tau1_g_tau2 ){ # ! tau_eq_inf ! nudt_over_thresh
    # hyperbolic-trigonometric function
    Sigma <- array( c(1,0, 0, Omega2) , c(2,2))
    fdt <- f*dt

    # function representation choice
    Exp <- exp(-fdt)
    Sin0 <- sinh(nudt)
    Sinc0 <- sinch(nudt,Sin0)
    Cos0 <- cosh(nudt)
    SincE <- Sinc0*Exp
    CosE <- Cos0*Exp

    c0 <- CosE + fdt*SincE
    c1 <- -(Omega2*dt)*SincE
    c2 <- -Omega2*(CosE - fdt*SincE)
    # end function representation

    Green <- array(c(
      c0, c1, -c1/Omega2, -c2/Omega2
    ) , c(2,2))

    # initially canceling terms
    CROSS <- fdt*Sinc0*Exp
    OUTER <- Cos0^2*dexp2(fdt,Exp) - CROSS^2
    CROSS <- 2*Cos0*Exp*CROSS
    Sin2 <- Sin0^2
    Sigma[1,1] <- OUTER - Sin2 - CROSS
    Sigma[2,2] <- (OUTER - Sin2 + CROSS) * Omega2

    c12 <- c1^2

    Sigma <- array(c(
      Sigma[1,1] - c12/Omega2,
      TT*c12,
      TT*c12,
      Sigma[2,2] - c12
    ), c(2,2))

    return(fix_dim(Green, Sigma))

  }
  if (K == 2 & dt_leq_inf ){ # ! tau_eq_inf ! nudt_over_thresh !tau1_g_tau2 
    #IID limit
    Sigma <- array( c(1,0, 0, Omega2) , c(2,2))
    fdt <- f*dt
    # function representation choice
    nudt <- nu*dt
    EXP <- (tau[1]>tau[2] && nudt>0.8813736)
    Exp <- exp(-fdt)
    # trigonometric
    Sin0 <- sin(nudt)
    Sinc0 <- sinc(nudt,Sin0)
    Cos0 <- cos(nudt)
    SincE <- Sinc0*Exp
    CosE <- Cos0*Exp

    c0 <- CosE + fdt*SincE
    c1 <- -(Omega2*dt)*SincE
    c2 <- -Omega2*(CosE - fdt*SincE)
    # end function representation

    Green <- array(c(
      c0, c1, -c1/Omega2, -c2/Omega2
    ) , c(2,2))

    CROSS <- fdt*Sinc0*Exp
    OUTER <- Cos0^2*dexp2(fdt,Exp) - CROSS^2
    CROSS <- 2*Cos0*Exp*CROSS
    Sin2 <- Sin0^2

    Sigma[1,1] <- OUTER + Sin2 - CROSS
    Sigma[2,2] <- (OUTER + Sin2 + CROSS) * Omega2

    # initially vanishing terms
    c12 <- c1^2

    Sigma <- array(c(
      Sigma[1,1] - c12/Omega2,
      TT*c12,
      TT*c12,
      Sigma[2,2] - c12
    ), c(2,2))
    # end OUF/OUO

    return(fix_dim(Green, Sigma))
  }

}

  for(i in 1:1000) run_over_inputs(pure_r_optim)
})
htmlwidgets::saveWidget(my_profile, "../ctmm/profile_refactored.html")
