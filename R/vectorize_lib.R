fix_dim <- function(Green, Sigma, sigma, DIM,K){
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

langevin_fn_factory <- function(K, tau, sigma, Omega2, f, nu, TT, DIM){
# K = 0, 1, 2
# K = 0 or 1 # IID-BM-OU
# K = 2 # IOU-OUF-OUO
if (K == 0) {
  return(function(dt, dt_leq_inf, nudt_over_thresh, n = length(dt)) { # ! dt_leq_inf ! tau_inf
    # IID limit
    Green <- array(0,c(n, 1,1))
    Sigma <- array(1,c(n, 1,1))
    return(list(Green=Green, Sigma=Sigma))
  })
}
tau_inf <- tau[1] == Inf
if (K == 1 & tau_inf ) {
  return(function(dt, dt_leq_inf, nudt_over_thresh, n = length(dt)) {
    Green <- array(as.integer(dt_leq_inf),c(n,1,1))
    Sigma <- array(2*dt ,c(n,1,1))
    return(list(Green=Green, Sigma=Sigma))
  })
}
if (K == 1) { # ! tau_inf
  return(function(dt, dt_leq_inf, nudt_over_thresh, n = length(dt)) {
    # reminder: may be worth precalculating 
    # I think K=length(tau). If that's wrong, this will break.
    dtau <- dt/tau 
    c0 <- exp(-dtau)
    Green <- array(
      ifelse(
        dt_leq_inf,
        # (BM,OU,IID]
        c0,
        # IID limit
        0
      ),
      c(n,1,1)
    )
    Sigma <- array(
      ifelse(
        dt_leq_inf,
        # (BM,OU,IID]
        dexp2(dtau,Exp=c0),
        # IID limit
        1
      ),
      c(n,1,1)
    )
    return(list(Green=Green, Sigma=Sigma))
  })
}
if (K == 2 & tau_inf) {
  return(function(dt, dt_leq_inf, nudt_over_thresh, n = length(dt)) {
    fdt <- f*dt
    dtau <- dt/tau[2]
    Exp <- exp(-dtau)
    DExp <- dexp1(dtau,Exp) # 1-exp(dt/tau[2])
    DExp2 <- DExp^2 # (1-exp(-dt/tau[2]))^2
    Green <- array(
      c(
         ifelse(dt_leq_inf, 1, 0),
         rep(0,n),
         ifelse(dt_leq_inf, tau[2] * DExp, 0),
         ifelse(dt_leq_inf, Exp, 0)
      ),
      c(n,2,2)        
    )

    s11 <- clamp( 2*dt - tau[2]*(2*DExp+DExp2) ,0,Inf) # does this still underflow?
    s22 <- dexp2(dtau,Exp)/tau[2]
    s23 <- clamp(DExp2,0, sqrt(s11*s22) ) 

    Sigma <- array(
      c(
        s11,
        ifelse(dt_leq_inf, s23, 0),
        ifelse(dt_leq_inf, s23, 0),
        s22
      ),        
      c(n,2,2)
    )
    return(list(Green=Green, Sigma=Sigma))
  })
}

tau1_g_tau2 <-  tau[1] > tau[2] 
if (K == 2 & tau1_g_tau2) { # !tau_inf
  return(function(dt, dt_leq_inf, nudt_over_thresh, n = length(dt)) {
    #TODO: only calculate variables for relevant case

    ### Variables for !nudt_over_thesh
    # function representation choice
    # n by k
    dtau <- array(dt, c(n,K))/array(rep(tau, each = n), c(n,K))
    #1
    dift <- diff(tau)
    # n by k
    Exp0 <- exp(-dtau)
    # n by k
    Exp_dift <- Exp0/dift

    ### Variables for !nudt_over_thesh
    # hyperbolic-trigonometric function
    # nx1
    fdt <- f*dt
    nudt <- nu*dt

    # function representation choice
    # nx1
    Exp_fdt <- exp(-fdt)
    Sin0 <- sinh(nudt)
    Sinc0 <- sinch(nudt,Sin0)
    Cos0 <- cosh(nudt)
    SincE <- Sinc0*Exp_fdt
    CosE <- Cos0*Exp_fdt

    # n by 1
    c0 <- ifelse(nudt_over_thresh,
      Exp_dift[,2]*tau[2] - Exp_dift[,1]*tau[1],
      CosE + fdt*SincE
    )
    c1 <-  ifelse(nudt_over_thresh,
      Exp_dift[,1] - Exp_dift[,2],
      -(Omega2*dt)*SincE
    )
    c2 <- ifelse(nudt_over_thresh,
      Exp_dift[,2]/tau[2] - Exp_dift[,1]/tau[1],
      -Omega2*(CosE - fdt*SincE)
    )

    # end function representation

    Green <- array(ifelse(rep(dt_leq_inf, 4), c(
      c0, c1, -c1/Omega2, -c2/Omega2
    ), 0) , c(n,2,2))

    ## Variables for nudt_over_thresh 
    # initially canceling terms
    # 1
    dift2 <- dift^2
    # k
    T2 <- tau^2
    # n
    S1 <- dexp2(dtau[,1],Exp0[,1])
    S2 <- dexp2(dtau[,2],Exp0[,2])
    S12 <- 2*tau[1]*tau[2]*dexp1(f*dt,Exp0[,1]*Exp0[,2])

    ## Variables for !nudt_over_thresh 
    CROSS <- fdt*Sinc0*Exp_fdt
    OUTER <- Cos0^2*dexp2(fdt,Exp_fdt) - CROSS^2
    CROSS <- 2*Cos0*Exp_fdt*CROSS
    Sin2 <- Sin0^2

    # initially vanishing terms
    s1 <- ifelse(
      nudt_over_thresh,
      (T2[1]*S1 - S12 + T2[2]*S2)/dift2,
      OUTER - Sin2 - CROSS
    )

    s2 <- ifelse(
      nudt_over_thresh,
      (T2[2]*S1 - S12 + T2[1]*S2)/dift2 * Omega2,
      (OUTER - Sin2 + CROSS) * Omega2
    )

    c12 <- c1^2

    Sigma <- array(ifelse(rep(dt_leq_inf, 4), c(
      s1 - c12/Omega2,
      TT*c12,
      TT*c12,
      s2 - c12
    # !dt_leq_inf case
    ), rep(c(1,0,0,Omega2), each=n)), c(n,2,2))

    return(list(Green=Green, Sigma=Sigma))
  })
}
if (K == 2) { # ! tau1_g_tau2
  return(function(dt, dt_leq_inf, nudt_over_thresh, n = length(dt)) {
    # ! tau_eq_inf ! nudt_over_thresh !tau1_g_tau2 
    fdt <- f*dt
    # function representation choice
    nudt <- nu*dt
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

    Green <- array(ifelse(rep(dt_leq_inf, 4), c(
      c0, c1, -c1/Omega2, -c2/Omega2
    ),0) , c(n,2,2))

    CROSS <- fdt*Sinc0*Exp
    OUTER <- Cos0^2*dexp2(fdt,Exp) - CROSS^2
    CROSS <- 2*Cos0*Exp*CROSS
    Sin2 <- Sin0^2

    s1 <- OUTER + Sin2 - CROSS
    s2 <- (OUTER + Sin2 + CROSS) * Omega2

    # initially vanishing terms
    c12 <- c1^2

    Sigma <- array(ifelse(rep(dt_leq_inf, 4), c(
      s1 - c12/Omega2,
      TT*c12,
      TT*c12,
      s2 - c12
    ), rep(c(1,0,0,Omega2), each=n)), c(n,2,2))
    # end OUF/OUO
    return(list(Green=Green, Sigma=Sigma))
  })
}
}

no_factory <- function(dt, K, tau, sigma, Omega2, f, nu, TT, DIM){

dt_leq_inf <- dt < Inf
nudt <- nu*dt
nudt_over_thresh <- nudt>0.8813736
n <- length(dt)

# K = 0, 1, 3
# K = 0 or 1 # IID-BM-OU
# K = 2 # IOU-OUF-OUO
if (K == 0) {
    # IID limit
    Green <- array(0,c(n, 1,1))
    Sigma <- array(1,c(n, 1,1))
    return(list(Green=Green, Sigma=Sigma))
}
tau_inf <- tau[1] == Inf
if (K == 1 & tau_inf ) {
    Green <- array(as.integer(dt_leq_inf),c(n,1,1))
    Sigma <- array(2*dt ,c(n,1,1))
    return(list(Green=Green, Sigma=Sigma))
}
if (K == 1) { # ! tau_inf
    # reminder: may be worth precalculating 
    # I think K=length(tau). If that's wrong, this will break.
    dtau <- dt/tau 
    c0 <- exp(-dtau)
    Green <- array(
      ifelse(
        dt_leq_inf,
        # (BM,OU,IID]
        c0,
        # IID limit
        0
      ),
      c(n,1,1)
    )
    Sigma <- array(
      ifelse(
        dt_leq_inf,
        # (BM,OU,IID]
        dexp2(dtau,Exp=c0),
        # IID limit
        1
      ),
      c(n,1,1)
    )
    return(list(Green=Green, Sigma=Sigma))
}
if (K == 2 & tau_inf) {
    fdt <- f*dt
    dtau <- dt/tau[2]
    Exp <- exp(-dtau)
    DExp <- dexp1(dtau,Exp) # 1-exp(dt/tau[2])
    DExp2 <- DExp^2 # (1-exp(-dt/tau[2]))^2
    Green <- array(
      c(
         ifelse(dt_leq_inf, 1, 0),
         rep(0,n),
         ifelse(dt_leq_inf, tau[2] * DExp, 0),
         ifelse(dt_leq_inf, Exp, 0)
      ),
      c(n,2,2)        
    )

    s11 <- clamp( 2*dt - tau[2]*(2*DExp+DExp2) ,0,Inf) # does this still underflow?
    s22 <- dexp2(dtau,Exp)/tau[2]
    s23 <- clamp(DExp2,0, sqrt(s11*s22) ) 

    Sigma <- array(
      c(
        s11,
        ifelse(dt_leq_inf, s23, 0),
        ifelse(dt_leq_inf, s23, 0),
        s22
      ),        
      c(n,2,2)
    )
    return(list(Green=Green, Sigma=Sigma))
}

tau1_g_tau2 <-  tau[1] > tau[2] 
if (K == 2 & tau1_g_tau2) { # !tau_inf
    #TODO: only calculate variables for relevant case

    ### Variables for !nudt_over_thesh
    # function representation choice
    # n by k
    dtau <- array(dt, c(n,K))/array(rep(tau, each = n), c(n,K))
    #1
    dift <- diff(tau)
    # n by k
    Exp0 <- exp(-dtau)
    # n by k
    Exp_dift <- Exp0/dift

    ### Variables for !nudt_over_thesh
    # hyperbolic-trigonometric function
    # nx1
    fdt <- f*dt
    nudt <- nu*dt

    # function representation choice
    # nx1
    Exp_fdt <- exp(-fdt)
    Sin0 <- sinh(nudt)
    Sinc0 <- sinch(nudt,Sin0)
    Cos0 <- cosh(nudt)
    SincE <- Sinc0*Exp_fdt
    CosE <- Cos0*Exp_fdt

    # n by 1
    c0 <- ifelse(nudt_over_thresh,
      Exp_dift[,2]*tau[2] - Exp_dift[,1]*tau[1],
      CosE + fdt*SincE
    )
    c1 <-  ifelse(nudt_over_thresh,
      Exp_dift[,1] - Exp_dift[,2],
      -(Omega2*dt)*SincE
    )
    c2 <- ifelse(nudt_over_thresh,
      Exp_dift[,2]/tau[2] - Exp_dift[,1]/tau[1],
      -Omega2*(CosE - fdt*SincE)
    )

    # end function representation

    Green <- array(ifelse(rep(dt_leq_inf, 4), c(
      c0, c1, -c1/Omega2, -c2/Omega2
    ), 0) , c(n,2,2))

    ## Variables for nudt_over_thresh 
    # initially canceling terms
    # 1
    dift2 <- dift^2
    # k
    T2 <- tau^2
    # n
    S1 <- dexp2(dtau[,1],Exp0[,1])
    S2 <- dexp2(dtau[,2],Exp0[,2])
    S12 <- 2*tau[1]*tau[2]*dexp1(f*dt,Exp0[,1]*Exp0[,2])

    ## Variables for !nudt_over_thresh 
    CROSS <- fdt*Sinc0*Exp_fdt
    OUTER <- Cos0^2*dexp2(fdt,Exp_fdt) - CROSS^2
    CROSS <- 2*Cos0*Exp_fdt*CROSS
    Sin2 <- Sin0^2

    # initially vanishing terms
    s1 <- ifelse(
      nudt_over_thresh,
      (T2[1]*S1 - S12 + T2[2]*S2)/dift2,
      OUTER - Sin2 - CROSS
    )

    s2 <- ifelse(
      nudt_over_thresh,
      (T2[2]*S1 - S12 + T2[1]*S2)/dift2 * Omega2,
      (OUTER - Sin2 + CROSS) * Omega2
    )

    c12 <- c1^2

    Sigma <- array(ifelse(rep(dt_leq_inf, 4), c(
      s1 - c12/Omega2,
      TT*c12,
      TT*c12,
      s2 - c12
    # !dt_leq_inf case
    ), rep(c(1,0,0,Omega2), each=n)), c(n,2,2))

    return(list(Green=Green, Sigma=Sigma))
}
if (K == 2) { # ! tau1_g_tau2
    # ! tau_eq_inf ! nudt_over_thresh !tau1_g_tau2 
    fdt <- f*dt
    # function representation choice
    nudt <- nu*dt
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

    Green <- array(ifelse(rep(dt_leq_inf, 4), c(
      c0, c1, -c1/Omega2, -c2/Omega2
    ),0) , c(n,2,2))

    CROSS <- fdt*Sinc0*Exp
    OUTER <- Cos0^2*dexp2(fdt,Exp) - CROSS^2
    CROSS <- 2*Cos0*Exp*CROSS
    Sin2 <- Sin0^2

    s1 <- OUTER + Sin2 - CROSS
    s2 <- (OUTER + Sin2 + CROSS) * Omega2

    # initially vanishing terms
    c12 <- c1^2

    Sigma <- array(ifelse(rep(dt_leq_inf, 4), c(
      s1 - c12/Omega2,
      TT*c12,
      TT*c12,
      s2 - c12
    ), rep(c(1,0,0,Omega2), each=n)), c(n,2,2))
    # end OUF/OUO
    return(list(Green=Green, Sigma=Sigma))
}
}
vect_from_orig <- function(dt, K, tau, sigma, Omega2, f, nu, TT, DIM){
  n <- length(dt)
  dt_inf <- !(dt<Inf)
  if(K<=1) # IID-BM-OU
  {
    # IID limit
    Green <- array(0,c(n,1,1))
    Sigma <- array(1,c(n,1,1))

    if(K)
    {
      if(tau[1]==Inf) # BM
      {
        Green[!dt_inf,1,1] <- 1
        # absorbing 1/tau into sigma # VAR -> Diffusion
        Sigma[,1,1] <- 2*dt
      }
      else # (BM,OU,IID]
      {
        dtau <- dt[!dt_inf]/tau
        c0 <- exp(-dtau)
        Green[!dt_inf,1,1] <- c0
        Sigma[!dt_inf,1,1] <- dexp2(dtau,Exp=c0)
      }
    } # >IID
  } # IID-BM-OU
  else if(K==2) # IOU-OUF-OUO
  {
    #IID limit
    Green <- array(0, c(n,2,2))
    Sigma <- array(rep(c(1,0,0,Omega2), each = n), c(n,2,2) )

    fdt <- f*dt

    if(tau[1]==Inf) # IOU
    {
      dtau <- dt/tau[2]
      Exp <- exp(-dtau)
      DExp <- dexp1(dtau,Exp) # 1-exp(dt/tau[2])

      Green[!dt_inf,1,1] <- 1
      Green[!dt_inf,1,2] <- tau[2]*DExp[!dt_inf]
      Green[!dt_inf,2,2] <- Exp[!dt_inf]

      # remember that sigma is D=sigma/tau[1]
      DExp2 <- DExp^2 # (1-exp(-dt/tau[2]))^2
      Sigma[,1,1] <- clamp( 2*dt - tau[2]*(2*DExp+DExp2) ,0,Inf) # does this still underflow?
      Sigma[,2,2] <- dexp2(dtau,Exp)/tau[2]

      # dt<Inf
      s1221 <- clamp(DExp2[!dt_inf],0, sqrt(Sigma[!dt_inf,1,1]*Sigma[!dt_inf,2,2]) ) # how does this get so far off?
      Sigma[!dt_inf,2,1] <- s1221
      Sigma[!dt_inf,1,2] <- s1221
      # 0 at dt=Inf
    } # END IOU
    else # if(dt<Inf) handled by subsetting # (IOU,OUF/OUO,IID]
    {
      # function representation choice
      nudt <- nu*dt[!dt_inf]
      # length(EXP) == sum(!dt_inf)
      EXP <- (tau[1]>tau[2] & nudt>0.8813736)

      c012 <- array(0,c(sum(!dt_inf), 3))

      # EXP exponential functions
      n_exp <- sum(EXP)
      #nxk
      dtau <- array(dt[!dt_inf][EXP], c(n_exp,K))/array(rep(tau, each = n_exp), c(n_exp,K))
      #1
      dift <- diff(tau)
      #nxk
      Exp0 <- exp(-dtau)
      #nxk
      Exp <- Exp0/dift
      #sum(!dt_inf)
      c012[EXP,1] <- Exp[,2]*tau[2] - Exp[,1]*tau[1]
      c012[EXP,2] <- Exp[,1] - Exp[,2] 
      c012[EXP,3] <- Exp[,2]/tau[2] - Exp[,1]/tau[1]

      # !EXP trigonometric and hyperbolic-trigonometric functions
      Exp <- exp(-fdt[!dt_inf][!EXP])

      if(tau[1]>tau[2]) # hyperbolic-trigonometric
      {
        Sin0 <- sinh(nudt[!EXP])
        Sinc0 <- sinch(nudt[!EXP],Sin0)
        Cos0 <- cosh(nudt[!EXP])
      }
      else # trigonometric
      {
        Sin0 <- sin(nudt[!EXP])
        Sinc0 <- sinc(nudt[!EXP],Sin0)
        Cos0 <- cos(nudt[!EXP])
      }
      SincE <- Sinc0*Exp
      CosE <- Cos0*Exp

      c012[!EXP,1] <- CosE + fdt[!dt_inf][!EXP]*SincE
      c012[!EXP,2] <- -(Omega2*dt[!dt_inf][!EXP])*SincE
      c012[!EXP,3] <- -Omega2*(CosE - fdt[!dt_inf][!EXP]*SincE)
      # end function representation

      Green[!dt_inf,1,1] <- c012[,1]
      Green[!dt_inf,2,1] <- c012[,2]
      Green[!dt_inf,1,2] <- -c012[,2]/Omega2
      Green[!dt_inf,2,2] <- -c012[,3]/Omega2

      ## initially canceling terms
      # EXP
      #1
      dift2 <- dift^2
      #k
      T2 <- tau^2
      #n_exp
      S1 <- dexp2(dtau[,1],Exp0[,1])
      S2 <- dexp2(dtau[,2],Exp0[,2])
      #n_exp
      S12 <- 2*tau[1]*tau[2]*dexp1(fdt,Exp0[,1]*Exp0[,2])
      Sigma[!dt_inf,1,1][EXP] <- (T2[1]*S1 - S12 + T2[2]*S2)/dift2
      Sigma[!dt_inf,2,2][EXP] <- (T2[2]*S1 - S12 + T2[1]*S2)/dift2 * Omega2

      # !EXP
      CROSS <- fdt[!dt_inf][!EXP]*Sinc0*Exp
      OUTER <- Cos0^2*dexp2(fdt[!dt_inf][!EXP],Exp) - CROSS^2
      CROSS <- 2*Cos0*Exp*CROSS
      Sin2 <- Sin0^2

      if(tau[1]>tau[2])
      {
        Sigma[!dt_inf,1,1][!EXP] <- OUTER - Sin2 - CROSS
        Sigma[!dt_inf,2,2][!EXP] <- (OUTER - Sin2 + CROSS) * Omega2
      }
      else
      {
        Sigma[!dt_inf,1,1][!EXP] <- OUTER + Sin2 - CROSS
        Sigma[!dt_inf,2,2][!EXP] <- (OUTER + Sin2 + CROSS) * Omega2
      }

      # back to dt < Inf, (both EXP and !EXP)
      # initially vanishing terms
      c12 <- c012[,2]^2
      Sigma[!dt_inf,1,1] <- Sigma[!dt_inf,1,1] - c12/Omega2
      Sigma[!dt_inf,1,2] <- TT*c12
      Sigma[!dt_inf,2,1] <- TT*c12
      Sigma[!dt_inf,2,2] <- Sigma[!dt_inf,2,2] - c12
    } # end OUF/OUO
  }
  # fix the dimension of the filter
  if(DIM==1) # 1D filter
  { Sigma <- sigma * Sigma }
  else # 2D filter
  {
    if(length(sigma)==1) { sigma <- diag(sigma,DIM) }
    K <- max(1,K)

    # TODO: Replace apply with matrix operations
    Sigma <- apply(Sigma, 1, function(Sigma) {
      Sigma <- outer(Sigma,sigma) # (k,k,d,d)
      Sigma <- aperm(Sigma,c(1,3,2,4)) # (k,k,d,d) -> (k,d,k,d)
      dim(Sigma) <- c(K*DIM,K*DIM)
      return(Sigma)
    })
    dim(Sigma) <- c(n,K*DIM,K*DIM)

    # BM/IOU prior fix
    NAN <- is.nan(Sigma)
    if(any(NAN)) { Sigma[NAN] <- 0 }

    Green <- apply(Green, 1, function(Green) {
      Green <- outer(Green,diag(DIM)) # (k,k,d,d)
      Green <- aperm(Green,c(1,3,2,4)) # (k,d,k,d)
      dim(Green) <- c(K*DIM,K*DIM)
      return(Green)
    })
    dim(Green) <- c(n,K*DIM,K*DIM)
  }
  return(list(Green=Green, Sigma=Sigma))

}