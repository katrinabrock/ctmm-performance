
###############################
# Propagator/Green's function and Two-time correlation from Langevin equation for Kalman filter and simulations
# random CTMM objects need to be run through get.taus() first, to precompute various parameters
v_langevin <- function(dt,CTMM,DIM=1)
{
  K <- CTMM$K
  tau <- CTMM$tau
  sigma <- methods::getDataPart(CTMM$sigma)

  n <- length(dt)
  dtf <- dt<Inf # true dt is finite
  adtf <- any(dtf)
  if(K<=1) # IID-BM-OU
  {
    # IID limit
    Green <- array(0,c(n,1,1))
    Sigma <- array(1,c(n,1,1))

    if(K)
    {
      if(tau[1]==Inf) # BM
      {
        Green[dtf,1,1] <- 1
        # absorbing 1/tau into sigma # VAR -> Diffusion
        Sigma[,1,1] <- 2*dt
      }
      else # (BM,OU,IID]
      {
        if(adtf) dtau <- dt[dtf]/tau
        c0 <- exp(-dtau)
        Green[dtf,1,1] <- c0
        Sigma[dtf,1,1] <- dexp2(dtau,Exp=c0)
      }
    } # >IID
  } # IID-BM-OU
  else if(K==2) # IOU-OUF-OUO
  {
    Omega2 <- CTMM$Omega2
    #IID limit
    Green <- array(0, c(n,2,2))
    Sigma <- array(rep(c(1,0,0,Omega2), each = n), c(n,2,2) )

    f <- CTMM$f.nu[1] # mean(f)
    nu <- CTMM$f.nu[2] # nu || omega
    TT <- CTMM$TfOmega2 # 2 f / Omega^2
    fdt <- f*dt

    if(tau[1]==Inf) # IOU
    {
      dtau <- dt/tau[2]
      Exp <- exp(-dtau)
      DExp <- dexp1(dtau,Exp) # 1-exp(dt/tau[2])
  
      if(adtf){
        Green[dtf,1,1] <- 1
        Green[dtf,1,2] <- tau[2]*DExp[dtf]
        Green[dtf,2,2] <- Exp[dtf]
      }

      # remember that sigma is D=sigma/tau[1]
      DExp2 <- DExp^2 # (1-exp(-dt/tau[2]))^2
      Sigma[,1,1] <- clamp( 2*dt - tau[2]*(2*DExp+DExp2) ,0,Inf) # does this still underflow?
      Sigma[,2,2] <- dexp2(dtau,Exp)/tau[2]

      # dt<Inf
      if(adtf){
        s1221 <- clamp(DExp2[dtf],0, sqrt(Sigma[dtf,1,1]*Sigma[dtf,2,2]) ) # how does this get so far off?
        Sigma[dtf,2,1] <- s1221
        Sigma[dtf,1,2] <- s1221
      }
      # 0 at dt=Inf
    } # END IOU
    else # (IOU,OUF/OUO,IID]
    {
      # function representation choice
      nudt <- nu*dt[dtf]
      # length(EXP) == sum(dtf)
      EXP <- (tau[1]>tau[2] & nudt>0.8813736)
      aEXP <- any(EXP)
      anEXP <- any(!EXP)

      c012 <- array(0,c(sum(dtf), 3))
      if(aEXP){
        # exponential functions
        n_exp <- sum(EXP)
        #nxk
        dtau <- array(dt[dtf][EXP], c(n_exp,K))/array(rep(tau, each = n_exp), c(n_exp,K))
        #1
        dift <- diff(tau)
        #nxk
        Exp0 <- exp(-dtau)
        #nxk
        Exp <- Exp0/dift
        #sum(dtf)
        c012[EXP,1] <- Exp[,2]*tau[2] - Exp[,1]*tau[1]
        c012[EXP,2] <- Exp[,1] - Exp[,2] 
        c012[EXP,3] <- Exp[,2]/tau[2] - Exp[,1]/tau[1]
      } 
      if(anEXP){
        # trigonometric and hyperbolic-trigonometric functions
        Exp <- exp(-fdt[dtf][!EXP])

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

        c012[!EXP,1] <- CosE + fdt[dtf][!EXP]*SincE
        c012[!EXP,2] <- -(Omega2*dt[dtf][!EXP])*SincE
        c012[!EXP,3] <- -Omega2*(CosE - fdt[dtf][!EXP]*SincE)
      }
      # end function representation

      Green[dtf,1,1] <- c012[,1]
      Green[dtf,2,1] <- c012[,2]
      Green[dtf,1,2] <- -c012[,2]/Omega2
      Green[dtf,2,2] <- -c012[,3]/Omega2

      ## initially canceling terms
      if(aEXP){
        #1
        dift2 <- dift^2
        #k
        T2 <- tau^2
        #n_exp
        S1 <- dexp2(dtau[,1],Exp0[,1])
        S2 <- dexp2(dtau[,2],Exp0[,2])
        #n_exp
        S12 <- 2*tau[1]*tau[2]*dexp1(fdt,Exp0[,1]*Exp0[,2])
        Sigma[dtf,1,1][EXP] <- (T2[1]*S1 - S12 + T2[2]*S2)/dift2
        Sigma[dtf,2,2][EXP] <- (T2[2]*S1 - S12 + T2[1]*S2)/dift2 * Omega2
      }
      if(anEXP){
        # !EXP
        CROSS <- fdt[dtf][!EXP]*Sinc0*Exp
        OUTER <- Cos0^2*dexp2(fdt[dtf][!EXP],Exp) - CROSS^2
        CROSS <- 2*Cos0*Exp*CROSS
        Sin2 <- Sin0^2

        if(tau[1]>tau[2])
        {
          Sigma[dtf,1,1][!EXP] <- OUTER - Sin2 - CROSS
          Sigma[dtf,2,2][!EXP] <- (OUTER - Sin2 + CROSS) * Omega2
        }
        else
        {
          Sigma[dtf,1,1][!EXP] <- OUTER + Sin2 - CROSS
          Sigma[dtf,2,2][!EXP] <- (OUTER + Sin2 + CROSS) * Omega2
        }
      }
      # initially vanishing terms
      c12 <- c012[,2]^2
      Sigma[dtf,1,1] <- Sigma[dtf,1,1] - c12/Omega2
      Sigma[dtf,1,2] <- TT*c12
      Sigma[dtf,2,1] <- TT*c12
      Sigma[dtf,2,2] <- Sigma[dtf,2,2] - c12
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


v_Langevin <- function(t,CTMM,DIM=1)
{
  # time-lag information from ctmm.prepare
  dt <- CTMM$dt
  dtl <- CTMM$dtl
  dti <- match(dt, dtl)

  n <- length(dt)
  tau <- CTMM$tau
  K <- max(1,length(tau))  # dimension of hidden state per spatial dimension

  # propagation information
  Green <- array(diag(K*DIM),c(K*DIM,K*DIM,n))
  Green <- aperm(Green,c(3,1,2)) # [n,K*DIM,K*DIM]
  Sigma <- array(0,c(n,K*DIM,K*DIM))

  dynamics <- CTMM$dynamics
  # default stationary process
  if(is.null(dynamics) || dynamics==FALSE || dynamics=="stationary")
  {
    nl <- length(dtl)
    Greenl <- array(0,c(nl,K*DIM,K*DIM))
    Sigmal <- array(0,c(nl,K*DIM,K*DIM))

    LANGEVIN <- v_langevin(dt=dtl,CTMM=CTMM,DIM=DIM)

    for(i in 1:n)
    {
      Green[i,,] <- LANGEVIN$Green[dti[i],,]
      Sigma[i,,] <- LANGEVIN$Sigma[dti[i],,]
    }
  }
  else if(dynamics=="change.point")
  {
    CP <- CTMM[[dynamics]] # change points
    CS <- get.states(CTMM) # states
    j <- 1
    for(i in 1:n)
    {
      while(t[i]>CP$stop[j]) # does this time increment cross a change point?
      {
        DT <- CP$stop[j] - max(t[i-1],CP$start[j]) # time until change
        LANGEVIN <- v_langevin(dt=DT,CTMM=CTMM[[CS[j]]],DIM=DIM)
        Green[i,,] <- LANGEVIN$Green %*% Green[i,,]
        Sigma[i,,] <- (LANGEVIN$Green %*% Sigma[i,,] %*% t(LANGEVIN$Green)) + LANGEVIN$Sigma
        dt[i] <- dt[i] - DT
        j <- j + 1
      }

      if(i>1 && CP$start[j]>t[i-1]) # did we cross a change point?
      {
        LANGEVIN <- v_langevin(dt=dt[i],CTMM=CTMM[[CS[j]]],DIM=DIM)
        Green[i,,] <- LANGEVIN$Green %*% Green[i,,]
        Sigma[i,,] <- (LANGEVIN$Green %*% Sigma[i,,] %*% t(LANGEVIN$Green)) + LANGEVIN$Sigma
      }
      else # we did not cross a change point
      {
        # do we need a fresh calculation?
        if(i==1 || dt[i] != dt[i-1]) { LANGEVIN <- v_langevin(dt=dt[i],CTMM=CTMM,DIM=DIM) }

        Green[i,,] <- LANGEVIN$Green # (K*DIM,K*DIM)
        Sigma[i,,] <- LANGEVIN$Sigma # (K*DIM,K*DIM)
      }
    } # end for i in 1:n
  } # end change.point dynamics

  R <- list(Green=Green,Sigma=Sigma)
  return(R)
}