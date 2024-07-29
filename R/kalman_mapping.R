iid_data <- readRDS('/home/TOP/kbrock/ctmm/kalman_in_out_iid.rds')
names(iid_data)<- c('z', 'u', 't', 'dt', 'CTMM', 'error', 'DIM', 'KALMAN')

str(iid_data[['z']])
str(iid_data[['u']]) # mean vector u
str(iid_data[['t']]) # timesteps
str(iid_data[['dt']]) # time diffs
str(iid_data[['CTMM']])
str(iid_data[['error']])
str(iid_data[['DIM']])

str(iid_data[['KALMAN']])

data2 <- iid_data

data2[['KALMAN']] <- NULL

data2[['KALMAN']] <- do.call(ctmm:::kalman, data2)

testthat::expect_equal(data2[['KALMAN']], iid_data[['KALMAN']])

fkf - a0 -  A vector giving the initial value/estimation of the state variable.
kf - B0 for the initial guess of the unobserved components
ctmm - u - ? 

fkf - P0 - A matrix giving the variance of a0.
kf - P0 for the initial guess of the unobserved components covariance matrix
ctmm - sFor[1,,] <- Sigma[1,,]

fkf - dt -  A matrix giving the intercept of the transition equation
kf - Dm -  for the constant in the state equation
ctmm - dt <- c(Inf,diff(t)) # time lags

fkf - ct - A matrix giving the intercept of the measurement equation
kf - Am -  for the constant in the observation equation

fkf - Tt - An array giving the factor of the transition equation
kf - Fm - for the state equation transition matrix

fkf - Zt - An array giving the factor of the measurement equation 
kf -  Hm - for the observation matrix in the observation equation
ctmm - P - observable state projection operator: t(P) %*% FULL -> OB

fkf - HHt - An array giving the variance of the innovations of the transition equation 
kf - Qm - for the covariance matrix of the errors in the state equation,

fkf - GGt - An array giving the variance of the disturbances of the measurement equation
kf - Rm - for the covariance matrix of the errors in the observation equation

fkf - yt -  A matrix containing the observations. “NA”-values are allowed
kf - yt - N x T matrix of data
ctmm - z - data z

ctmm
zFor = forcast / estimate
sFor = error in est








