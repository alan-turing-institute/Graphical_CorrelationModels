
library(INLA)
library(mvtnorm)
library(corGraphs)
inla.setOption(num.threads="12:1")
set.seed(1)
# data simulation
nsujet=30 # number of indivduals

# Biomarker 1
b1_0=0.2 # intercept
b1_1=-0.1 # slope
b1_4=0.1 # continuous covariate
b1_5=-0.2 # binary covariate
b1_e=0.01 # residual error

# Biomarker 2
b2_0=0.2 # intercept
b2_1=-0.1 # slope
b2_4=0.1 # continuous covariate
b2_5=-0.2 # binary covariate
b2_e=0.01 # residual error

# Biomarker 3
b3_0=0.2 # intercept
b3_1=-0.1 # slope
b3_4=0.1 # continuous covariate
b3_5=-0.2 # binary covariate
b3_e=0.01 # residual error

# Biomarker 4
b4_0=0.2 # intercept
b4_1=-0.1 # slope
b4_4=0.1 # continuous covariate
b4_5=-0.2 # binary covariate
b4_e=0.01 # residual error

gap=1
followup=3 # follow-up time

# random effects variance and covariance
c1=1 # random intercept biomarker 1
c2=0.5 # random slope biomarker 1
c3=1 # random quadratic slope biomarker 1
c4=0.5 # random cubic slope biomarker 1

c5=1 # random intercept biomarker 2
c6=0.5 # random slope biomarker 2
c7=1 # random quadratic slope biomarker 2
c8=0.5 # random cubic slope biomarker 2

# highest correlation
c1c2=0.9
c3c4=0.9
c5c6=0.9
c7c8=0.9

# second level
c1c3=0.8
c1c4=0.8
c2c3=0.8
c2c4=0.8

c5c7=0.8
c5c8=0.8
c6c7=0.8
c6c8=0.8

# third level
c1c5=0.6
c1c6=0.6
c2c5=0.6
c2c6=0.6
c3c5=0.6
c3c6=0.6
c4c5=0.6
c4c6=0.6

c1c7=0.6
c1c8=0.6
c2c7=0.6
c2c8=0.6
c3c7=0.6
c3c8=0.6
c4c7=0.6
c4c8=0.6

SDCOR <- matrix(c(c1, c1c2, c1c3, c1c4, c1c5, c1c6, c1c7, c1c8,
                  c1c2, c2, c2c3, c2c4, c2c5, c2c6, c2c7, c2c8,
                  c1c3, c2c3, c3, c3c4, c3c5, c3c6, c3c7, c3c8,
                  c1c4, c2c4, c3c4, c4, c4c5, c4c6, c4c7, c4c8,
                  c1c5, c2c5, c3c5, c4c5, c5, c5c6, c5c7, c5c8,
                  c1c6, c2c6, c3c6, c4c6, c5c6, c6, c6c7, c6c8,
                  c1c7, c2c7, c3c7, c4c7, c5c7, c6c7, c7, c7c8,
                  c1c8, c2c8, c3c8, c4c8, c5c8, c6c8, c7c8, c8), ncol=8)
R <- SDCOR
diag(R) <- 1
S <- diag(SDCOR)
Sigma <- diag(S) %*% R %*% diag(S)

# Sigma <- CatDyn::cor2cov(SDCOR)

mestime=seq(0,followup,gap) # measurement times
timej=rep(mestime, nsujet) # time column
nmesindiv=followup/gap+1 # number of individual measurements

nmesy= nmesindiv*nsujet # number of longi measurements
idY<-as.factor(rep(1:nsujet, each=nmesindiv)) # id

### begin data generation
nsim=100
res_sim <- NULL
CT <- NULL
DIC <- NULL
WAIC <- NULL
for(simi in 1:nsim){
  set.seed(simi)
  # random effects generation
  MVnorm <- rmvnorm(nsujet, rep(0, 8), Sigma)

  b1_int = MVnorm[,1]
  b1_intY <- rep(b1_int, each=nmesindiv)
  b1_slo = MVnorm[,2]
  b1_sloY <- rep(b1_slo, each=nmesindiv)
  b2_int = MVnorm[,3]
  b2_intY <- rep(b2_int, each=nmesindiv)
  b2_slo = MVnorm[,4]
  b2_sloY <- rep(b2_slo, each=nmesindiv)
  b3_int = MVnorm[,5]
  b3_intY <- rep(b3_int, each=nmesindiv)
  b3_slo = MVnorm[,6]
  b3_sloY <- rep(b3_slo, each=nmesindiv)
  b4_int = MVnorm[,7]
  b4_intY <- rep(b4_int, each=nmesindiv)
  b4_slo = MVnorm[,8]
  b4_sloY <- rep(b4_slo, each=nmesindiv)


  binX=rbinom(nsujet,1, 0.5) # binary covariate
  binXY=rep(binX, each=nmesindiv)

  ctsX=round(rnorm(nsujet,1, 0.5),2) # continuous covariate
  ctsXY=rep(ctsX, each=nmesindiv)

  # linear predictors
  linPred1 <- b1_0+b1_intY+(b1_1+b1_sloY)*timej+
    b1_4*ctsXY+b1_5*binXY
  linPred2 <- b2_0+b2_intY+(b2_1+b2_sloY)*timej+
    b2_4*ctsXY+b2_5*binXY
  linPred3 <- b3_0+b3_intY+(b3_1+b3_sloY)*timej+
    b3_4*ctsXY+b3_5*binXY
  linPred4 <- b4_0+b4_intY+(b4_1+b4_sloY)*timej+
    b4_4*ctsXY+b4_5*binXY
  Y1 <- rnorm(nmesy, linPred1, b1_e)
  Y2 <- rnorm(nmesy, linPred2, b2_e)
  Y3 <- rnorm(nmesy, linPred3, b3_e)
  Y4 <- rnorm(nmesy, linPred4, b4_e)

  # longitudinal dataset
  id <- as.integer(idY)
  longDat <- data.frame(id, timej, binXY, ctsXY, Y1, Y2, Y3, Y4)

  library(INLAjoint)
  A <- Sys.time()
  verif <- joint(formLong = list(Y1 ~ timej + binXY + ctsXY + (1+timej|id),
                                 Y2 ~ timej + binXY + ctsXY + (1+timej|id),
                                 Y3 ~ timej + binXY + ctsXY + (1+timej|id),
                                 Y4 ~ timej + binXY + ctsXY + (1+timej|id)),
                 corLong=T, id="id", timeVar="timej", dataLong = longDat,
                 family=c("gaussian", "gaussian", "gaussian", "gaussian"))
  B <- Sys.time()
  CT2 <- difftime(B, A, units="secs")
  sverif <- summary(verif, sdcor=T)
  res_sim <- rbind(res_sim,
                   c(round(sverif$ReffList[[1]][, 4][c(1:8, 9, 22, 31, 36, 10, 32, 12, 14, 23, 25)], 3),
                     round(verif$summary.fixed$`0.5quant`, 3)))

  CT <- c(CT, CT2)
  DIC <- c(DIC, verif$dic$dic)
  WAIC <- c(WAIC, verif$waic$waic)
}
print(apply(res_sim, 2, mean))
print(apply(res_sim, 2, function(x) quantile(x, c(0.025))))
print(apply(res_sim, 2, function(x) quantile(x, c(0.975))))

save(res_sim, CT, DIC, WAIC, file="./4Y_30_FULL.RData")


























