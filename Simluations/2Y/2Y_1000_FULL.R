
library(INLA)
library(mvtnorm)
library(corGraphs)
inla.setOption(num.threads="12:1")
set.seed(1)
# data simulation
nsujet=1000 # number of indivduals

# Biomarker 1
b1_0=0.2 # intercept
b1_1=-0.1 # slope
b1_e=0.01 # residual error

# Biomarker 2
b2_0=0.2 # intercept
b2_1=-0.1 # slope
b2_e=0.01 # residual error

gap=0.5
followup=5 # follow-up time

S <- list(p1 ~ p2 + p3, p2 ~ c1 + c2, p3 ~ c3 + c4)
SP <- GraphDens(S)
SP_plot = GraphPlot(S, base=0) # base argument
#plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)


# random effects variance and covariance
c1=1 # random intercept biomarker 1
c2=0.2 # random slope biomarker 1
c3=0.1 # random quadratic slope biomarker 1
c4=0.5 # random cubic slope biomarker 1

# highest correlation
c1c2=0.9
c3c4=0.9

# second level
c1c3=0.8
c1c4=0.8
c2c3=0.8
c2c4=0.8

SDCOR <- matrix(c(c1, c1c2, c1c3, c1c4,
                  c1c2, c2, c2c3, c2c4,
                  c1c3, c2c3, c3, c3c4,
                  c1c4, c2c4, c3c4, c4), ncol=4)
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
  MVnorm <- rmvnorm(nsujet, rep(0, 4), Sigma)

  b1_int = MVnorm[,1]
  b1_intY <- rep(b1_int, each=nmesindiv)
  b1_slo = MVnorm[,2]
  b1_sloY <- rep(b1_slo, each=nmesindiv)
  b2_int = MVnorm[,3]
  b2_intY <- rep(b2_int, each=nmesindiv)
  b2_slo = MVnorm[,4]
  b2_sloY <- rep(b2_slo, each=nmesindiv)

  # linear predictors
  linPred1 <- b1_0+b1_intY+(b1_1+b1_sloY)*timej
  linPred2 <- b2_0+b2_intY+(b2_1+b2_sloY)*timej
  Y1 <- rnorm(nmesy, linPred1, b1_e)
  Y2 <- rnorm(nmesy, linPred2, b2_e)

  # longitudinal dataset
  id <- as.integer(idY)
  longDat <- data.frame(id, timej, Y1, Y2)

  library(INLAjoint)
  A <- Sys.time()
  verif <- joint(formLong = list(Y1 ~ timej + (1+timej|id),
                                 Y2 ~ timej + (1+timej|id)),
                 corLong=T, id="id", timeVar="timej", dataLong = longDat,
                 family=c("gaussian", "gaussian"))
  B <- Sys.time()
  CT2 <- difftime(B, A, units="secs")
  sverif <- summary(verif, sdcor=T)
  res_sim <- rbind(res_sim,
                   c(round(sverif$ReffList[[1]][, 4][c(1:4, 5, 10, 6)], 3),
                     round(verif$summary.fixed$`0.5quant`, 3)))

  CT <- c(CT, CT2)
  DIC <- c(DIC, verif$dic$dic)
  WAIC <- c(WAIC, verif$waic$waic)
}
print(apply(res_sim, 2, mean))
print(apply(res_sim, 2, function(x) quantile(x, c(0.025))))
print(apply(res_sim, 2, function(x) quantile(x, c(0.975))))

save(res_sim, CT, DIC, WAIC, file="./2Y_1000_FULL.RData")


























