
library(INLA)
library(mvtnorm)
library(corGraphs)
inla.setOption(num.threads="12:1")
set.seed(1234)
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


  binX=rbinom(nsujet,1, 0.5) # binary covariate
  binXY=rep(binX, each=nmesindiv)

  ctsX=round(rnorm(nsujet,1, 0.5),2) # continuous covariate
  ctsXY=rep(ctsX, each=nmesindiv)

  # linear predictors
  linPred1 <- b1_0+b1_intY+(b1_1+b1_sloY)*timej
  linPred2 <- b2_0+b2_intY+(b2_1+b2_sloY)*timej
  Y1 <- rnorm(nmesy, linPred1, b1_e)
  Y2 <- rnorm(nmesy, linPred2, b2_e)

  # longitudinal dataset
  id <- as.integer(idY)
  longDat <- data.frame(id, timej, binXY, ctsXY, Y1, Y2)


  # INLA
  NL1 <- dim(longDat)[1]
  Y1 <- c(longDat$Y1, rep(NA, NL1))
  Y2 <- c(rep(NA, NL1), longDat$Y2)
  Yjoint <- list(Y1,Y2)

  linear.covariate <- data.frame(
    InteY1 = c(rep(1,NL1), rep(NA, NL1)), # intercept Y1
    TIMEY1 = c(longDat$timej,rep(NA, NL1)), # time Y1
    InteY2 = c(rep(NA,NL1), rep(1, NL1)), # intercept Y2
    TIMEY2 = c(rep(NA,NL1), longDat$timej)) # time Y2
  random.covariate<-data.frame(idx=c(rep(1, NL1), rep(3, NL1)), # random intercept Y1 and Y2
                               idx2=c(rep(2, NL1), rep(4, NL1)), # random slope Y1 and Y2
                               Widx=c(longDat$timej, longDat$timej),
                               idxR=c(longDat$id, longDat$id))
  jointdf = append(linear.covariate, random.covariate)
  joint.data <- append(jointdf, list(Yjoint=Yjoint))


  S <- list(p1 ~ p2 + p3, p2 ~ c1 + c2, p3 ~ c3 + c4)
  SP <- GraphDens(S)

  gmodel <- dcg_model(
    dcg = S,
    sigma.prior.reference = rep(5, SP$NC),
    sigma.prior.probability = rep(0.2, SP$NC),
    lambda = 3,
    iprior = 3,
    useINLAprecomp = FALSE,
    debug =  0
  )

  formulaJ= Yjoint ~ -1 + InteY1 + TIMEY1 +
    InteY2 + TIMEY2 +
    f(idx, model = gmodel, replicate = idxR)+#, vb.correct=FALSE) + # random intercepts
    f(idx2, Widx, copy="idx", replicate = idxR)#, vb.correct=FALSE) # random slopes




  # SP_plot = GraphPlot(S, base=0) # base argument
  # plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)
  #
  #
  # Argm <- list(S=S, lambda=lambda, SP=SP, Tdist=Tdist, GraphPrior=GraphPrior, init=-1)
  # CorGraphs = inla.rgeneric.define(CorGraphs.model, Argm=Argm)

  JMinla <- inla(formulaJ,family = c("gaussian", "gaussian"),
                 data=joint.data,verbose = F, safe=T,
                 control.compute=list(config=TRUE, dic=T, waic=T),
                 control.inla=list(h=0.001, tolerance=0.001))

  thet <- c(JMinla$summary.hyperpar$`0.5quant`[grep("idx", rownames(JMinla$summary.hyperpar))])
  ESTcor <- ThetaCor(SP, thet)
  colnames(ESTcor) = rownames(ESTcor) = SP$STR[[1]][(SP$NP+1):(SP$NP+SP$NC)]
  res_sim <- rbind(res_sim,
                   c(round(exp(thet[1:4]), 3),
                     round(ESTcor, 3)[upper.tri(round(ESTcor, 3))][c(1, 6, 2)],
                     round(JMinla$summary.fixed$`0.5quant`, 3)))
  CT <- c(CT, JMinla$cpu.used[4])
  DIC <- c(DIC, JMinla$dic$dic)
  WAIC <- c(WAIC, JMinla$waic$waic)
}
print(apply(res_sim, 2, mean))
print(apply(res_sim, 2, function(x) quantile(x, c(0.025))))
print(apply(res_sim, 2, function(x) quantile(x, c(0.975))))
save(res_sim, CT, DIC, WAIC, file="./2Y_1000_L3.RData")



# uncertainty
# NS <- 1000
# MATH <- NULL
# HSMP <- inla.hyperpar.sample(NS, JMinla)
# for(i in 1:NS){
#   MATH <-rbind(MATH, c(ThetaCor(SP, HSMP[i, 3:9])))
# }
# LO_ESTcor <- apply(MATH, 2, function(x) quantile(x, 0.025))
# UP_ESTcor <- apply(MATH, 2, function(x) quantile(x, 0.975))
# LO_ESTcorM <- matrix(LO_ESTcor, ncol=4, nrow=4)
# UP_ESTcorM <- matrix(UP_ESTcor, ncol=4, nrow=4)
# round(LO_ESTcorM, 3)
# round(UP_ESTcorM, 3)

































