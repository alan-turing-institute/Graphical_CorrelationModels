
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

  #### Plot
  # library("suddengains")
  # library(tidyr)
  # data_wide <- spread(longDat, timej,Y2)
  #
  # plot_sg_trajectories(data =data_wide,
  #                      id_var = "id",
  #                      #select_id_list = c("2", "4", "5", "9", "100"),
  #                      var_list = c("0" ,    "0.5" ,  "1"   ,  "1.5"  , "2" ,    "2.5"  , "3"    , "3.5"   ,"4"  ,   "4.5"  , "5"),
  #                      show_id = TRUE,
  #                      id_label_size = 4,
  #                      label.padding = .2,
  #                      show_legend = FALSE,
  #                      colour = "viridis",
  #                      viridis_option = "D",
  #                      viridis_begin = 0,
  #                      viridis_end = .8,
  #                      connect_missing = FALSE,
  #                      scale_x_num = TRUE,
  #                      scale_x_num_start = 1,
  #                      apaish = TRUE,
  #                      xlab = "Time",
  #                      ylab = "Biomarker Y2")


  # INLA
  NL1 <- dim(longDat)[1]
  Y1 <- c(longDat$Y1, rep(NA, NL1*3))
  Y2 <- c(rep(NA, NL1), longDat$Y2, rep(NA, NL1*2))
  Y3 <- c(rep(NA, NL1*2), longDat$Y3, rep(NA, NL1))
  Y4 <- c(rep(NA, NL1*3), longDat$Y4)
  Yjoint <- list(Y1, Y2, Y3, Y4)

  linear.covariate <- data.frame(
    InteY1 = c(rep(1,NL1), rep(NA, NL1*3)), # intercept Y1
    TIMEY1 = c(longDat$timej,rep(NA, NL1*3)), # time Y1
    binXY1 = c(longDat$binX, rep(NA, NL1*3)), # binX Y1
    ctsXY1 = c(longDat$ctsX, rep(NA, NL1*3)), # ctsX Y1
    InteY2 = c(rep(NA,NL1), rep(1, NL1), rep(NA, NL1*2)), # intercept Y2
    TIMEY2 = c(rep(NA,NL1), longDat$timej, rep(NA, NL1*2)), # time Y2
    binXY2 = c(rep(NA,NL1), longDat$binX, rep(NA, NL1*2)), # binX Y2
    ctsXY2 = c(rep(NA,NL1), longDat$ctsX, rep(NA, NL1*2)), # ctsX Y2
    InteY3 = c(rep(NA, NL1*2), rep(1,NL1), rep(NA, NL1)), # intercept Y3
    TIMEY3 = c(rep(NA, NL1*2), longDat$timej,rep(NA, NL1)), # time Y3
    binXY3 = c(rep(NA, NL1*2), longDat$binX, rep(NA, NL1)), # binX Y3
    ctsXY3 = c(rep(NA, NL1*2), longDat$ctsX, rep(NA, NL1)), # ctsX Y3
    InteY4 = c(rep(NA, NL1*3), rep(1,NL1)), # intercept Y4
    TIMEY4 = c(rep(NA, NL1*3), longDat$timej), # time Y4
    binXY4 = c(rep(NA, NL1*3), longDat$binX), # binX Y4
    ctsXY4 = c(rep(NA, NL1*3), longDat$ctsX)) # ctsX Y4
  random.covariate<-data.frame(idx=c(rep(1, NL1), rep(3, NL1), rep(5, NL1), rep(7, NL1)), # random intercepts
                               idx2=c(rep(2, NL1), rep(4, NL1), rep(6, NL1), rep(8, NL1)), # random slopes
                               Widx=c(longDat$timej, longDat$timej, longDat$timej, longDat$timej),
                               idxR=c(longDat$id, longDat$id, longDat$id, longDat$id))
  jointdf = append(linear.covariate, random.covariate)
  joint.data <- append(jointdf, list(Yjoint=Yjoint))

  S <- list(p1 ~ p2 + p3, p2 ~ p4 + p5, p3 ~ p6 + p7,
            p4 ~ c1 + c2, p5 ~ c3 + c4,
            p6 ~ c5 + c6, p7 ~ c7 + c8)
  SP <- GraphDens(S)
  # SP_plot = GraphPlot(S, base=0) # base argument
  # plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)

  gmodel <- dcg_model(
    dcg = S,
    sigma.prior.reference = rep(5, SP$NC),
    sigma.prior.probability = rep(0.2, SP$NC),
    lambda = 1,
    iprior = 3,
    useINLAprecomp = FALSE,
    debug =  0
  )

  formulaJ= Yjoint ~ -1 + InteY1 + TIMEY1+ binXY1 + ctsXY1  +
    InteY2 + TIMEY2 + binXY2 + ctsXY2 +
    InteY3 + TIMEY3 + binXY3 + ctsXY3 +
    InteY4 + TIMEY4 + binXY4 + ctsXY4 +
    f(idx, model = gmodel, replicate = idxR) + # random intercepts
    f(idx2, Widx, copy="idx", replicate = idxR) # random slopes

  # inla.setOption(num.threads="7:1")
  JMinla <- inla(formulaJ,family = c("gaussian", "gaussian", "gaussian", "gaussian"),
                 data=joint.data,verbose = F, safe=T,
                 control.compute=list(config=TRUE, dic=T, waic=T),
                 control.inla=list(h=0.001, tolerance=0.001))
  # summary(JMinla)
  # JMinla <- inla.rerun(JMinla) # rerun to get more accurate estimates
  # save(JMinla, file="RandomEff_Example_Lam7.RData")

  # print(summary(JMinla))

  thet <- c(JMinla$summary.hyperpar$`0.5quant`[grep("idx", rownames(JMinla$summary.hyperpar))])
  ESTcor <- ThetaCor(SP, thet)
  colnames(ESTcor) = rownames(ESTcor) = SP$STR[[1]][(SP$NP+1):(SP$NP+SP$NC)]

  res_sim <- rbind(res_sim,
  c(round(exp(thet[1:8]), 3),
  round(ESTcor, 3)[upper.tri(round(ESTcor, 3))][c(1 ,6, 15, 28,
                                                  2, 20, 7, 16,
                                                  9, 18)],
  round(JMinla$summary.fixed$`0.5quant`, 3)))
  CT <- c(CT, JMinla$cpu.used[4])
  DIC <- c(DIC, JMinla$dic$dic)
  WAIC <- c(WAIC, JMinla$waic$waic)
}
print(apply(res_sim, 2, mean))
print(apply(res_sim, 2, function(x) quantile(x, c(0.025))))
print(apply(res_sim, 2, function(x) quantile(x, c(0.975))))
save(res_sim, DIC, WAIC, CT, file="./4Y_30_L1.RData")


# uncertainty
# NS <- 1000
# MATH <- NULL
# HSMP <- inla.hyperpar.sample(NS, JMinla)
# for(i in 1:NS){
#   MATH <-rbind(MATH, c(ThetaCor(SP, HSMP[i, 3:17])))
# }
# LO_ESTcor <- apply(MATH, 2, function(x) quantile(x, 0.025))
# UP_ESTcor <- apply(MATH, 2, function(x) quantile(x, 0.975))
# LO_ESTcorM <- matrix(LO_ESTcor, ncol=8, nrow=8)
# UP_ESTcorM <- matrix(UP_ESTcor, ncol=8, nrow=8)
# round(LO_ESTcorM, 3)
# round(UP_ESTcorM, 3)

# library(INLAjoint)
# verif <- joint(formLong = list(Y1 ~ timej + binXY + ctsXY + (1+timej|id),
#                                Y2 ~ timej + binXY + ctsXY + (1+timej|id),
#                                Y3 ~ timej + binXY + ctsXY + (1+timej|id),
#                                Y4 ~ timej + binXY + ctsXY + (1+timej|id)),
#                corLong=T, id="id", timeVar="timej", dataLong = longDat,
#                family=c("gaussian", "gaussian", "gaussian", "gaussian"))
# print(summary(verif, sdcor=T))
#



























