
rm(list=ls())
set.seed(1234)
############################################
#Chunk code:  source
############################################

source("RGEN.R")

############################################
#Chunk code:  libraries and inla options
############################################
inla.setOption(inla.call=NULL)
library(INLA)
library("suddengains")
library(tidyr)
library(mvtnorm)
inla.setOption(pardiso.license = "~/pardiso.lic")
inla.setOption(inla.mode="experimental")



############################################
#Chunk code:  # data simulation
############################################

nsujet=500 # number of indivduals

# Biomarker 1
b1_0=0.2 # intercept
b1_1=-0.1 # slope
b1_2=-0.1 # quadratic slope
b1_3=0.1 # cubic slope
b1_4=0.1 # continuous covariate
b1_5=-0.2 # binary covariate
b1_e=0.01 # residual error

# Biomarker 2
b2_0=0.2 # intercept
b2_1=-0.1 # slope
b2_2=-0.1 # quadratic slope
b2_3=0.1 # cubic slope
b2_4=0.1 # continuous covariate
b2_5=-0.2 # binary covariate
b2_e=0.01 # residual error

gap=0.5
followup=5 # follow-up time

############################################
#Chunk code:  graphs 
############################################
S <- list(p1 ~ p2 + p3, p2 ~ p4 + p5, p3 ~ p6 + p7,
          p4 ~ c1 + c2, p5 ~ c3 + c4,
          p6 ~ c5 + c6, p7 ~ c7 + c8)
SP <- GraphDens(S)
SP_plot = GraphPlot(S, base=0) # base argument
plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)

############################################
#Chunk code:  random effects
############################################
# random effects variance and covariance
c1=1 # random intercept biomarker 1
c2=0.2 # random slope biomarker 1
c3=0.1 # random quadratic slope biomarker 1
c4=0.5 # random cubic slope biomarker 1

c5=3 # random intercept biomarker 2
c6=4 # random slope biomarker 2
c7=0.5 # random quadratic slope biomarker 2
c8=0.2 # random cubic slope biomarker 2

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

c5c7=0.7
c5c8=0.7
c6c7=0.7
c6c8=0.7

# third level
c1c5=0.6
c1c6=0.6
c2c5=0.6
c2c6=0.6
c3c5=0.6
c3c6=0.6
c4c5=0.6
c4c6=0.6

c1c7=0.5
c1c8=0.5
c2c7=0.5
c2c8=0.5
c3c7=0.5
c3c8=0.5
c4c7=0.5
c4c8=0.5

############################################
#Chunk code: matrix std dv and cor
############################################
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

############################################
#Chunk code:  time setting
############################################


mestime=seq(0,followup,gap) # measurement times
timej=rep(mestime, nsujet) # time column
nmesindiv=followup/gap+1 # number of individual measurements

nmesy= nmesindiv*nsujet # number of longi measurements
idY<-as.factor(rep(1:nsujet, each=nmesindiv)) # id

### begin data generation
# random effects generation
MVnorm <- rmvnorm(nsujet, rep(0, 8), Sigma)

b1_int = MVnorm[,1]
b1_intY <- rep(b1_int, each=nmesindiv)
b1_slo = MVnorm[,2]
b1_sloY <- rep(b1_slo, each=nmesindiv)
b1_slo2 = MVnorm[,3]
b1_slo2Y <- rep(b1_slo2, each=nmesindiv)
b1_slo3 = MVnorm[,4]
b1_slo3Y <- rep(b1_slo3, each=nmesindiv)
b2_int = MVnorm[,5]
b2_intY <- rep(b2_int, each=nmesindiv)
b2_slo = MVnorm[,6]
b2_sloY <- rep(b2_slo, each=nmesindiv)
b2_slo2 = MVnorm[,7]
b2_slo2Y <- rep(b2_slo2, each=nmesindiv)
b2_slo3 = MVnorm[,8]
b2_slo3Y <- rep(b2_slo3, each=nmesindiv)


binX=rbinom(nsujet,1, 0.5) # binary covariate
binXY=rep(binX, each=nmesindiv)

ctsX=round(rnorm(nsujet,1, 0.5),2) # continuous covariate
ctsXY=rep(ctsX, each=nmesindiv)


############################################
#Chunk code: linear predictors
############################################

linPred1 <- b1_0+b1_intY+(b1_1+b1_sloY)*timej+
            (b1_2+b1_slo2Y)*timej^2+(b1_3+b1_slo3Y)*timej^3+
            b1_4*ctsXY+b1_5*binXY
linPred2 <- b2_0+b2_intY+(b2_1+b2_sloY)*timej+
            (b2_2+b2_slo2Y)*timej^2+(b2_3+b2_slo3Y)*timej^3+
            b2_4*ctsXY+b2_5*binXY
Y1 <- rnorm(nmesy, linPred1, b1_e)
Y2 <- rnorm(nmesy, linPred2, b2_e)



############################################
#Chunk code:  longitudinal dataset
############################################

id <- as.integer(idY)
longDat <- data.frame(id, timej, binXY, ctsXY, Y1, Y2)
print(head(longDat, 20))
print(str(longDat))
print(summary(longDat))


############################################
#Chunk code:  longitudinal dataset
############################################


data_wide <- spread(longDat, timej,Y2)

plot_sg_trajectories(data =data_wide,
                     id_var = "id",
                     #select_id_list = c("2", "4", "5", "9", "100"),
                     var_list = c("0" ,    "0.5" ,  "1"   ,  "1.5"  , "2" ,    "2.5"  , "3"    , "3.5"   ,"4"  ,   "4.5"  , "5"),
                     show_id = TRUE,
                     id_label_size = 4,
                     label.padding = .2,
                     show_legend = FALSE,
                     colour = "viridis",
                     viridis_option = "D",
                     viridis_begin = 0,
                     viridis_end = .8,
                     connect_missing = FALSE,
                     scale_x_num = TRUE,
                     scale_x_num_start = 1,
                     apaish = TRUE,
                     xlab = "Time",
                     ylab = "Biomarker Y2")



############################################
#Chunk code:  dataset organization 
############################################

NL1 <- dim(longDat)[1]
Y1 <- c(longDat$Y1, rep(NA, NL1))
Y2 <- c(rep(NA, NL1), longDat$Y2)
Yjoint <- list(Y1,Y2)



############################################
#Chunk code:  linear cov dataframe 
############################################


linear.covariate <- data.frame(
  InteY1 = c(rep(1,NL1), rep(NA, NL1)), # intercept Y1
  TIMEY1 = c(longDat$timej,rep(NA, NL1)), # time Y1
  TIME2Y1 = c(longDat$timej^2,rep(NA, NL1)), # time Y1
  TIME3Y1 = c(longDat$timej^3,rep(NA, NL1)), # time Y1
  binXY1 = c(longDat$binX, rep(NA, NL1)), # binX Y1
  ctsXY1 = c(longDat$ctsX, rep(NA, NL1)), # ctsX Y1
  InteY2 = c(rep(NA,NL1), rep(1, NL1)), # intercept Y2
  TIMEY2 = c(rep(NA,NL1), longDat$timej), # time Y2
  TIMEY22 = c(rep(NA,NL1), longDat$timej^2), # time Y2
  TIMEY32 = c(rep(NA,NL1), longDat$timej^3), # time Y2
  binXY2 = c(rep(NA,NL1), longDat$binX), # binX Y2
  ctsXY2 = c(rep(NA,NL1), longDat$ctsX)) # ctsX Y2
random.covariate<-data.frame(idx=c(rep(1, NL1), rep(5, NL1)), # random intercept Y1 and Y2
                             idx2=c(rep(2, NL1), rep(6, NL1)), # random slope Y1 and Y2
                             idx3=c(rep(3, NL1), rep(7, NL1)), # random quadratic slope
                             idx4=c(rep(4, NL1), rep(8, NL1)), # random cubic slope
                             Widx=c(longDat$timej, longDat$timej),
                             Widx2=c(longDat$timej^2, longDat$timej^2),
                             Widx3=c(longDat$timej^3, longDat$timej^3),
                             idxR=c(longDat$id, longDat$id))
jointdf = append(linear.covariate, random.covariate)
joint.data <- append(jointdf, list(Yjoint=Yjoint))


############################################
#Chunk code:   forumla for INLA
############################################

formulaJ= Yjoint ~ -1 + InteY1 + TIMEY1+ TIME2Y1 + TIME3Y1 + binXY1 + ctsXY1  +
  InteY2 + TIMEY2 + TIMEY22 + TIMEY32 +binXY2 + ctsXY2 +
  f(idx, model = CorGraphs, replicate = idxR) + # random intercepts
  f(idx2, Widx, copy="idx", replicate = idxR) + # random slopes
  f(idx3, Widx2, copy="idx", replicate = idxR) + # random quadratic slopes
  f(idx4, Widx3, copy="idx", replicate = idxR) # random cubic slopes



############################################
#Chunk code:   graph and prior
############################################
S <- list(p1 ~ p2 + p3, p2 ~ p4 + p5, p3 ~ p6 + p7,
          p4 ~ c1 + c2, p5 ~ c3 + c4,
          p6 ~ c5 + c6, p7 ~ c7 + c8)
lambda=3

SP_plot = GraphPlot(S, base=0) # base argument
plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)


SP <- GraphDens(S)
Argm <- list(S=S, lambda=lambda, SP=SP, Tdist=Tdist, GraphPrior=GraphPrior, init=-1)
CorGraphs = inla.rgeneric.define(CorGraphs.model, Argm=Argm)



############################################
#Chunk code:   INLA fit
############################################
# inla.setOption(num.threads="7:1")
JMinla <- inla(formulaJ,family = c("gaussian", "gaussian"),
               data=joint.data,verbose = T, safe=T,
               control.compute=list(config=TRUE))



############################################
#Chunk code:   results
############################################
summary(JMinla)
# JMinla <- inla.rerun(JMinla) # rerun to get more accurate estimates
 # save(JMinla, file="RandomEff_Example_Lam7.RData")

print(summary(JMinla))
thet <- c(JMinla$summary.hyperpar$mean[grep("idx", rownames(JMinla$summary.hyperpar))])
ESTcor <- ThetaCor(SP, thet)
colnames(ESTcor) = rownames(ESTcor) = SP$STR[[1]][(SP$NP+1):(SP$NP+SP$NC)]
print(ESTcor) # estimated sd/correlation
print(round(ESTcor, 3)) # estimated sd/correlation
SGM <- cov2cor(Sigma)
diag(SGM) <- sqrt(diag(Sigma))
print(SGM) # true correlation

ESTcov <- ThetaCor(SP, thet, COV=T)
colnames(ESTcov) = rownames(ESTcov) = SP$STR[[1]][(SP$NP+1):(SP$NP+SP$NC)]
print(ESTcov) # estimated variance-covariance
print(Sigma) # true variance-covariance


############################################
#Chunk code:    uncertainty
############################################


NS <- 1000
MATH <- NULL
HSMP <- inla.hyperpar.sample(NS, JMinla)
for(i in 1:NS){
  MATH <-rbind(MATH, c(ThetaCor(SP, HSMP[i, 3:17])))
}
LO_ESTcor <- apply(MATH, 2, function(x) quantile(x, 0.025))
UP_ESTcor <- apply(MATH, 2, function(x) quantile(x, 0.975))
LO_ESTcorM <- matrix(LO_ESTcor, ncol=8, nrow=8)
UP_ESTcorM <- matrix(UP_ESTcor, ncol=8, nrow=8)
round(LO_ESTcorM, 3)
round(UP_ESTcorM, 3)
































