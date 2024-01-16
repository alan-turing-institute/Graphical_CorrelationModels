rm(list=ls())
set.seed(1234)
#############################
# Chunk Code: settings
#############################


#inla.setOption(inla.call=NULL)
library(INLA)
#inla.setOption(pardiso.license = "~/pardiso.lic")
#inla.setOption(inla.mode="experimental")
set.seed(1234)

#############################
# Chunk Code: libraries
#############################

library(mvtnorm)
library(tidyverse)
library(dplyr)
library(GGally)
library(xtable)
#############################
# Chunk Code: RGEN source
#############################

source('RGEN.R')

#############################
# Chunk Code: generate data
#############################

mean1=c(0,0,0,0)



R <- matrix(c(56, 50,3,2,
              50, 62, 5,6,
              3, 5,27,13,
              2,6,13,24), 
            nrow = 4, ncol = 4, byrow = TRUE)


mydf=as.data.frame(MASS:::mvrnorm(n=215, mu=mean1, Sigma=R))


#############################
# Chunk Code: correlations
#############################
mydf%>%cor()
apply (mydf,2, sd)


mydf%>%ggcorr(label=TRUE)
names(mydf)=c('Horvath', 'Hannum','PhenoAge', 'GrimAge')
#############################
# Chunk Code: graph
#############################

S <- list(p1~ p2 + p3 , p2~c1+c2, p3~c3+c4)
# S <- list(p1~ p2 + c3+c4, p2~ c1+c2)
lambda=7
SP <- GraphDens(S)
SP_plot = GraphPlot(S, base=0)
plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)

#############################
# Chunk Code: data organization
#############################
n=nrow(mydf)
idx=c(rep(1,n), rep(2,n), rep(3,n), rep(4,n))
idxR=rep(1:n, 4)

##########################################
# Chunk Code: Outcomes
#############################

Y <- c(mydf$Horvath, mydf$Hannum, mydf$PhenoAge, mydf$GrimAge)
#############################
# Chunk Code: parameters settings
#############################
lambda=5
Argm <- list(S=S, lambda=lambda, SP=SP, Tdist=Tdist, GraphPrior=GraphPrior, init=-1)
CorGraphs = inla.rgeneric.define(CorGraphs.model, Argm=Argm)

#############################
# Chunk Code: data in a list
#############################
mydata=list(idx=idx, idxR=idxR, Y=Y)

#############################
# Chunk Code: fit
#############################
fit <- inla(Y ~ -1+ f(idx, model = CorGraphs, replicate = idxR),
            family ="gaussian",
              data=mydata, control.family = list(hyper =  list(prec = list(initial = 10,fixed = TRUE))),
              verbose=T)
summary(fit)

#############################
# Chunk Code: results
#############################

thet <- c(fit$summary.hyperpar$mean)
ESTcor <- ThetaCor(SP, thet)
colnames(ESTcor) = rownames(ESTcor) = SP$STR[[1]][(SP$NP+1):(SP$NP+SP$NC)]
print(ESTcor) # cov mat


# uncertainty
NS <- 1000
MATH <- NULL
HSMP <- inla.hyperpar.sample(NS, fit)
for(i in 1:NS){
  MATH <-rbind(MATH, c(ThetaCor(SP, HSMP[i,])))
}
LO_ESTcor <- apply(MATH, 2, function(x) quantile(x, 0.025))
UP_ESTcor <- apply(MATH, 2, function(x) quantile(x, 0.975))
LO_ESTcorM <- matrix(LO_ESTcor, ncol=4, nrow=4)
UP_ESTcorM <- matrix(UP_ESTcor, ncol=4, nrow=4)
round(LO_ESTcorM, 3)
round(UP_ESTcorM, 3)













