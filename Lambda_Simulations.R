############################################
#Chunk code:  library 
############################################

library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(grid)

############################################
#Chunk code:  PC prior functions 
############################################

source("./RGEN.R")

############################################
#Chunk code:   one parent 2 children
############################################
SDsim <- exp(seq(-5,log(5), len=30))
S <- list(p1 ~ c1 + c2)
GR <- GraphDens(S)
lDi <- NULL
for(tt in 1:length(SDsim)){
  lDi <- c(lDi, Tdist(GR$JD[[1]], SD1 = SDsim[tt], STR=GR$STR[[1]] , NC = GR$NC))
}
plot(SDsim, lDi, type="o", log="xy")
MapF <- splinefun(lDi, SDsim)
#curve(MapF, from=-2, to=2)
Dsim <- MapF(1)

lambda <- c(1,3, 5, 7, 9, 11)# c(1,4,5,6,7,10)
#par(mfrow=c(2,3), mar=c(2,2,1,1))
rst=list()
rst <- vector("list", length=length(lambda))
for(lam in 1:length(lambda)){
  # generate distance as exp function of lambda
  DistSim <- rexp(n = 300, rate = lambda[lam])
  # SDsim from distance
  ThetasSim <- MapF(DistSim)
  # cor from SDsim
  rst[[lam]] <- sapply(ThetasSim, function(x) ThetaCor(GR, c(10,10,log(x)))[1,2])
}

############################################
#Chunk code:    plot  
############################################
df= as.data.frame(do.call(c, rst))
names(df)='value'
df$q_values=rep(lambda, each=300)
df$q_values <- factor(df$q_values, 
                        labels=c('lambda==1','lambda==3', 'lambda==5', 'lambda==7', 'lambda==9','lambda==11'))




p6=ggplot(df, aes(x = value))+ 
  geom_histogram(bins=50)+
  facet_wrap(.~q_values,labeller = label_parsed,ncol=3, scales='free')+
  theme_bw()+
  theme(legend.position='bottom', legend.title = element_blank(), panel.spacing = unit(0.5, "lines"),
        axis.title=element_blank())

p6


ggsave(p6, file='lambda_1parent_10x3.pdf',dpi=300)



  
  

############################################
#Chunk code:   (from 1 parent to two parent), then the prior is conditional on the value of the parameter for the first latent:
############################################


SDsim0 <- c(0.1, 1, 2)
lambda <- c(3, 5, 7)
hst <- vector("list", length=length(lambda)*length(SDsim0))
it=1


for(lam in 1:length(lambda)){
  for(tet in 1:length(SDsim0)){
    SDsim <- seq(0.01,5, len=50)
    S <- list(p1 ~ p2 + c1, p2 ~ c2 + c3)
    GR <- GraphDens(S)
    lDi <- NULL
    for(tt in 1:length(SDsim)){
      lDi <- c(lDi, Tdist(GR$JD[[1]], GR$JD[[2]], SD1 = SDsim[tt], SD0 = SDsim0[tet], STR=GR$STR[[1]], NC = GR$NC))
    }
    #plot(SDsim, lDi, type= o", main=paste0("SD (base)=", SDsim0[tet]), xlab="SD (new latent)")
    MapF <- splinefun(lDi, SDsim)
    #curve(MapF, from=0, to=2)
    Dsim <- MapF(1)
    
    # generate distance as exp function of lambda
    DistSim <- rexp(n = 200, rate = lambda[lam])
    # SDsim from distance
    ThetasSim <- MapF(DistSim)
    # cor from SDsim
    rst<- sapply(ThetasSim, function(x) ThetaCor(GR, c(10,10, 10,log(x), log(SDsim0[tet])))[1,2])
    hst[[it]] <- rst
    #hist(rst, main= paste0("lambda=", lambda[lam], " / SD (base)=", SDsim0[tet]),breaks=50, xlab="correlation", xlim=c(0,1))
    it <- it+1
    print(it)
  }
}


############################################
#Chunk code:   plot 
############################################


df= as.data.frame(do.call(c, hst))
names(df)='value'
df$q_values=rep(1:9, each=200)

df$q_values <- factor(df$q_values, 
                      labels=c('lambda==3~sigma[base]==0.1', 'lambda==3~sigma[base]==1',  'lambda==3~sigma[base]==2', 
                                'lambda==5~sigma[base]==0.1','lambda==5~sigma[base]==1', 'lambda==5~sigma[base]==2', 
                               'lambda==7~sigma[base]==0.1','lambda==7~sigma[base]==1', 'lambda==7~sigma[base]==2'))
                               
                            
p7=ggplot(df, aes(x = value))+ 
  geom_histogram(bins = 50)+
  facet_wrap(.~q_values,labeller = label_parsed,ncol=3, scales='free')+
  theme_bw()+
  theme(legend.position='bottom', legend.title = element_blank(), panel.spacing = unit(0.5, "lines"),
        axis.title=element_blank())

p7


ggsave(p7, file='lambda_2parents_10x45.pdf',dpi=300)

