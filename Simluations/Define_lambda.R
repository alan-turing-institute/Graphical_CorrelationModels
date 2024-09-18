# ## PC prior
# source("./RGEN.R")
library(corGraphs)
# The base model is defined by setting the last latent's standard deviation to zero.
#
# We use function Tdist() that computes the distance between two models:


# We need to define lambda parameter for the exponential prior for the distance.
#
# We can sample distances and map them to correlations for a given lambda.
#
# An appropriate choice is a lambda that contracts towards 0.
#
# We can start with a simple 1 parent / 2 children situation, where the base model is 2 independent children:
# our parameter is scaled as log(Standard deviation of the latent) in rgeneric
SDsim <- exp(seq(-5,log(5), len=300))
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

lambda <- c(1, 3, 5)# c(1,4,5,6,7,10)
par(mfrow=c(1,3), mar=c(2,2,2,1))
rst <- vector("list", length=length(lambda))
for(lam in 1:length(lambda)){
  # generate distance as exp function of lambda
  DistSim <- rexp(n = 1000, rate = lambda[lam])
  # SDsim from distance
  ThetasSim <- MapF(DistSim)
  # cor from SDsim
  rst[[lam]] <- sapply(ThetasSim, function(x) ThetaCor(GR, c(10,10,log(x)))[1,2])
  # hist(rst[[lam]], main=paste0("lambda=", lambda[lam]), breaks=50)
  plot(density(rst[[lam]]), main=bquote(lambda == .(lambda[lam])), xlim=c(0,1))
}

# We can do similar, but for the next step (from 1 parent to two parent), then the prior is conditional on the value of the parameter for the first latent:

SDsim0 <- c(0.1, 1, 2)
lambda <- c(1, 3, 5)
hst <- vector("list", length=length(lambda)*length(SDsim0))
it=1
par(mfrow=c(3,3), mar=c(2,2,2,1))
for(lam in 1:length(lambda)){
  for(tet in 1:length(SDsim0)){
    SDsim <- seq(0.01,5, len=50)
    S <- list(p1 ~ p2 + c1, p2 ~ c2 + c3)
    GR <- GraphDens(S)
    lDi <- NULL
    for(tt in 1:length(SDsim)){
      lDi <- c(lDi, Tdist(GR$JD[[1]], GR$JD[[2]], SD1 = SDsim[tt], SD0 = SDsim0[tet], STR=GR$STR[[1]], NC = GR$NC))
    }
    #plot(SDsim, lDi, type="o", main=paste0("SD (base)=", SDsim0[tet]), xlab="SD (new latent)")
    MapF <- splinefun(lDi, SDsim)
    #curve(MapF, from=0, to=2)
    Dsim <- MapF(1)


    # generate distance as exp function of lambda
    DistSim <- rexp(n = 1000, rate = lambda[lam])
    # SDsim from distance
    ThetasSim <- MapF(DistSim)
    # cor from SDsim
    rst <- sapply(ThetasSim, function(x) ThetaCor(GR, c(10,10, 10,log(x), log(SDsim0[tet])))[1,2])
    # hst[[it]] <- hist(rst, main= paste0("lambda=", lambda[lam], " / SD (base)=", SDsim0[tet]),
    #                   breaks=50, xlab="correlation", xlim=c(0,1))
    hst[[it]] <- plot(density(rst), main= bquote(lambda == .(lambda[lam]) ~ sigma[base] == .(SDsim0[tet])),
                      breaks=50, xlab="correlation", xlim=c(0,1))
    it <- it+1
  }
}


































