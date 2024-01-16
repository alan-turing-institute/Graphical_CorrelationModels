# raw old  functions
library(INLA)
library(graph)
library(numDeriv)



### GraphDens: function that returns density for a given model structure
## input:
# S: model stucture given as a formula
## output:
# FL: formula for negative logarithm of density
# JD: same but as a function
# STR: model structure
# Pmat: precision matrix
# NP: number of parents
# NC: number of children
GraphDens <- function(S){
  parents <- sapply(S, function(x) strsplit(as.character(x), split="~")[[2]])
  children <- sapply(S, function(x) strsplit(as.character(x), split="~")[[3]])
  childrenU <- sapply(children, function(x) strsplit(gsub("\\s", "", as.character(x)), split="\\+"))
  NPi <- NP <- length(parents) # total number of parents
  FL <- vector("list", length(parents))
  JD <- vector("list", length(parents))
  STR <- vector("list", length(parents))
  for(i in 1:(length(parents))){
    if(i!=1){
      #base = i-1
      parentsrm <- parents[NP] # remove those parents
      parents <- parents[1:(NP-1)] # remaining parents
      for(j in 1:length(parentsrm)){
        rmitem <- which(sapply(childrenU, function(x) parentsrm[j]%in%x)) # remove parentrm j
        childrenU[[rmitem]] <- childrenU[[rmitem]][-which(childrenU[[rmitem]]==parentsrm[j])]
        OldC <- childrenU[[rmitem]][grep("c", childrenU[[rmitem]])]
        NewC <- childrenU[[length(childrenU)]][grep("c", childrenU[[length(childrenU)]])]
        NewChildren <- c(childrenU[[rmitem]], childrenU[[length(childrenU)]])
        childrenU[[rmitem]] <- c(sort(NewChildren[grep("p", NewChildren)]),
                                 sort(NewChildren[grep("c", NewChildren)]))
        childrenU <- childrenU[-length(childrenU)]
      }
      NP <- NP-1 # total number of parents
    }
    STR0_i <- unique(c(parents, unlist(childrenU)))
    #MODIFICATION OF ORDER (for double digits)
    STR_i <- c(STR0_i[grep("p", STR0_i)][order(as.integer(substr(STR0_i[grep("p", STR0_i)], 2, nchar(STR0_i[grep("p", STR0_i)]))))],
               STR0_i[grep("c", STR0_i)][order(as.integer(substr(STR0_i[grep("c", STR0_i)], 2, nchar(STR0_i[grep("c", STR0_i)]))))])
    # STR_i <- c(sort(STR0_i[grep("p", STR0_i)]), sort(STR0_i[grep("c", STR0_i)]))
    edg <- vector("list", length=length(STR_i))
    names(edg) <- STR_i
    NC <- length(STR_i)-NP # total number of children
    fam <- c(STR_i[-c(1:NP)], STR_i[1:NP]) # reorder: children - parents
    FL_i <- "0.5*(" #formula
    if(i==1) Pmat <- matrix(NA, nrow=length(STR_i), ncol=length(STR_i)) # Precision matrix
    if(i==1) colnames(Pmat) <- rownames(Pmat) <- fam
    if(i==1) Pmat[1:NC, 1:NC] <- diag(NC) # children
    for(k in 1:length(parents)){
      for(j in 1:length(childrenU[[k]])){
        edg[[k]] <- list(edges=childrenU[[k]]) #graph
      }
      for(j in 1:NC){ # children - parents
        if(fam[j] %in% childrenU[[k]]){
          if(i==1) Pmat[NC+k, j] <- Pmat[j, NC+k] <- 1 # diagonal
          FL_i <- paste0(FL_i, "(x[",j,"]-x[",NC+k,"])^2+") # formula for joint density - children
        }else{
          if(i==1) Pmat[NC+k, j] <- Pmat[j, NC+k] <- 0 # diagonal
        }
      }
      NCNP <- length(grep("c", childrenU[[k]]))  # number of children not parents themselves
      if(i==1) Pmat[parents[k],parents[k]] <- 1
      if(length(parents)>1){
        for(j in 2:length(parents)){
          if(parents[j]%in%childrenU[[k]]){
            if(i==1) Pmat[parents[k],parents[k]] <- 1 # diagonal for parents
            if(i==1) Pmat[parents[k],parents[j]] <- Pmat[parents[j],parents[k]] <- 1 # off diagonal for parents
            FL_i <- paste0(FL_i, "((x[",NC+j,"]-x[",NC+k,"])^2)/SDev[",j,"]^2+") # formula for joint density - parents
          }else if(j>k){
            if(i==1) Pmat[parents[k],parents[j]] <- Pmat[parents[j],parents[k]] <- 0 # off diagonal for parents
          }
        }
      }
    }
    FL_i <- paste0(FL_i, "(x[",NC+1,"]^2)/SDev[1]^2)") # add common ancestor
    JD_i <- function(x, SDev) {} # build joint density function
    body(JD_i) <- parse(text=FL_i)

    FL[[i]] <- FL_i
    JD[[i]] <- JD_i
    STR[[i]] <- STR_i
  }
  return(list(FL=FL, JD=JD, STR=STR, Pmat=Pmat, NP=NPi, NC=length(grep("c", unlist(childrenU))))) # do we need S in the output?
}









### Tdist: function that returns the distance between two models
## input:
# f1: negative logarithm of density function for model 1
# f0: negative logarithm of density function for model 0 (NULL if model 0 is independent, i.e., last step)
# SD1: standard deviation for model 1
# SD0: standard deviation for model 0 (NULL if model 0 is independent, i.e., last step)
# STR: model structure (i.e., for the number of nodes)
# NC: number of children
# output: distance
Tdist <- function(f1, f0=NULL, SD1, SD0=NULL, STR, NC){ # f negative logarithm of density
  Q1 <- hessian(f1, rep(0, length(STR)), SDev=c(SD0, SD1)) # precision 1
  diag(Q1) <- diag(Q1)+1e-12
  L1 <- chol(Q1)
  VC1 <- chol2inv(L1) # variance-covariance 1
  LUB1 <- VC1[1:NC,1:NC] # left-upper block 1
  Cor1 <- diag(diag(LUB1)^(-1/2)) %*% LUB1 %*% diag(diag(LUB1)^(-1/2)) # correlation 1
  lc1 <- chol(Cor1)
  Dt1 <- 2*sum(log(diag(lc1))) # determinant 1
  if(length(SD0)>0) { # only if model 0 is not identity
    Q0 <- hessian(f0, rep(0, length(STR)-1), SDev=SD0) # precision 0
    diag(Q0) <- diag(Q0)+1e-12
    L0 <- chol(Q0)
    VC0 <- chol2inv(L0) # variance-covariance 0
    LUB0 <- VC0[1:NC,1:NC] # left-upper block 0
    Cor0 <- diag(diag(LUB0)^(-1/2)) %*% LUB0 %*% diag(diag(LUB0)^(-1/2)) # correlation 0
    lc0 <- chol(Cor0)
    Dt0 <- 2*sum(log(diag(lc0))) # determinant 0
    KLD <- 0.5*(sum(diag(solve(Cor0)%*%Cor1))-length(diag(Cor1))-Dt1+Dt0) # KLD
    Dst <- sqrt(2*KLD) # distance
  }else{
    Dst <- sqrt(max(0, -Dt1))##-log(min(Dt1, 1))) # distance
  }
  return(Dst)
}





### GraphPrior: computes the prior distribution for a given model structure
## input:
# S: model stucture given as a formula
# lat: latent parameter
# lambda: parameter of the exponential prior on the distance
# SP: GraphDens applied to S (avoids to recompute it)
# Tdist: function that returns the distance between two models
## output:
# PRT: prior for theta
GraphPrior <- function(S, lat, lambda, SP, Tdist){
  #lambda=lambda*sqrt(SP$NP)
  DIS <- vector("list", length=SP$NP)
  PRT <- vector("list", length=SP$NP) # prior theta
  if(SP$NP>1){
    for(i in 1:(SP$NP-1)){
      SDBase <- exp(lat[-c((length(lat)-i+1):length(lat))]) # standard deviation of model 0
      SDNew <- exp(lat[length(lat)-i+1]) # standard deviation of model 1
      DIS[[i]] <- Tdist(f1=SP$JD[[i]], f0=SP$JD[[i+1]], SD1=SDNew, SD0=SDBase, STR=SP$STR[[i]], NC=SP$NC)
      fq <- function(SDNew){
        return(Tdist(f1=SP$JD[[i]], f0=SP$JD[[i+1]], SD1=SDNew, SD0=SDBase, STR=SP$STR[[i]], NC=SP$NC))
      }
      if(log(SDNew)<(-4.5/sqrt(SP$NP-i+1))){
        SPF <- splinefun(c(-sqrt(c(4,9))/sqrt(SP$NP-i+1)), c(log(abs(grad(fq,exp(-sqrt(4)/sqrt(SP$NP-i+1))))), log(abs(grad(fq,exp(-sqrt(9)/sqrt(SP$NP-i+1)))))))
        LGD <- SPF(log(SDNew))
      }else{
        LGD <- log(abs(grad(fq,SDNew)))
      }
      # LGD <- log(abs(grad(fq,SDNew)))
      # N=1000
      # GP <- NULL
      # tt=seq(-5,5, len=N)
      # GP <- sapply(tt, function(x) log(lambda) - lambda*DIS[[i]] +log(abs(grad(fq,exp(x)))))
      # dev.new()
      # plot(tt, GP, pch=19)
      # dev.off()
      # browser()
      Prd <- log(lambda) - lambda*DIS[[i]] + LGD # prior for distance (log scale)
      PRT[[i]] <- Prd # prior for theta
    }
  }
  DIS[[SP$NP]] = Tdist(f1=SP$JD[[SP$NP]], SD1=exp(lat[SP$NP]), STR=SP$STR[[SP$NP]], NC=SP$NC)
  fqbase <- function(SDNew){ # last step for identity model (special case)
    return(Tdist(f1=SP$JD[[SP$NP]], SD1=SDNew, STR=SP$STR[[SP$NP]], NC=SP$NC))
  }
  if(lat[SP$NP]< (-4.5)){
    SPF <- splinefun(c(-2, -4), c(log(abs(grad(fqbase,exp(-2)))), log(abs(grad(fqbase,exp(-4))))))
    LGD <- SPF(lat[SP$NP])
  }else{
    LGD <- log(abs(grad(fqbase,exp(lat[SP$NP]))))
  }
  Prd <- log(lambda) - lambda*DIS[[SP$NP]] + LGD
  PRT[[SP$NP]] <- Prd
  return(PRT)
}


######R generic

CorGraphs.model <- function (cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                     "log.prior", "quit"), theta = NULL)
{
  envir <- parent.env(environment())
  library(Rgraphviz)
  library(numDeriv)
  if (exists("GS", envir=envir)) {
    GS <- get("GS", envir=envir) # need to get or already in environment?
  } else {
    GS <- Argm$SP # graph structure
    assign("GS", GS, envir=envir)
  }

  graph <- function() {
    return(Q()) # does not depend on theta
  }

  Q <- function() {
    q <- exp(theta[-(1:Argm$SP$NC)])
    Q <- hessian(GS$JD[[1]], rep(0, length(GS$STR[[1]])), SDev=q)*GS$Pmat # precision matrix
    diag(Q) <- diag(Q)+1e-12
    VC <- solve(Q) # variance-covariance
    LUB <- VC[1:Argm$SP$NC,1:Argm$SP$NC] # left-upper block
    Cor <- diag(diag(LUB)^(-1/2)) %*% LUB %*% diag(diag(LUB)^(-1/2)) # correlation
    SD <- diag(exp(-1/2*theta[1:Argm$SP$NC]))
    COV <- SD %*% Cor %*% SD
    return (solve(COV))
  }

  mu <- function() {
    return(numeric(0))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    q <- theta[-(1:Argm$SP$NC)]
    val <- (sum(dgamma(exp(theta[1:Argm$SP$NC]), shape = 1, rate = 0.01, log = TRUE)) + sum(theta[1:Argm$SP$NC]) +
              sum(unlist(Argm$GraphPrior(Argm$S, q, Argm$lambda, Argm$SP, Argm$Tdist)))+sum(theta[-(1:Argm$SP$NC)]))
    return(val)
  }

  initial <- function() {
    return(c(rep(4, Argm$SP$NC), rep(Argm$init, Argm$SP$NP)))
  }

  quit <- function() {
    return(invisible())
  }

  if (!length(theta)) {
    theta <- initial()
  }

  val <- do.call(match.arg(cmd), args = list())
  return(val)
}
















# function to map from theta to correlation matrix
ThetaCor <- function(SP, lat, COV=F){
  # remove round
  Q <- hessian(SP$JD[[1]], rep(0, length(SP$STR[[1]])), SDev=exp(lat[-(1:SP$NC)]))#*SP$Pmat # precision matrix
  diag(Q) <- diag(Q)+1e-12
  VC <- solve(Q) # variance-covariance #chol2solve?
  LUB <- VC[1:SP$NC,1:SP$NC] # left-upper block
  Cor <- diag(diag(LUB)^(-1/2)) %*% LUB %*% diag(diag(LUB)^(-1/2)) # correlation
  SD <- diag(exp(-1/2*lat[1:SP$NC]))
  if(COV){
    MatOut <- SD %*% Cor %*% SD
  }else{
    MatOut <- Cor
    diag(MatOut) <- diag(SD)
  }
  return(MatOut)
}



# function to plot a graph
GraphPlot <- function(S, base=0, fontsize=c(14, 14), width=c(0.75, 0.75), height=c(0.5,0.5), col=c("red", "limegreen", "lightsalmon", "lightgreen")){
  parents <- sapply(S, function(x) strsplit(as.character(x), split="~")[[2]])
  children <- sapply(S, function(x) strsplit(as.character(x), split="~")[[3]])
  childrenU <- sapply(children, function(x) strsplit(gsub("\\s", "", as.character(x)), split="\\+"))
  NP <- length(parents) # total number of parents
  if(base!=0){
    parentsrm <- parents[NP:(NP-base+1)] # remove those parents
    parents <- parents[1:(NP-base)] # remaining parents
    for(i in 1:length(parentsrm)){
      # remove parentrm i
      rmitem <- which(sapply(childrenU, function(x) parentsrm[i]%in%x))
      childrenU[[rmitem]] <- childrenU[[rmitem]][-which(childrenU[[rmitem]]==parentsrm[i])]
      NewChildren <- c(childrenU[[rmitem]], childrenU[[length(childrenU)]])
      childrenU[[rmitem]] <- c(sort(NewChildren[grep("p", NewChildren)]), sort(NewChildren[grep("c", NewChildren)]))
      childrenU <- childrenU[-length(childrenU)]
    }
    NP <- NP-base # total number of parents
  }
  nod0 <- unique(c(parents, unlist(childrenU)))
  nod <- c(nod0[grep("p", nod0)][order(as.integer(substr(nod0[grep("p", nod0)], 2, nchar(nod0[grep("p", nod0)]))))],
           nod0[grep("c", nod0)][order(as.integer(substr(nod0[grep("c", nod0)], 2, nchar(nod0[grep("c", nod0)]))))])
  edg <- vector("list", length=length(nod))
  names(edg) <- nod
  nAttrs <- NULL # for plots
  nAttrs$color <- c(sapply(nod[grep("p", nod)], function(x) x=col[1]), sapply(nod[grep("c", nod)], function(x) x=col[2]))
  nAttrs$fillcolor <- c(sapply(nod[grep("p", nod)], function(x) x=col[3]), sapply(nod[grep("c", nod)], function(x) x=col[4]))
  nAttrs$shape <- c(sapply(nod[grep("p", nod)], function(x) x="box"), sapply(nod[grep("c", nod)], function(x) x="circle"))
  nAttrs$height <- c(sapply(nod[grep("p", nod)], function(x) x=height[1]), sapply(nod[grep("c", nod)], function(x) x=height[2]))
  nAttrs$width <- c(sapply(nod[grep("p", nod)], function(x) x=width[1]), sapply(nod[grep("c", nod)], function(x) x=width[2]))
  nAttrs$fontsize <- c(sapply(nod[grep("p", nod)], function(x) x=fontsize[1]), sapply(nod[grep("c", nod)], function(x) x=fontsize[2]))
  NC <- length(nod)-NP # total number of children
  # build precision matrix
  fam <- c(nod[-c(1:NP)], nod[1:NP]) # reorder: children - parents
  for(i in 1:length(parents)){
    for(j in 1:length(childrenU[[i]])){
      edg[[i]] <- list(edges=childrenU[[i]]) #graph
    }
  }
  gr <- graphNEL(nodes=nod, edgeL=edg, edgemode='directed')
  return(list(gr=gr, nAttrs=nAttrs))
}


GraphPlotPrior <- function(S, fontsize=c(14, 14), width=c(0.75, 0.75), height=c(0.5,0.5), fontsizeLast=14, widthLast=0.75, heightLast=0.5){
  parents <- sapply(S, function(x) strsplit(as.character(x), split="~")[[2]])
  children <- sapply(S, function(x) strsplit(as.character(x), split="~")[[3]])
  childrenU <- sapply(children, function(x) strsplit(gsub("\\s", "", as.character(x)), split="\\+"))
  nod0 <- unique(c(parents, unlist(childrenU)))
  nod <- c(nod0[grep("p", nod0)][order(as.integer(substr(nod0[grep("p", nod0)], 2, nchar(nod0[grep("p", nod0)]))))],
           nod0[grep("c", nod0)][order(as.integer(substr(nod0[grep("c", nod0)], 2, nchar(nod0[grep("c", nod0)]))))])
  GRAPHS <- vector("list", length=length(parents)+1)
  for(i in 0:(length(parents)-1)){
      SP <- GraphPlot(S, base=i, fontsize=fontsize, width=width, height=height)
      GRAPHS[[i+1]] <- list("gr"=SP$gr, "nAttrs"=SP$nAttrs)
  }
  nAttrs <- NULL # for plots
  nAttrs$color <- sapply(nod[grep("c", nod)], function(x) x="limegreen")
  nAttrs$fillcolor <- sapply(nod[grep("c", nod)], function(x) x="lightgreen")
  nAttrs$shape <- sapply(nod[grep("c", nod)], function(x) x="circle")
  nAttrs$height <- sapply(nod[grep("c", nod)], function(x) x=heightLast)
  nAttrs$width <- sapply(nod[grep("c", nod)], function(x) x=widthLast)
  nAttrs$fontsize <- sapply(nod[grep("c", nod)], function(x) x=fontsizeLast)
  GRAPHS[[length(parents)+1]] <- list("gr"=graphNEL(nodes=nod[(length(parents)+1):length(nod)]), "nAttrs"=nAttrs)
  return(GRAPHS)
}











