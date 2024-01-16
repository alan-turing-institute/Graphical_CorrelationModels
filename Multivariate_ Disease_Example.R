############################################
#Chunk code:  library 
############################################
library(sf)
library(spdep)
library(ggplot2)
library(patchwork)
library(INLA)

############################################
#Chunk code:  settings
############################################
### inlasetup
inla.setOption(
    num.threads = "6:1",
    safe = FALSE
)

### inla setup controls
ctrc <- list(
    dic = TRUE, waic = TRUE, cpo = TRUE,
    control.gcpo = list(
        enable = TRUE,
        num.level.sets = 5,
        strategy = "posterior")
)
ctri <- list(
    int.strategy = "ccd"
)




############################################
#Chunk code:  sourcing extra functions 
############################################

source('RGEN.R')
source("mbesag.R")




############################################
#Chunk code:  load data
############################################


load('GermanyData.Rdata')

############################################
#Chunk code: dat
############################################

## vnames
Vnames <- c("Oral", "Osph", "Lary", "Lung")
names(Vnames) <- Vnames
table(mydatalong$idx <- factor(gsub("Osp", "Osph", mydatalong$idx), Vnames))
(vnames <- tolower(Vnames))
(K <- length(vnames))

### Summary
tapply(mydatalong$Obs / mydatalong$E,
       mydatalong$idx, summary)

### limitis for ploting SMR
r.smr <- c(0, 2.75)

### observed SMR covariance, SD and correlation 
scc.fn <- function(x, mat = FALSE) {
    if(mat) {
        r <- x
    } else {
        r <- cov(x)
    }
    r[upper.tri(r)] <- cov2cor(r)[upper.tri(r)]
    diag(r) <- sqrt(diag(r))
    return(r)
}

n <- nrow(myshp)
sdCovCor.obs <- scc.fn(matrix(mydatalong$Obs / mydatalong$E, n))
cc.dn <- list(Vnames, Vnames)
dimnames(sdCovCor.obs) <- cc.dn
round(sdCovCor.obs, 2) ## upper.tri as in Table 2 of Held et. al. (2005)

### creating graph
nbl <- poly2nb(myshp)
stopifnot(n == length(nn <- card(nbl)))
graph <- sparseMatrix(
    i = rep(1:n, nn),
    j = unlist(nbl[nn>0]),
    x = 1L)

### fix columnames
names(myshp@data) <- gsub(
    "osp_", "osph_", gsub("exp_osp", "osp_exp", names(myshp@data)))
names(myshp)

### Models 1: scaled BYM model for each disease
pcprec <- list(prior = "pc.prec", param = c(1, 0.01))
m1f <- Obs ~
    f(id.area, model = "bym2", graph = graph)

lres.m1 <-lapply(vnames, function(v) {
    inla(formula = m1f, 
         family = "poisson",
         data = data.frame(
             Obs = myshp@data[, paste0(v, "_obs")],
             id.area = myshp$id,
             expected = myshp@data[, paste0(v, "_exp")]), 
         E = expected, 
         control.inla = ctri,
         control.compute = ctrc
         )
})

sdCovCor.m1 <- scc.fn(sapply(
    lres.m1, function(r) r$summary.random$id$mean))
dimnames(sdCovCor.m1) <- cc.dn
round(sdCovCor.m1, 2)

## Model 2: multivariate iid model 
mcov.model <- inla.rgeneric.define( ## use LKJ prior
    rgeneric.MBesag, 
    Args = list(NComps = K, Structure = Diagonal(n),
                LKJeta = 1, pcprec.lambda = -log(0.01)/1)
)

m2f <- Obs ~ 0 + idx +
    f(id, model = mcov.model) ##"iidkd", order = K, n = n*K)

res.m2 <- inla(
    formula = m2f,
    family = "poisson",
    data = mydatalong, 
    E = E,
    verbose = !TRUE,
    control.compute = ctrc,
    control.inla = ctri
)

res.m2$cpu

sdCovCor.m2 <- scc.fn(interpret.ptheta(
    res.m2$mode$theta, K)$covar, mat = TRUE)
dimnames(sdCovCor.m2) <- cc.dn
round(sdCovCor.m2, 2)

############################################
#Chunk code:  Model 2: multivariate Besag model
############################################


R.unscaled <- Diagonal(n, nn) - graph
Rscaled <- inla.scale.model(
    Q = R.unscaled,
    constr = list(A = matrix(1, 1, n), e = 0))

mcar.model <- inla.rgeneric.define(
    rgeneric.MBesag, 
    Args = list(NComps = K, Structure = Rscaled,
                LKJeta = 1, pcprec.lambda = -log(0.01)/1)
)

m3f <- Obs ~ 0 + idx +
    f(id, model = mcar.model, 
      extraconstr = list(
          A = kronecker(diag(K), matrix(1, 1, n)),
          e = rep(0, K))
      ) 

### fit the model
res.m3 <- inla(
    formula = m3f, 
    family = "poisson",
    data = mydatalong, 
    E = E,
    verbose = !TRUE,
    control.mode = list(theta = res.m2$mode$theta, restart = TRUE),
    control.compute = ctrc,
    control.inla = ctri
)

res.m3$cpu

sdCovCor.m3 <- scc.fn(interpret.ptheta(
    res.m3$mode$theta, K)$covar, mat = TRUE)
dimnames(sdCovCor.m3) <- cc.dn
round(sdCovCor.m3, 2)

#############################
### Model 4: graph model 
#############################

S <- list(p1 ~ p2 + c4, p2 ~ p3 + c2 + c3, p3 ~ c1)

SP <- GraphDens(S)
SP_plot <- GraphPlot(S, base=0)

par(mar = c(1, 1, 1, 1))
plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)

#############################
# Chunk Code: parameters settings
#############################
lambda <- 5
Args <- list(S=S, lambda=lambda, SP=SP,
             Tdist=Tdist, GraphPrior=GraphPrior, init=-1)
Args$SP$Structure <- Rscaled

GBesag <- inla.rgeneric.define(
    rgeneric.CorGraphsBesag, Argm=Args)

### fit
m4f <- Obs ~ 0 + idx + 
    f(id, model = GBesag, 
      extraconstr = list(
          A = kronecker(diag(K), matrix(1, 1, n)),
          e = rep(0, K))
      ) 

### fit the model
res.m4 <- inla(
    formula = m4f, 
    family = "poisson",
    data = mydatalong, 
    E = E,
    verbose = !TRUE, 
    control.compute = ctrc,
    control.inla = ctri
)

res.m4$cpu

sdCovCor.m4 <- scc.fn(ThetaCor(
    SP, res.m4$mode$theta, COV = TRUE), mat = TRUE)
dimnames(sdCovCor.m4) <- cc.dn
round(sdCovCor.m4, 2)

### correlation of the posterior mean
round(scc.fn(matrix(res.m4$summary.random$id$mean, n)), 2)

#######################################
### Model 5: simpler graph model 
#######################################

S2 <- list(p1 ~ p2 + c4, p2 ~ c1 + c2 + c3)
SP2 <- GraphDens(S2)
SP2_plot <- GraphPlot(S2, base=0)

par(mar = c(1, 1, 1, 1))
plot(SP2_plot$gr, nodeAttrs = SP2_plot$nAttrs)

#############################
# Chunk Code: parameters settings
#############################
lambda2 <- lambda
Args2 <- list(S=S2, lambda=lambda2, SP=SP2,
             Tdist=Tdist, GraphPrior=GraphPrior, init=-1)
Args2$SP$Structure <- Rscaled

GBesag2 <- inla.rgeneric.define(
    rgeneric.CorGraphsBesag, Argm=Args2)

### fit
m5f <- Obs ~ 0 + idx + 
    f(id, model = GBesag2, 
      extraconstr = list(
          A = kronecker(diag(K), matrix(1, 1, n)),
          e = rep(0, K))
      ) 

### fit the model
res.m5 <- inla(
    formula = m5f, 
    family = "poisson",
    data = mydatalong, 
    E = E,
    verbose = !TRUE, 
    control.compute = ctrc,
    control.inla = ctri
)

res.m5$mode$theta
res.m5$cpu

sdCovCor.m5 <- scc.fn(ThetaCor(
    SP2, res.m5$mode$theta, COV = TRUE), mat = TRUE)
dimnames(sdCovCor.m5) <- cc.dn
round(sdCovCor.m5, 2)

#############################
# Chunk Code:  DIC, WAIC, LOO-CPO and LGO-CPO 
#############################

utab <- data.frame(
    DIC = c(MA = sum(sapply(lres.m1, function(r) r$dic$dic)),
            apply(cbind(MB = res.m2$dic$local.dic, 
                        MC = res.m3$dic$local.dic, 
                        MD = res.m4$dic$local.dic,
                        ME = res.m5$dic$local.dic), 2, sum)),
    WAIC = c(MA = sum(sapply(lres.m1, function(r) r$waic$waic)),
             apply(cbind(MB = res.m2$waic$local.waic, 
                         MC = res.m3$waic$local.waic, 
                         MD = res.m4$waic$local.waic,
                         ME = res.m5$waic$local.waic), 2, sum)),
    CPO = c(-sum(sapply(lres.m1, function(r)
        log(r$cpo$cpo))), 
        -sum(log(res.m2$cpo$cpo)), 
        -sum(log(res.m3$cpo$cpo)), 
        -sum(log(res.m4$cpo$cpo)),
        -sum(log(res.m5$cpo$cpo))),
    GCPO = c(-sum(sapply(lres.m1, function(r)
        log(r$gcpo$gcpo))), 
        -sum(log(res.m2$gcpo$gcpo)), 
        -sum(log(res.m3$gcpo$gcpo)), 
        -sum(log(res.m4$gcpo$gcpo)),
        -sum(log(res.m5$gcpo$gcpo))))
utab

xtable::xtable(utab,
       caption = "DIC, WAIC, CPO and GCPO for the four models fitted.")

    
#############################
# Chunk Code:  tables
#############################

dic.tab <- t(cbind(
    MA = sapply(lres.m1, function(r) r$dic$dic),
    apply(cbind(MB = res.m2$dic$local.dic, 
                MC = res.m3$dic$local.dic, 
                MD = res.m4$dic$local.dic,
                ME = res.m5$dic$local.dic), 2,
          tapply, rep(1:K, each = n), sum)))

dic.tab

waic.tab <- t(cbind(
    MA = sapply(lres.m1, function(r) r$waic$waic),
    apply(cbind(MB = res.m2$waic$local.waic, 
                MC = res.m3$waic$local.waic, 
                MD = res.m4$waic$local.waic,
                ME = res.m5$waic$local.waic), 2,
          tapply, rep(1:K, each = n), sum)))

waic.tab

cpo.tab <- t(cbind(
    MA= sapply(lres.m1, function(r)
        -sum(log(r$cpo$cpo))), 
    apply(cbind(MB = res.m2$cpo$cpo, 
                MC = res.m3$cpo$cpo, 
                MD = res.m4$cpo$cpo,
                ME = res.m5$cpo$cpo), 2,
          tapply, rep(1:K, each = n),
          function(x) -sum(log(x)))))

cpo.tab

gcpo.tab <- t(cbind(
    MA= sapply(lres.m1, function(r)
        -sum(log(r$gcpo$gcpo))), 
    apply(cbind(MB = res.m2$gcpo$gcpo,
                MC = res.m3$gcpo$gcpo, 
                MD = res.m4$gcpo$gcpo,
                ME = res.m5$gcpo$gcpo), 2,
          tapply, rep(1:K, each = n),
          function(x) -sum(log(x)))))

gcpo.tab

#############################
# Chunk Code:  join tables 
#############################

tabs <- rbind(cbind(dic.tab, waic.tab),
              cbind(cpo.tab, gcpo.tab))
nm <- nrow(dic.tab)
rownames(tabs) <- paste0(
    "$",
    paste0("M_", c(LETTERS[1:nm],
                   letters[1:nm])), 
    "$")

## xtable
xtb <- xtable::xtable(tabs, digits = rep(0:1, c(1, 8)),
                      align = "lrrrrrrrr",
                      caption = "DIC, WAIC, CPO and GCPO for the four models fitted.")



## xtable print for the paper (change some rownames and hline position)
print(xtb, hline.after = c(-1, 0, c(1,2) * nm),
      add.to.row = list(
          pos = list(0, nm),
          command = c(paste(paste0(
              "& \\multicolumn{", K, "}{c}{", c("DIC", "WAIC"), '}',
              collapse=''), '\\\\\n'), 
              paste(paste0("& \\multicolumn{", K, "}{c}{", c("CPO", "GCPO"), '}',
                           collapse=''), '\\\\\n'))),
      sanitize.text.function = function(x) x
      )

tab.models <- data.frame(
    Disease = c("Independent", "Unstructured correlation",
                "Unstructured correlation", "Graph model",
                "Simplified graph model"),
    Spatial = c("Correlated (BYM)", "Independent",
                rep("Correlated (Besag)", 3)),
    Parameters = c(2*4, 4+6, 4+6, 4+3, 4+2),
    row.names = LETTERS[1:5])

xtable::xtable(tab.models,
               caption = "Synthetic description of the fitted models.")

#############################
# Chunk Code:  correlations
#############################
nsampl <- 3e3
th.m4.samples <- t(inla.hyperpar.sample(
    n = nsampl, res.m4, intern = TRUE, improve = TRUE))

th.m5.samples <- t(inla.hyperpar.sample(
    n = nsampl, res.m5, intern = TRUE, improve = TRUE))

system.time(scc.m4.s <- t(apply(th.m4.samples, 2, function(x) 
    scc.fn(ThetaCor(SP, x, COV = TRUE), mat = TRUE))))

system.time(scc.m5.s <- t(apply(th.m5.samples, 2, function(x) 
    scc.fn(ThetaCor(SP2, x, COV = TRUE), mat = TRUE))))

dim(scc.m4.s)
dim(scc.m5.s)

scc4.li <- matrix(apply(scc.m4.s, 2, quantile, 0.025), 4)
scc4.ls <- matrix(apply(scc.m4.s, 2, quantile, 0.975), 4)

round(sdCovCor.m4, 2)
round(scc4.li, 2)
round(scc4.ls, 2)

scc5.li <- matrix(apply(scc.m5.s, 2, quantile, 0.025), 4)
scc5.ls <- matrix(apply(scc.m5.s, 2, quantile, 0.975), 4)

round(sdCovCor.m5, 2)
round(scc5.li, 2)
round(scc5.ls, 2)


#############################
# Chunk Code:  plot prep
#############################

### Convert to sf format 
selcolumns <- c("id", paste0(rep(vnames, each = 3),
                             c("_obs", "_exp", "_smr")))
dmap <- st_as_sf(myshp[, selcolumns])

### map visualize
ggmap0 <- ggplot(data = dmap) +
    theme_minimal() + xlab("") + ylab("") +
    theme(legend.position = c(0.95, 0.25),
          axis.text.x=element_blank(),  ##remove x axis labels
          axis.text.y=element_blank())
gg.smr <- list(
    ggmap0 + geom_sf(aes(fill = oral_smr),
                     col = "transparent")+
    scale_fill_viridis_c(
        name = "Obs.\nOral\nSMR",
        limits = c(0, 2.75)
    ),
    ggmap0 + geom_sf(aes(fill = osph_smr),
                     col = "transparent")+
    scale_fill_viridis_c(
        name = "Obs.\nOsph.\nSMR",
        limits = c(0, 2.75)
    ),
    ggmap0 + geom_sf(aes(fill = lary_smr),
                     col = "transparent")+
    scale_fill_viridis_c(
        name = "Obs.\nLary.\nSMR",
        limits = c(0, 2.75)
    ),
    ggmap0 + geom_sf(aes(fill = lung_smr),
                     col = "transparent")+
    scale_fill_viridis_c(
        name = "Obs.\nLung\nSMR",
        limits = c(0.5, 2)
    ))
               
wrap_plots(gg.smr, ncol = 2)

lsmr.m4est <- matrix(res.m4$summary.random$id$mean, n)
summary(lsmr.m4est)
summary(exp(lsmr.m4est))

dmap$oral_est <- exp(lsmr.m4est[, 1])
dmap$osph_est <- exp(lsmr.m4est[, 2])
dmap$lary_est <- exp(lsmr.m4est[, 3])
dmap$lung_est <- exp(lsmr.m4est[, 4])

ggmap1 <- ggplot(data = dmap) +
    theme_minimal() + xlab("") + ylab("") +
    theme(legend.position = c(0.95, 0.25),
          axis.text.x=element_blank(),  ##remove x axis labels
          axis.text.y=element_blank()) 

gg.smr.est <- list(
    ggmap1 + geom_sf(aes(fill = oral_est),
                     col = "transparent") +
    scale_fill_viridis_c(
        name = "Fitted\nOral\nSMR",
        limits = c(1/1.9, 1.9)),
    ggmap1 + geom_sf(aes(fill = osph_est),
                     col = "transparent")+
    scale_fill_viridis_c(
        name = "Fitted\nOsph.\nSMR",
        limits = c(1/1.9, 1.9)
    ),
    ggmap1 + geom_sf(aes(fill = lary_est),
                     col = "transparent")+
    scale_fill_viridis_c(
        name = "Fitted\nLary.\nSMR",
        limits = c(1/1.9, 1.9)
    ),
    ggmap1 + geom_sf(aes(fill = lung_est),
                     col = "transparent")+
    scale_fill_viridis_c(
        name = "Fitted\nLung\nSMR",
        limits = c(1/1.9, 1.9)
    ))

png("smrmaps.png",
    width = 4500,
    height = 3000,
    res = 300)

wrap_plots(c(gg.smr, gg.smr.est), ncol = 4)

dev.off()

if(FALSE)
    system("eog smrmaps.png &")


#############################
# Chunk Code: plot in ggplot
#############################
library(ggpubr)
library(patchwork)

 
d1=ggmap0 + geom_sf(aes(fill = oral_smr),
                 col = "transparent")+
  scale_fill_viridis_c(
    name = "SMR",
    limits = c(0, 2.75))+ggtitle('Oral')+
   ylab('SMR')


d2=ggmap0 + geom_sf(aes(fill = osph_smr),
                 col = "transparent")+
  scale_fill_viridis_c(
    name = "Obs.\nOsph.\nSMR",
    limits = c(0, 2.75))+ggtitle('Osph')
 
d3=ggmap0 + geom_sf(aes(fill = lary_smr),
                 col = "transparent")+
  scale_fill_viridis_c(
    name = "Obs.\nLary.\nSMR",
    limits = c(0, 2.75))+ggtitle('Lary')

d4=ggmap0 + geom_sf(aes(fill = lung_smr),
                 col = "transparent")+
  scale_fill_viridis_c(
    name = "Obs.\nLung\nSMR",
    limits = c(0, 2.75))+ggtitle('Lung')



g1=ggarrange(d1, d2, d3, d4, ncol=4, nrow=1, common.legend = TRUE, legend="right")
g1


p1=ggmap1 + geom_sf(aes(fill = oral_est),
                 col = "transparent") +
  scale_fill_viridis_c(
    name = "Est. \nSMR",
    limits = c(0,2.75))+ylab('Estimated SMR')

p2=ggmap1 + geom_sf(aes(fill = osph_est),
                 col = "transparent")+
  scale_fill_viridis_c(
    name = "Fitted\nOsph.\nSMR",
    limits = c(0, 2.75))

p3=ggmap1 + geom_sf(aes(fill = lary_est),
                 col = "transparent")+
  scale_fill_viridis_c(
    name = "Fitted\nLary.\nSMR",
    limits = c(0, 2.75))
  
p4=ggmap1 + geom_sf(aes(fill = lung_est),
                 col = "transparent")+
  scale_fill_viridis_c(
    name = "Fitted\nLung\nSMR",
    limits = c(0, 2.75))



g2=ggarrange(d1, d2, d3, d4,p1, p2, p3, p4, ncol=4, nrow=2, common.legend = TRUE, legend="right", label.y ='TOP')
g2


ggsave(file='smrmaps.png' , g2, dpi=600)



combined <- d1 + d2+d3+d4 +p1+p2+p3+p4 & theme(legend.position = "right")
combined + plot_layout(ncol=4, guides = "collect")

