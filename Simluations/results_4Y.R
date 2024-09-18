# print results of simulations
# setwd()

load("./4Y/4Y_30_L1.RData")
DIC_ref <- median(DIC)
WAIC_ref <- median(WAIC)
res_4Y_30_L1 <- rbind(cbind(apply(res_sim, 2, mean),
                            apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                            apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                      quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                      quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                      quantile(CT, c(0.5, 0.025, 0.975)))
load("./4Y/4Y_30_L3.RData")
res_4Y_30_L3 <- rbind(cbind(apply(res_sim, 2, mean),
                            apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                            apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                      quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                      quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                      quantile(CT, c(0.5, 0.025, 0.975)))
load("./4Y/4Y_30_L5.RData")
res_4Y_30_L5 <- rbind(cbind(apply(res_sim, 2, mean),
                            apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                            apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                      quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                      quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                      quantile(CT, c(0.5, 0.025, 0.975)))
load("./4Y/4Y_30_FULL.RData")
res_4Y_30_FULL <- unname(rbind(cbind(apply(res_sim, 2, mean),
                                     apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                                     apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                               quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                               quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                               quantile(CT, c(0.5, 0.025, 0.975))))
load("./4Y/4Y_1000_L1.RData")
DIC_ref <- median(DIC)
WAIC_ref <- median(WAIC)
res_4Y_1000_L1 <- rbind(cbind(apply(res_sim, 2, mean),
                        apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                        apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                        quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                        quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                        quantile(CT, c(0.5, 0.025, 0.975)))
load("./4Y/4Y_1000_L3.RData")
res_4Y_1000_L3 <- rbind(cbind(apply(res_sim, 2, mean),
                       apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                       apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                       quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                       quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                       quantile(CT, c(0.5, 0.025, 0.975)))
load("./4Y/4Y_1000_L5.RData")
res_4Y_1000_L5 <- rbind(cbind(apply(res_sim, 2, mean),
                        apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                        apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                        quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                        quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                        quantile(CT, c(0.5, 0.025, 0.975)))
load("./4Y/4Y_1000_FULL.RData")
res_4Y_1000_FULL <- unname(rbind(cbind(apply(res_sim, 2, mean),
                              apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                              apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                        quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                        quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                        quantile(CT, c(0.5, 0.025, 0.975))))

True_4Y <- c("$\\sigma_{c1}$=1 & ",
             "$\\sigma_{c2}$=0.5 & ",
             "$\\sigma_{c3}$=1 & ",
             "$\\sigma_{c4}$=0.5 & ",
             "$\\sigma_{c5}$=1 & ",
             "$\\sigma_{c6}$=0.5 & ",
             "$\\sigma_{c7}$=1 & ",
             "$\\sigma_{c8}$=0.5 & ",
             "$\\rho_1$=0.9 & ",
             "$\\rho_2$=0.9 & ",
             "$\\rho_3$=0.9 & ",
             "$\\rho_4$=0.9 & ",
             "$\\rho_5$=0.8 & ",
             "$\\rho_6$=0.8 & ",
             "$\\rho_7$=0.6 & ",
             "$\\rho_8$=0.6 & ",
             "$\\rho_9$=0.6 & ",
             "$\\rho_10$=0.6 & ",
             "\\hdashline $\\beta_{01}$=0.2 & ",
             "$\\beta_{11}$=-0.1 & ",
             "$\\beta_{21}$=-0.2 & ",
             "$\\beta_{31}$=0.1 & ",
             "$\\beta_{02}$=0.2 & ",
             "$\\beta_{12}$=-0.1 & ",
             "$\\beta_{22}$=-0.2 & ",
             "$\\beta_{32}$=0.1 & ",
             "$\\beta_{03}$=0.2 & ",
             "$\\beta_{13}$=-0.1 & ",
             "$\\beta_{23}$=-0.2 & ",
             "$\\beta_{33}$=0.1 & ",
             "$\\beta_{04}$=0.2 & ",
             "$\\beta_{14}$=-0.1 & ",
             "$\\beta_{24}$=-0.2 & ",
             "$\\beta_{34}$=0.1 & ",
             "\\hdashline DIC & ",
             "WAIC & ",
             "Comp. time (s) & ")

res_4Y_30_L1F <- apply(round(res_4Y_30_L1, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_4Y_30_L3F <- apply(round(res_4Y_30_L3, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_4Y_30_L5F <- apply(round(res_4Y_30_L5, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_4Y_1000_L1F <- apply(round(res_4Y_1000_L1, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_4Y_1000_L3F <- apply(round(res_4Y_1000_L3, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_4Y_1000_L5F <- apply(round(res_4Y_1000_L5, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_4Y_30_FULLF <- apply(round(res_4Y_30_FULL, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_4Y_1000_FULLF <- apply(round(res_4Y_1000_FULL, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))

tab_4Y_30 <- cbind(paste0(True_4Y, res_4Y_30_L1F, " & ",res_4Y_30_L3F,
                          " & ", res_4Y_30_L5F, " & ", res_4Y_30_FULLF, " & \\\\\n"))
cat(tab_4Y_30)

tab_4Y_1000 <- cbind(paste0(True_4Y, res_4Y_1000_L1F, " & ",res_4Y_1000_L3F,
                            " & ", res_4Y_1000_L5F, " & ", res_4Y_1000_FULLF, " & \\\\\n"))
cat(tab_4Y_1000)








