# print results of simulations
# setwd() 

load("./2Y/2Y_30_L1.RData")
DIC_ref <- median(DIC)
WAIC_ref <- median(WAIC)
res_2Y_30_L1 <- rbind(cbind(apply(res_sim, 2, mean),
                            apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                            apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                      quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                      quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                      quantile(CT, c(0.5, 0.025, 0.975)))
load("./2Y/2Y_30_L3.RData")
res_2Y_30_L3 <- rbind(cbind(apply(res_sim, 2, mean),
                            apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                            apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                      quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                      quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                      quantile(CT, c(0.5, 0.025, 0.975)))
load("./2Y/2Y_30_L5.RData")
res_2Y_30_L5 <- rbind(cbind(apply(res_sim, 2, mean),
                            apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                            apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                      quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                      quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                      quantile(CT, c(0.5, 0.025, 0.975)))
load("./2Y/2Y_30_L3_MIS.RData")
res_2Y_30_L3_MIS <- rbind(cbind(apply(res_sim, 2, mean),
                            apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                            apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                      quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                      quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                      quantile(CT, c(0.5, 0.025, 0.975)))
load("./2Y/2Y_30_FULL.RData")
res_2Y_30_FULL <- unname(rbind(cbind(apply(res_sim, 2, mean),
                                     apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                                     apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                               quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                               quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                               quantile(CT, c(0.5, 0.025, 0.975))))
load("./2Y/2Y_1000_L1.RData")
DIC_ref <- median(DIC)
WAIC_ref <- median(WAIC)
res_2Y_1000_L1 <- rbind(cbind(apply(res_sim, 2, mean),
                              apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                              apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                        quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                        quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                        quantile(CT, c(0.5, 0.025, 0.975)))
load("./2Y/2Y_1000_L3.RData")
res_2Y_1000_L3 <- rbind(cbind(apply(res_sim, 2, mean),
                              apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                              apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                        quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                        quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                        quantile(CT, c(0.5, 0.025, 0.975)))
load("./2Y/2Y_1000_L5.RData")
res_2Y_1000_L5 <- rbind(cbind(apply(res_sim, 2, mean),
                              apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                              apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                        quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                        quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                        quantile(CT, c(0.5, 0.025, 0.975)))
load("./2Y/2Y_1000_L3_MIS.RData")
res_2Y_1000_L3_MIS <- rbind(cbind(apply(res_sim, 2, mean),
                                  apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                                  apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                            quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                            quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                            quantile(CT, c(0.5, 0.025, 0.975)))
load("./2Y/2Y_1000_FULL.RData")
res_2Y_1000_FULL <- unname(rbind(cbind(apply(res_sim, 2, mean),
                                apply(res_sim, 2, function(x) quantile(x, c(0.025))),
                                apply(res_sim, 2, function(x) quantile(x, c(0.975)))),
                                quantile(DIC, c(0.5, 0.025, 0.975))-DIC_ref,
                                quantile(WAIC, c(0.5, 0.025, 0.975))-WAIC_ref,
                          quantile(CT, c(0.5, 0.025, 0.975))))

True_2Y <- c("$\\sigma_{c1}$=1 & ",
             "$\\sigma_{c2}$=0.2 & ",
             "$\\sigma_{c3}$=0.1 & ",
             "$\\sigma_{c4}$=0.5 & ",
             "$\\rho_1$=0.9 & ",
             "$\\rho_2$=0.9 & ",
             "$\\rho_3$=0.8 & ",
             "\\hdashline $\\beta_{01}$=0.2 & ",
             "$\\beta_{11}$=-0.1 & ",
             "$\\beta_{02}$=0.2 & ",
             "$\\beta_{12}$=-0.1 & ",
             "\\hdashline DIC & ",
             "WAIC & ",
             "Comp. time (s) & ")

res_2Y_30_L1F <- apply(round(res_2Y_30_L1, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_2Y_30_L3F <- apply(round(res_2Y_30_L3, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_2Y_30_L5F <- apply(round(res_2Y_30_L5, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_2Y_30_L3_MISF <- apply(round(res_2Y_30_L3_MIS, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_2Y_30_FULLF <- apply(round(res_2Y_30_FULL, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_2Y_1000_L1F <- apply(round(res_2Y_1000_L1, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_2Y_1000_L3F <- apply(round(res_2Y_1000_L3, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_2Y_1000_L5F <- apply(round(res_2Y_1000_L5, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_2Y_1000_L3_MISF <- apply(round(res_2Y_1000_L3_MIS, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))
res_2Y_1000_FULLF <- apply(round(res_2Y_1000_FULL, 2), 1, function(x) paste0(x[1], " [", x[2], ", ", x[3], "]"))

tab_2Y_30 <- cbind(paste0(True_2Y, res_2Y_30_L1F, " & ",res_2Y_30_L3F,
                          " & ", res_2Y_30_L5F, " & ", res_2Y_30_L3_MISF,
                          " & ", res_2Y_30_FULLF," & \\\\\n"))
cat(tab_2Y_30)

tab_2Y_1000 <- cbind(paste0(True_2Y, res_2Y_1000_L1F, " & ", res_2Y_1000_L3F,
                            " & ", res_2Y_1000_L5F, " & ", res_2Y_1000_L3_MISF, " & ",
                            res_2Y_1000_FULLF," & \\\\\n"))
cat(tab_2Y_1000)







