
#Combined program for running simulation and generating plots

#This program accepts a RANGE OF PRIOR distributions for delta by choice of delta.l, delta.u and tails. It then constructs the optimal CPZ decision rule and the optimal PPZ decision rule based on the specified cpmin=ppmin and specified cpmax=ppmax. It plots the decision rules and the CP and PP functions. Also, it simulates the designs to get unconditional power and unconditional sample size at different values of delta in the range of delta.vec. Finally, it also produces zonewise power and expected sample size.

#Footnotes:
#(1)Note that the optimal CPZ decision rule depends on the specified delta.min. But the optimal PPZ decision rule is free of dependence on delta.min. It does however depend on the choice of prior for delta. The prior is specified by three parameters, delta_l, delta_u and tails, such that P0(delta<delta_l)=P0(delta>delta_u)=tails. This specification automatically induces the prior mean to be m0=(delta_l+ldelta_u)/2 and the prior variance to be s0=(delta_u - delta_l)/(qnorm(1-tails)-qnorm(tails))

#(2)If the prior is informative (tails value is small, say 0.0001), then the posterior mean of delta will be more influenced by m0 and less influenced by the interim z1. Thus we would expect the PPZ decision rule to be very sensitive to the choice of m0. And this would also spill over into the unconditional operating characteristics of the PPZ design. For example if m0 is equal to delta.min, and (delta.l, delta.u) is a small interval with m0 at its center, and tails is 0.0001, the PPZ design should closely resemble the CPZ design with a small loss of power to account for having spread the smallest clinically meaningful delta across (delta.l, delta.min) instead of placing a point mass at delta.min.

#(3)But if the prior on delta is non-informative, (tails value is large, say 0.45, or alternatively delta_l and delta_u are very far apart, say 6*sigma apart) then the posterior mean will be more influenced by z1 and less by m0. Thus we would expect the PPZ decision rule to be much less sensitive to the choice of m0 and we would expect it to lose more unconditional power because the design is not being optimized at delta=delta.min but the uncertainty about delta.min is being spread out over the region (delta.l, delta.u)

#(4) In this formulation "tails" refers to the area to the left of delta_L OR to the right of delta_U but NOT to the sum of the areas

#(5) n.sim is the parameter for selecting the number of simulations. The plots and tables use n.sim=100000. But if experimenting with the code, use n.sim=1000 to save time
rm(list=ls())

## load Packages
library(tidyverse)

#--------------------------------------------------------------------------------------------------------------------------------

# To prevent R from creating an Rplots.pdf file when this script is executed non-interactively.
options(device = function(...) {
  grDevices::pdf(NULL)
})

## Paths and directories
# Sourcing the helper functions from the tools directory
tryCatch(
  suppressWarnings(source(file.path("tools", "helper.R"), local = .GlobalEnv)),
  error = function(first.error) {
    source(
      file.path("..", "..", "tools", "helper.R"),
      local = .GlobalEnv
    )
  }
)

# Setting up the base directory for executing this code
base.dir <- get.base.dir(
  cur.dir = getwd(),
  sec.dir = "Chapter 5/Sec 5.3/"
)

prog.in <- file.path(base.dir, "cpz_ppz_functions.r")
out.dir <- base.dir

# root <- getwd()

# base.dir <- file.path(root, "Adaptive with Bayesian", "CPZ_PPZ_comparison")
# prog.in  <- file.path(base.dir, "Programs", "cpz_ppz_functions.r")
# out.dir  <-file.path(base.dir, "Results.Ch4") # file.path(base.dir, "Results")


# base.dir <- file.path("D:\\Cyrus-Office\\Papers\\Adaptive\\MehtaPapers\\Optimal Promising Zone Designs\\Software\\CPZ_PPZ_comparison")
# prog.in <- file.path(base.dir, "Programs", "cpz_ppz_functions.R")
# out.dir <- file.path(base.dir, "Results")

#--------------------------------------------------------------------------------------------------------------------------------

## Saving Results
Simulation_by_delta_CPZ<- data.frame("True_delta"= numeric(0), "Average_SS_CPZ" =numeric(0),
                                     "Power_CPZ" = numeric(0), "Min_CP_PP" = numeric(0),"Max_CP_PP" = numeric(0))

Simulation_by_delta_PPZ<- data.frame("SimID"=numeric(0), "delta.l"= numeric(0), "delta.u"= numeric(0),
                                     "tails"= numeric(0), "True_delta"= numeric(0), "Average_SS_PPZ" = numeric(0),
                                     "Power_PPZ" = numeric(0), "Min_CP_PP" = numeric(0),"Max_CP_PP" = numeric(0))


Simulation_Zonewise_CPZ<-data.frame("Min_CP_PP" = numeric(0),"Max_CP_PP" = numeric(0),
                                    "True_delta"= numeric(0),
                                    "Zone"= character(0), "Prob_Zonewise_CPZ"=numeric(0),
                                    "Average_SS_Zonewise_CPZ" = numeric(0), "Power_Zonewise_CPZ" = numeric(0))


Simulation_Zonewise_PPZ<-data.frame("SimID"=numeric(0),"Min_CP_PP" = numeric(0),"Max_CP_PP" = numeric(0),
                                    "delta.l"= numeric(0), "delta.u"= numeric(0),
                                    "tails"= numeric(0), "True_delta"= numeric(0),
                                    "Zone"= character(0), "Prob_Zonewise_PPZ"=numeric(0),
                                    "Average_SS_Zonewise_PPZ" = numeric(0), "Power_Zonewise_PPZ" = numeric(0))

ppz.boundary   <-data.frame("SimID"=numeric(0),
                            "delta.l"= numeric(0), "delta.u"= numeric(0),
                            "tails"= numeric(0), "cf" = numeric(0),
                            "z.l" = numeric(0), "z.c" = numeric(0),
                            "z.u" = numeric(0), "ce" = numeric(0))

#--------------------------------------------------------------------------------------------------------------------------------

## Input files
source(prog.in)

#--------------------------------------------------------------------------------------------------------------------------------

## Specifying design parameters
n1 <- 85
n2 <- 170
k <- 2
n2min <- n2                                      # minimum final sample size after SSR
n2max <- n2*k                                    # maximum final sample size after SSR
w1 <- sqrt(n1/n2)                                # stage 1 weight for CHW combination statistic
w2 <- sqrt((n2 - n1)/n2)                         # stage 2 weight for CHW combination statistic

z.alpha <- 1.96859564                            # critical boundary at second look
ce <- 2.96258805                                 # early efficacy boundary 
cf <- -1                                         # early futility boundary

cpmin <-ppmin <- 0.8
cpmax <- ppmax <- 0.9

# simulations parameters
#delta.vec <- seq(0, 10, 1)
delta.vec <- seq(7, 10, 1)
delta.min <-7
n.sim <- 100000

# parameters for prior density of delta
sigma <- 20
delta.l <-  (delta.min) - c(1, 7, 6*sigma) #c(1, 3, 5, 7, sigma, 3*sigma, 6*sigma) #delta.min - c(1, 3, 5, 7, sigma, 3*sigma, 6*sigma) 
delta.u <-  (delta.min) + c(1, 7, 6*sigma) #c(1, 3, 5, 7, sigma, 3*sigma, 6*sigma) #delta.min + c(1, 3, 5, 7, sigma, 3*sigma, 6*sigma)

tails <- 0.0001
s0 <- rep(0, length(delta.l))
m0 <- rep(0, length(delta.l))
lb <- rep(0, length(delta.l))
ub <- rep(0, length(delta.l))

######----------------------------------------------------------------------------------------------------------

########################################################################

### CPZ Simulations (100000 sims only) (currently changed to 1000 to save time)

CPZ.Sim.out<- OC.Power.SS.deltas.CPZ(delta.vec, delta=delta.min,
                                     n1=n1, n2=n2, k=k, ce=ce, cf=cf, z.alpha=z.alpha, 
                                     cpmin=cpmin, cpmax=cpmax, n.sim=n.sim, seed=NULL)

## CPZ design boundaries
cpz.boundry <- CPZ.Sim.out$boundaries

## CPZ Summary across deltas
CPZ_OC_out <- CPZ.Sim.out$delta_wise_OC
Simulation_by_delta_CPZ_temp <-cbind(CPZ_OC_out[,-4], cpmin, cpmax)
names(Simulation_by_delta_CPZ_temp) <- names(Simulation_by_delta_CPZ)
Simulation_by_delta_CPZ <- rbind(Simulation_by_delta_CPZ,Simulation_by_delta_CPZ_temp)


## CPZ Zonewise Summaries
delta.true <-7   ### delta to obtain zonewise summaries, i.e, response generation delta
CPZ.SSR.out <-SSR.Power.Zonewise.CPZ(cpmin=cpmin, cpmax=cpmax, n1=n1, n2=n2, k=k,
                                     ce=ce, cf=cf, z.alpha=z.alpha,
                                     delta=delta.min, delta.true=delta.true,
                                     sigma, n.sim=n.sim, seed=NULL)

ssr.cpz <- CPZ.SSR.out$ssr.power.df     ## SSR-CP data.frame not saved to csv due to the size

ssr.cpz.zonewise <- cbind(cpmin, cpmax, CPZ.SSR.out$zonewise_summaries[,-6]) 
names(ssr.cpz.zonewise) <- names(Simulation_Zonewise_CPZ)
Simulation_Zonewise_CPZ <- rbind(Simulation_Zonewise_CPZ, ssr.cpz.zonewise)

## Saving all CPZ summaries to csv files
# write.table(cpz.boundry,
#             file.path(out.dir, "CPZ Summaries", "CPZ_Boundaries.csv"),
#             sep=",", row.names=FALSE)

# write.table(Simulation_by_delta_CPZ,
#             file.path(out.dir, "CPZ Summaries","CPZ_Operating_Characteristics_across_deltas.csv"),
#             sep=",", row.names=FALSE)

# write.table(Simulation_Zonewise_CPZ,
#             file.path(out.dir, "CPZ Summaries", paste0("CPZ_Operating_Characteristics_Zonewise_at_delta= ",
#                                                        delta.true, ".csv")),
#             sep=",", row.names=FALSE)

write.table(cpz.boundry,
            file.path(out.dir, "CPZ_Boundaries.csv"),
            sep=",", row.names=FALSE)

write.table(Simulation_by_delta_CPZ,
            file.path(out.dir, "CPZ_Operating_Characteristics_across_deltas.csv"),
            sep=",", row.names=FALSE)

write.table(Simulation_Zonewise_CPZ,
            file.path(out.dir, paste0("CPZ_Operating_Characteristics_Zonewise_at_delta= ",
                                                       delta.true, ".csv")),
            sep=",", row.names=FALSE)

## CPZ Data for plotting SSR and CP
z1vec<- seq(-2,4, 0.01)

z.l.cpz <- get.z1.cp(n1, n2, n2x=n2max, z.alpha, delta=delta.min, sigma, target.cp=cpmin, ce, cf) 
z.c.cpz <- get.z1.cp(n1, n2, n2x=n2max, z.alpha, delta=delta.min, sigma, target.cp=cpmax, ce, cf)
z.u.cpz <- get.z1.cp(n1, n2, n2x=n2, z.alpha, delta=delta.min, sigma, target.cp=cpmax, ce, cf)

ss.cpz.out<-sapply(z1vec, ss.cpz.function, zl=z.l.cpz, zc=z.c.cpz, zu=z.u.cpz,
                   ce=ce, cf=cf, z.alpha=z.alpha, target.cp=cpmax,
                   n1=n1, sigma=sigma, delta=delta.min, n2=n2, n2x=n2max)

cp.cpz.out<- mapply(FUN = cp.star, z1=z1vec,  n2 = n2, n2star=ss.cpz.out, n1=n1, 
                    delta=delta.min, sigma=sigma, z.alpha=z.alpha, ce=ce, cf=cf)

#####################################################################################################


start.time <- proc.time()
for(i in 1:length(delta.l)){
 
   ## PPZ Simulations  
  s0[i] <- (delta.u[i] - delta.l[i])/(qnorm(1-tails)-qnorm(tails))
  m0[i] <- delta.l[i]-s0[i]*qnorm(tails)
  n0 <- (sigma/s0[i])^2

# Make n.sim=1000 for experimenting. Make n.sim=100000 for final  
  PPZ.Sim.out<- OC.Power.SS.deltas.PPZ(ppmin=ppmin, ppmax=ppmax, m0=m0[i], s0=s0[i],
                                       n1=n1, n2=n2, k=k, ce=ce, cf=cf, z.alpha=z.alpha,
                                       delta.vec, sigma,
                                       n.sim=n.sim, seed=NULL)
  
  
  ## PPZ design boundaries
  ppz.bdry.temp <-cbind(i, delta.l[i], delta.u[i], tails, PPZ.Sim.out$boundaries)
  
  names(ppz.bdry.temp) <- names(ppz.boundary)
  ppz.boundary <- rbind(ppz.boundary, ppz.bdry.temp)

  ## PPZ Summary across deltas
  PPZ_OC_out <- PPZ.Sim.out$delta_wise_OC
  Simulation_by_delta_PPZ_temp<-cbind(i, delta.l[i], delta.u[i], tails, PPZ_OC_out[, -4], ppmin, ppmax)
  names(Simulation_by_delta_PPZ_temp)<-names(Simulation_by_delta_PPZ)
  Simulation_by_delta_PPZ <- rbind(Simulation_by_delta_PPZ, Simulation_by_delta_PPZ_temp)
  
  
  ## PPZ Zonewise Summary
  delta.true <-7   ### delta to obtain zonewise summaries, i.e, response generation delta
  PPZ.SSR.out <-SSR.Power.Zonewise.PPZ(ppmin=ppmin, ppmax=ppmax, m0=m0[i], s0= s0[i],
                                       n1=n1, n2=n2, k=k, ce=ce, cf=cf, z.alpha=z.alpha,
                                       delta.true=delta.true, sigma,
                                       n.sim=n.sim, seed=NULL)
  
  ssr.ppz.zonewise <- cbind(i, ppmin, ppmax, delta.l[i],
                            delta.u[i], tails, 
                            PPZ.SSR.out$zonewise_summaries[,-6]) 
  names(ssr.ppz.zonewise) <- names(Simulation_Zonewise_PPZ)
  Simulation_Zonewise_PPZ<- rbind(Simulation_Zonewise_PPZ, ssr.ppz.zonewise )
  
  
  ssr.ppz <- PPZ.SSR.out$ssr.power.df ## SSR-PP data.frame not saved to csv due to the size
  
  
  
  
  # Plotting Operating characteristics across different values of delta
 
  lb[i]<- delta.l[i]
  ub[i]<- delta.u[i]
  deltas <- unique(PPZ_OC_out$True_delta)
  xlim.l <- floor(min(PPZ_OC_out$Average_SS, CPZ_OC_out$Average_SS)/20)*20
  xlim.u <- ceiling(max(PPZ_OC_out$Average_SS, CPZ_OC_out$Average_SS)/20)*20
  std.delta  <- deltas/sigma
    
  
  obj = list(dl=lb[i], du=ub[i])
  
#*******************************************************************************************************
##BLOCK 1. Plotting Expected Sample Size and Unconditional Power vs delta/sigma 
  
# png(filename = file.path(out.dir, "ESS_Power_vs_delta",
#                            paste0("ESS_Power_CPZ_PPZ_delta.u=", ub[i],
#                                   "_delta.l=", lb[i], "_tails=", tails,
#                                   "_cp_pp_min=", cpmin, "_cp_pp_max=", cpmax, ".png")),
#       width = 1000, height = 600, units = "px")
png(filename = file.path(out.dir,
                           paste0("ESS_Power_CPZ_PPZ_delta.u=", ub[i],
                                  "_delta.l=", lb[i], "_tails=", tails,
                                  "_cp_pp_min=", cpmin, "_cp_pp_max=", cpmax, ".png")),
      width = 1000, height = 600, units = "px")
  
  par(mfrow=c(1,2))
  
  ## Expected Sample Size vs delta plot 
  plot(1, type = "n",
       xlim = c(min(std.delta), max(std.delta)), ylim = c(xlim.l, xlim.u),
       xlab = expression(bold(paste(delta,"/", sigma))), 
       ylab = "Expected Sample Size",
       main = bquote(bold("Expected Sample Size at "~
                            delta[l]== .(obj$dl)~";"~ delta[u]== .(obj$du))),
       
       sub = paste0 ("CP Min= PP Min= ", cpmin,
                     " and CP Max= PP Max= ", cpmax),
       font.main=2, font.lab=2, font.axis=1, cex.main=1.4, cex.lab=1.2, cex.axis=1.2)
  
  abline(h=seq(xlim.l,xlim.u,10), v=seq(min(std.delta), max(std.delta), 0.05),
         col = "gray90", lwd = 2, lty = 1)
  lines(CPZ_OC_out$True_delta/sigma, CPZ_OC_out$Average_SS, type = "l",
        col = "red", lwd = 3, lty = 1)
  lines(PPZ_OC_out$True_delta/sigma, PPZ_OC_out$Average_SS, type = "l",
        col = "black", lwd = 2, lty = 2)
  legend("topleft", c("PPZ","CPZ"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)
  
  #........................................................................................
  
  ## Unconditional Power vs delta plot
  plot(1, type = "n",
       xlim = c(min(std.delta), max(std.delta)), ylim = c(0, 1),
       xlab = expression(bold(paste(delta,"/", sigma))), 
       ylab = "Power",
       main = bquote(bold("Unconditional Power at "~
                            delta[l]== .(obj$dl)~ ";"~ delta[u]== .(obj$du))),
       sub = paste0 ("CP Min= PP Min= ", cpmin,
                     " and CP Max= PP Max= ", cpmax),
       font.main=2, font.lab=2, font.axis=1, cex.main=1.4, cex.lab=1.2, cex.axis=1.2 )
  
  # Conditional power without SSR
  abline(h=seq(0,1,0.1), v=seq(min(std.delta), max(std.delta), 0.05),
         col = "gray90", lwd = 2, lty = 1)
  lines(CPZ_OC_out$True_delta/sigma, CPZ_OC_out$Power/100, type = "l",
        col = "red", lwd = 3, lty = 1)
  lines(PPZ_OC_out$True_delta/sigma, PPZ_OC_out$Power/100, type = "l",
        col = "black", lwd = 2, lty = 2)
  legend("topleft", c("PPZ","CPZ"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)
  
  dev.off()

#######***************************************************************************************************
  ## PPZ Data for plotting SSR and CP
  z.l.ppz <- get.z1.pp(target.pp=ppmin, m0=m0[i], s0= s0[i], n1=n1, sigma=sigma, 
                       n2=n2, n2x=n2max, z.alpha=z.alpha, ce=ce, cf=cf)
  z.u.ppz <- get.z1.pp(target.pp=ppmax, m0=m0[i], s0= s0[i], n1=n1, sigma=sigma,
                       n2=n2, n2x=n2, z.alpha=z.alpha, ce=ce, cf=cf)
  z.c.ppz <- get.z1.pp(target.pp=ppmax, m0=m0[i], s0= s0[i], n1=n1, sigma=sigma, 
                       n2=n2, n2x=n2max, z.alpha=z.alpha, ce=ce, cf=cf)
  
  ss.ppz.out<-sapply(z1vec, ss.ppz.function, zl=z.l.ppz, zc=z.c.ppz, zu=z.u.ppz, 
                     ce=ce, cf=cf, z.alpha=z.alpha, target.pp=ppmax, m0=m0[i], s0= s0[i],
                     n1=n1,  n2=n2, n2x=n2max, sigma=sigma)
  
  pp.ppz.out<- mapply(FUN = pp, z1=z1vec, n2x=ss.ppz.out, m0=m0[i], s0= s0[i],
                      n1=n1, sigma=sigma, n2=n2, 
                      z.alpha=z.alpha, ce=ce, cf=cf)
  
  xlim.l1 <- min(z1vec)
  xlim.u1 <- max(z1vec)
  
  ylim.l1<- floor(min(ss.ppz.out, ss.cpz.out)/20)*20
  ylim.u1<- ceiling(max(ss.ppz.out, ss.cpz.out)/20)*20
  
##***************************************************************************************
#   
# ### BLOCK 2a.Plotting Sample Size Rule and Conditional Power vs Z1 for PPZ and CPZ (not needed for book. 2x2 grid plots suffice)
#   
#   png(filename = file.path(out.dir,"SSR and CP Plots",
#                            paste0("SSR_CP.PP_CPZ_PPZ_delta.u=", delta.u[i],
#                                   "_delta.l=", delta.l[i], "_tails=", tails,
#                                   "_cp_pp_min=", cpmin, "_cp_pp_max=", cpmax,".png")),
#       width = 1000, height = 600, units = "px")
#   
#   par(mfrow=c(1,2))
#   
#   # Sample Size Rule
#   plot(1, type = "n",
#        xlim = c(xlim.l1, xlim.u1), ylim = c(ylim.l1, ylim.u1),
#        xlab = "Z-Statistic at Interim Analysis", ylab = "Sample Size Rule",
#        main = bquote(bold("Sample Size Rule at: "~ delta~"/"~sigma==~0.35 ~";"~
#                             delta[l]== .(obj$dl)~";"~ delta[u]== .(obj$du))),
#        sub = paste0("PP Min= ",cpmin, " and PP Max= ", cpmax),
#        #sub = bquote("For Prior P(" ~ .(obj$dl) ~ " < " ~ delta ~ " < " ~ .(obj$du) ~ ") = 0.9999"),
#        font.main=2, font.lab=2, font.axis=1, cex.main=1.4, cex.lab=1.2, cex.axis=1.2)
#   
#   abline(h=seq(ylim.l1,ylim.u1, 20), v=seq(xlim.l1, xlim.u1, 0.5), col = "gray90", lwd = 2, lty = 1)
#   lines(z1vec, ss.cpz.out, type = "l", col = "red", lwd = 3, lty = 1)
#   lines(z1vec, ss.ppz.out, type = "l", col = "black", lwd = 2, lty = 2)
#   legend("topleft", c("PPZ","CPZ"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)
#   
#   #........................................................................................
#   
#   # Conditional power without SSR
#   
#   plot(1, type = "n",
#        xlim = c(xlim.l1, xlim.u1), ylim = c(0, 1),
#        xlab = "Z-Statistic at Interim Analysis", ylab = "Conditional Power or Predictive Power",
#        main = bquote(bold("Conditional Power at: "~ delta~"/"~sigma==~ 0.35 ~";"~
#                             delta[l]== .(obj$dl)~";"~ delta[u]== .(obj$du))),
#        sub = paste0("CP Min= PP Min= ",cpmin, " and CP Max= PP Max= ", cpmax),
#        #sub = bquote("For Prior P(" ~ .(obj$dl) ~ " < " ~ delta ~ " < " ~ .(obj$du) ~ ") = 0.9999"),
#        font.main=2, font.lab=2, font.axis=1, cex.main=1.4, cex.lab=1.2, cex.axis=1.2)
#  
#   abline(h=seq(0,1,0.1), v=seq(xlim.l1, xlim.u1, 0.5), col = "gray90", lwd = 2, lty = 1)
#   lines(z1vec, cp.cpz.out, type = "l", col = "red", lwd = 3, lty = 1)
#   lines(z1vec, pp.ppz.out, type = "l", col = "black", lwd = 2, lty = 2)
#   legend("topleft", c("PPZ","CPZ"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=0.5)
#   
#   dev.off()


#######***********************************************************************************
  ### BLOCK 2b.Plotting Sample Size Rule and Conditional Power vs Z1 for PPZ only (needed in book chapter)
i=2
  
#   png(filename = file.path(out.dir,"SSR and CP Plots",
#                            paste0("SSR_CP.PP_CPZ_PPZ_delta.u=", delta.u[i],
#                                   "_delta.l=", delta.l[i], "_tails=", tails,
#                                   "_cp_pp_min=", cpmin, "_cp_pp_max=", cpmax,".png")),
#       width = 1000, height = 600, units = "px")
  png(filename = file.path(out.dir,
                           paste0("SSR_CP.PP_CPZ_PPZ_delta.u=", delta.u[i],
                                  "_delta.l=", delta.l[i], "_tails=", tails,
                                  "_cp_pp_min=", cpmin, "_cp_pp_max=", cpmax,".png")),
      width = 1000, height = 600, units = "px")
  
  par(mfrow=c(1,2), mar = c(5, 6.5, 4, 2) + 0.1)
  
  # Predictive Power
  
  plot(1, type = "n",
       xlim = c(xlim.l1, xlim.u1), ylim = c(0, 1),
       xlab = "Z-Statistic at Interim Analysis", 
       ylab = bquote(bold("Predictive Power if " ~ delta == .(delta.true))),
       main = "Predictive Power",
       #sub = bquote("For Prior P(" ~ .(obj$dl) ~ " < " ~ delta ~ " < " ~ .(obj$du) ~ ") = 0.9999"),
       font.main=2, font.lab=2, font.axis=1, cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
  
  abline(h=seq(0,1,0.1), v=seq(xlim.l1, xlim.u1, 0.5), col = "gray90", lwd = 2, lty = 1)
  
  abline(v = z.l.ppz, lty = 3, lwd = 2, col = "red")
  abline(v = z.u.ppz, lty = 3, lwd = 2, col = "red")
  
  lines(z1vec, pp.ppz.out, type = "l", col = "black", lwd = 2, lty = 1)
  
  x.mid <- (z.l.ppz + z.u.ppz)/2
  
  text(x.mid+0.6, 0.15,
       labels = "Promising Zone",
       col = "black", cex = 1.8, font = 2)
  
  text(x.mid+0.4, 0.09,
       labels = paste0("Start = ", round(z.l.ppz, 2),
                       "   End = ", round(z.u.ppz, 2)),
       col = "black", cex = 1.75, font = 1)
  
  
  # Sample Size Rule
  
  plot(1, type = "n",
       xlim = c(xlim.l1, xlim.u1), ylim = c(ylim.l1, ylim.u1),
       xlab = "Z-Statistic at Interim Analysis", ylab = "Sample Size Rule",
       main = "Sample Size Rule",
       #sub = bquote("For Prior P(" ~ .(obj$dl) ~ " < " ~ delta ~ " < " ~ .(obj$du) ~ ") = 0.9999"),
       font.main=2, font.lab=2, font.axis=1, cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
  
  abline(h = seq(50, 1000, 50), v = seq(xlim.l1, xlim.u1, 0.5), col = "gray90", lwd = 2, lty = 1)
  
  abline(v = z.l.ppz, lty = 3, lwd = 2, col = "red")
  abline(v = z.u.ppz, lty = 3, lwd = 2, col = "red")
  
  lines(z1vec, ss.ppz.out, type = "l", col = "black", lwd = 2, lty = 1)
  
  text(x.mid-1, 120,
       labels = "Promising Zone",
       col = "black", cex = 1.8, font = 2)
  
  text(x.mid-1, 105,
       labels = paste0("Start = ", round(z.l.ppz, 2),
                       "   End = ", round(z.u.ppz, 2)),
       col = "black", cex = 1.8, font = 1)
  
  dev.off()
  

  

#########################################################################################################
##### BLOCK 3.  2x2 Grid Plots (png plots)
  
  xlim.l2 <- floor(min(PPZ_OC_out$Average_SS, CPZ_OC_out$Average_SS)/20)*20
  xlim.u2 <- ceiling(max(PPZ_OC_out$Average_SS, CPZ_OC_out$Average_SS)/20)*20

#   png(filename = file.path(out.dir, "2_by_2_grid",
#                            paste0("ESS_Power_SSR_CP_CPZ_PPZ_delta.u=", ub[i],
#                                   "_delta.l=", lb[i], "_tails=", tails,
#                                   "_cp_pp_min=", cpmin, "_cp_pp_max=", cpmax, ".png")),
#       width = 1000, height = 1000, units = "px")
  png(filename = file.path(out.dir,
                           paste0("ESS_Power_SSR_CP_CPZ_PPZ_delta.u=", ub[i],
                                  "_delta.l=", lb[i], "_tails=", tails,
                                  "_cp_pp_min=", cpmin, "_cp_pp_max=", cpmax, ".png")),
      width = 1000, height = 1000, units = "px")

  par(mfrow=c(2,2))

  ## Unconditional Power vs delta plot
  plot(1, type = "n",
       xlim = c(min(std.delta), max(std.delta)), ylim = c(0, 1),
       xlab = expression(bold(paste("Effect Size= ", delta,"/", sigma))),
       ylab = "Power",
       main = bquote(bold("Power if Prior P("~ .(obj$dl)~"< " ~ delta ~
                            "< " ~.(obj$du)~") = 0.9999 ")),
       sub = paste0 ("CP Min= PP Min= ", cpmin,
                     " and CP Max= PP Max= ", cpmax),
       font.main=2, font.lab=2, font.axis=1, cex.main=1.4, cex.lab=1.2, cex.axis=1.2 )

  abline(h=seq(0,1,0.1), v=seq(min(std.delta), max(std.delta), 0.025),
         col = "gray90", lwd = 2, lty = 1)
  lines(CPZ_OC_out$True_delta/sigma, CPZ_OC_out$Power/100, type = "l",
        col = "red", lwd = 3, lty = 1)
  lines(PPZ_OC_out$True_delta/sigma, PPZ_OC_out$Power/100, type = "l",
        col = "black", lwd = 2, lty = 2)
  legend("topleft", c("PPZ","CPZ"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)

  #........................................................................................

  ## Expected Sample Size vs delta plot
  plot(1, type = "n",
       xlim = c(min(std.delta), max(std.delta)), ylim = c(xlim.l2, xlim.u2),
       xlab = expression(bold(paste("Effect Size= ", delta,"/", sigma))),
       ylab = "Expected Sample Size",
       main = bquote(bold("Expected Sample Size if Prior P("~ .(obj$dl)~"< " ~ delta ~
                            "< " ~.(obj$du)~") = 0.9999 ")),
       sub = paste0 ("CP Min= PP Min= ", cpmin,
                     " and CP Max= PP Max= ", cpmax),
      font.main=2, font.lab=2, font.axis=1, cex.main=1.4, cex.lab=1.2, cex.axis=1.2)

  abline(h=seq(xlim.l2,xlim.u2,10), v=seq(min(std.delta), max(std.delta), 0.025), col = "gray90", lwd = 2, lty = 1)
  lines(CPZ_OC_out$True_delta/sigma, CPZ_OC_out$Average_SS, type = "l",
        col = "red", lwd = 3, lty = 1)
  lines(PPZ_OC_out$True_delta/sigma, PPZ_OC_out$Average_SS, type = "l",
        col = "black", lwd = 2, lty = 2)
  legend("topleft", c("PPZ","CPZ"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)


  ## SSR Plot
  xlim.l3 <- min(z1vec)
  xlim.u3 <- max(z1vec)

  ylim.l3<- floor(min(ss.ppz.out, ss.cpz.out)/20)*20
  ylim.u3<- ceiling(max(ss.ppz.out, ss.cpz.out)/20)*20

  plot(1, type = "n",
       xlim = c(xlim.l3, xlim.u3), ylim = c(0, 1),
       xlab = "Z-Statistic at Interim Analysis", 
       ylab = bquote("Pred/Cond Power if " ~ delta == .(delta.true)),
       main = bquote(bold("Conditional Power if Prior P("~ .(obj$dl)~"< " ~ delta ~
                            "< " ~.(obj$du)~") = 0.9999 ")),
       sub = paste0("CP Min= PP Min= ",cpmin, " and CP Max= PP Max= ", cpmax),

       font.main=2, font.lab=2, font.axis=1, cex.main=1.4, cex.lab=1.2, cex.axis=1.2)

  abline(h=seq(0,1,0.1), v=seq(xlim.l3, xlim.u3, 0.5), col = "gray90", lwd = 2, lty = 1)
  lines(z1vec, cp.cpz.out, type = "l", col = "red", lwd = 3, lty = 1)
  lines(z1vec, pp.ppz.out, type = "l", col = "black", lwd = 2, lty = 2)
  legend("topleft", c("PPZ","CPZ"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)

  #........................................................................................

  # Conditional and Predictive power plot
  plot(1, type = "n",
       xlim = c(xlim.l3, xlim.u3), ylim = c(ylim.l3, ylim.u3),
       xlab = "Z-Statistic at Interim Analysis", ylab = "Sample Size Rule",
       main = bquote(bold("Sample Size Rule if Prior P("~ .(obj$dl)~"< " ~ delta ~
                            "< " ~.(obj$du)~") = 0.9999 ")),
       sub = paste0("CP Min= PP Min= ",cpmin, " and CP Max= PP Max= ", cpmax),
       # main = expression(bold(paste("Sample Size Rule Prior P(", 6,"<", delta, "<",8,")= 0.9999"))),
       font.main=2, font.lab=2, font.axis=1, cex.main=1.4, cex.lab=1.2, cex.axis=1.2)

  # Sample Size Rule
  abline(h=seq(ylim.l3,ylim.u3, 20), v=seq(xlim.l3, xlim.u3, 0.5), col = "gray90", lwd = 2, lty = 1)
  lines(z1vec, ss.cpz.out, type = "l", col = "red", lwd = 3, lty = 1)
  lines(z1vec, ss.ppz.out, type = "l", col = "black", lwd = 2, lty = 2)
  legend("topleft", c("PPZ","CPZ"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)

  dev.off()

###############################################################################################  

###### BLOCK 4.  2x2 Grid Plots (PDF plots)
  
  xlim.l2 <- floor(min(PPZ_OC_out$Average_SS, CPZ_OC_out$Average_SS) / 20) * 20
  xlim.u2 <- ceiling(max(PPZ_OC_out$Average_SS, CPZ_OC_out$Average_SS) / 20) * 20
  
  xlim.l3 <- min(z1vec)
  xlim.u3 <- max(z1vec)
  
  ylim.l3 <- floor(min(ss.ppz.out, ss.cpz.out) / 20) * 20
  ylim.u3 <- ceiling(max(ss.ppz.out, ss.cpz.out) / 20) * 20
  
  # Revised x-axis limits for top row
  x.min <- 0.35
  x.max <- 0.50
  
#   plot.dir <- file.path(out.dir, "2_by_2_grid")
  plot.dir <- out.dir
  dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
  
  outfile <- file.path(
    plot.dir,
    paste0("ESS_Power_SSR_CP_CPZ_PPZ_delta.u=", ub[i],
           "_delta.l=", lb[i],
           "_tails=", tails,
           "_cp_pp_min=", cpmin,
           "_cp_pp_max=", cpmax,
           ".pdf")
  )
  
  pdf(file = outfile, width = 5, height = 7)
  
  par(mfrow = c(2, 2),
      mar = c(4.2, 4.2, 3.0, 1.2),
      mgp = c(2.2, 0.7, 0),
      tcl = -0.25,
      oma = c(0.2, 0.2, 0.2, 0.2))
  
  # ----------------------------
  # Common style settings
  # ----------------------------
  main.cex   <- 1.00
  sub.cex    <- 0.75
  lab.cex    <- 0.95
  axis.cex   <- 0.82
  legend.cex <- 0.78
  
  grid.col <- "gray85"
  grid.lwd <- 0.7
  ppz.lwd  <- 1.
  cpz.lwd  <- 1.
  axis.lwd <- 0.8
  box.lwd  <- 0.8
  
  prior.sub <- bquote("For Prior P(" ~ .(obj$dl) ~ " < " ~ delta ~ " < " ~ .(obj$du) ~ ") = 0.9999")
  
  #*******************************************************************************************************
  # 1. Power
  
  plot(1, type = "n",
       xlim = c(x.min, x.max),
       xaxt = "n",   # <-- suppress default x-axis
       ylim = c(0, 1),
       xlab = expression(bold("Effect Size = " ~ delta/sigma)),
       ylab = "Unconditional Power",
       main = "Unconditional Power",
       font.main = 2, font.lab = 2, font.axis = 1,
       cex.main = main.cex, cex.lab = lab.cex, cex.axis = axis.cex,
       lwd = axis.lwd)
  
  title(sub = prior.sub, cex.sub = sub.cex)
  
  abline(h = seq(0, 1, 0.1),
         v = seq(x.min, x.max, 0.025),
         col = grid.col, lwd = grid.lwd, lty = 1)
  
  lines(CPZ_OC_out$True_delta / sigma, CPZ_OC_out$Power / 100,
        type = "l", col = "red", lwd = cpz.lwd, lty = 1)
  lines(PPZ_OC_out$True_delta / sigma, PPZ_OC_out$Power / 100,
        type = "l", col = "black", lwd = ppz.lwd, lty = 2)
  
  axis(1, at = seq(x.min, x.max, by = 0.025), cex.axis = axis.cex)
  
  legend("topleft", c("PPZ", "CPZ"),
         col = c("black", "red"),
         lty = c(2, 1),
         lwd = c(ppz.lwd, cpz.lwd),
         cex = legend.cex,
         bty = "n")
  
  box(lwd = box.lwd)
  
  #*******************************************************************************************************
  # 2. Expected Sample Size
  
  plot(1, type = "n",
       xlim = c(x.min, x.max),
       xaxt = "n",   # <-- suppress default x-axis
       ylim = c(xlim.l2, xlim.u2),
       xlab = expression(bold("Effect Size = " ~ delta/sigma)),
       ylab = "Expected Sample Size",
       main = "Expected Sample Size",
       font.main = 2, font.lab = 2, font.axis = 1,
       cex.main = main.cex, cex.lab = lab.cex, cex.axis = axis.cex,
       lwd = axis.lwd)
  
  title(sub = prior.sub, cex.sub = sub.cex)
  
  abline(h = seq(xlim.l2, xlim.u2, 10),
         v = seq(x.min, x.max, 0.025),
         col = grid.col, lwd = grid.lwd, lty = 1)
  
  lines(CPZ_OC_out$True_delta / sigma, CPZ_OC_out$Average_SS,
        type = "l", col = "red", lwd = cpz.lwd, lty = 1)
  lines(PPZ_OC_out$True_delta / sigma, PPZ_OC_out$Average_SS,
        type = "l", col = "black", lwd = ppz.lwd, lty = 2)
  
  axis(1, at = seq(x.min, x.max, by = 0.025), cex.axis = axis.cex)
  
  legend("topleft", c("PPZ", "CPZ"),
         col = c("black", "red"),
         lty = c(2, 1),
         lwd = c(ppz.lwd, cpz.lwd),
         cex = legend.cex,
         bty = "n")
  
  box(lwd = box.lwd)
  
  #*******************************************************************************************************
  # 3. Conditional/Predictive Power
  
  plot(1, type = "n",
       xlim = c(xlim.l3, xlim.u3),
       ylim = c(0, 1),
       xlab = "Z-Statistic at Interim Analysis",
       ylab = bquote("Cond/Pred Power if " ~ delta == .(delta.true)),
       main = "Cond or Pred Power",
       font.main = 2, font.lab = 2, font.axis = 1,
       cex.main = main.cex, cex.lab = lab.cex, cex.axis = axis.cex,
       lwd = axis.lwd)
  
  title(sub = prior.sub, cex.sub = sub.cex)
  
  abline(h = seq(0, 1, 0.1),
         v = seq(xlim.l3, xlim.u3, 0.5),
         col = grid.col, lwd = grid.lwd, lty = 1)
  
  lines(z1vec, cp.cpz.out, type = "l",
        col = "red", lwd = cpz.lwd, lty = 1)
  lines(z1vec, pp.ppz.out, type = "l",
        col = "black", lwd = ppz.lwd, lty = 2)
  
  legend("topleft", c("PPZ", "CPZ"),
         col = c("black", "red"),
         lty = c(2, 1),
         lwd = c(ppz.lwd, cpz.lwd),
         cex = legend.cex,
         bty = "n")
  
  box(lwd = box.lwd)
  
  #*******************************************************************************************************
  # 4. Sample Size Rule
  
  plot(1, type = "n",
       xlim = c(xlim.l3, xlim.u3),
       ylim = c(ylim.l3, ylim.u3),
       xlab = "Z-Statistic at Interim Analysis",
       ylab = "Sample Size Rule",
       main = "Sample Size Rule",
       font.main = 2, font.lab = 2, font.axis = 1,
       cex.main = main.cex, cex.lab = lab.cex, cex.axis = axis.cex,
       lwd = axis.lwd)
  
  title(sub = prior.sub, cex.sub = sub.cex)
  
  abline(h = seq(ylim.l3, ylim.u3, 20),
         v = seq(xlim.l3, xlim.u3, 0.5),
         col = grid.col, lwd = grid.lwd, lty = 1)
  
  lines(z1vec, ss.cpz.out, type = "l",
        col = "red", lwd = cpz.lwd, lty = 1)
  lines(z1vec, ss.ppz.out, type = "l",
        col = "black", lwd = ppz.lwd, lty = 2)
  
  legend("topleft", c("PPZ", "CPZ"),
         col = c("black", "red"),
         lty = c(2, 1),
         lwd = c(ppz.lwd, cpz.lwd),
         cex = legend.cex,
         bty = "n")
  
  box(lwd = box.lwd)
  
  dev.off()
}



#####################################################################################
end.time <- proc.time()
elpased.time <- end.time - start.time

# ####################################
# write.table(Simulation_by_delta_PPZ,
#             file.path(out.dir, "PPZ Summaries", "PPZ_Operating_Characteristics_across_deltas.csv"),
#             sep=",", row.names=FALSE)

# write.table(ppz.boundary,
#             file.path(out.dir, "PPZ Summaries", "PPZ_Boundaries.csv"),
#             sep=",", row.names=FALSE)

# write.table(Simulation_Zonewise_PPZ,
#             file.path(out.dir, "PPZ Summaries",paste0("PPZ_Operating_Characteristics_Zonewise_at_delta= ",
#                                                       delta.true, ".csv")),
#             sep="," , row.names=FALSE)

# ###################################


