
## Run this code AFTER you have run "Run_Unconstrained_JT_Pilz_Comparison.r"
## and created Delta_Power_ESS_JT.csv and Delta_Power_ESS_Pilz.csv files

## This code  plots the JT optimal design used in Ch4 of book
#Figure 1: SSR rule and CP of JT design that maximizes conditional power for Z1 in (-4, +4) with uncon power=0.75 and delta/sigma = 0.35
#Figure 2: Unconditional Power and ESS of these optimal JT and Pilz designs for delta/sigma in (0.35 to 0.5)
                                                                                  
## Additionally, this code produces the actual dataframe of z1vec versus SSR_JT used in the conditional power plot

## Footnote1: for future reviewers of this work: Figure 1 is constructed by plotting SSR_JT_vs_z1vec.csv. 
## That csv file shows the promising zone starting at z1=-0.364. This is 
## a slight approximation. If you do a more detailed study using the code in "marginal.Rev-Cost.r" you will see that
## in fact the promising zone starts at -0.6. But the impact is negligible. The discrepancy is a result of intergerization of the optimal sample size.

## Footnote2: Ignore everything that references the Pilz design including SS.df. 
## Thatdataframe  was created for matching power with JT. Here we are concerned only with a JT design having a desired expected sample size   
### However there are references to the csv file Delta_Power_ESS_Pilz.csv that was created initially when running
### Run_Unconstrained_JT_Pilz_Comparison.R. Although that file is no longer in use it has to be present in the out.dir
### As the ouput from some earlier (now unnecessary) run of Run_Unconstrained_JT_Pilz_Comparison.R. So make sure its there,


rm(list=ls())

library(tidyverse)

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

#--------------------------------------------------------------------------------------
# base.dir <- file.path("D:", "Apurva", "CIC", "CIC", "OPZ", "Software",
#                       "Optimizing_TwoStage_Design", "All Optimization Designs", 
#                        "June 2021", "Unconstrained-Pilz-JT-ComparisonPrograms_June2021")

# root <- getwd()
# base.dir<-file.path(root, "Adaptive with PILZ-JT comparison")
# prog.dir <- file.path(base.dir, "Programs")
# out.dir <- file.path(base.dir, "Results.cf.-1")
#out.dir <- file.path("D:/Cyrus-Office/BookProject/All-Chapters/screens")

# Setting up the base directory for executing this code
base.dir <- get.base.dir(
  cur.dir = getwd(),
  sec.dir = "Chapter 5/Sec 5.4/"
)

prog.dir <- base.dir
out.dir <- base.dir

prog.in1<- file.path(prog.dir, "optimal_pz_fns_v4(match power).R")
prog.in2<- file.path(prog.dir, "Optimization_functions_PilzDesigns.R")  
data.in1 <- file.path(out.dir, "Delta_Power_ESS_JT.csv")
data.in2 <- file.path(out.dir, "Delta_Power_ESS_Pilz.csv")


#--------------------------------------------------------------------------------------
source(prog.in1)
source(prog.in2)
JT.dat <- read.csv(data.in1)
Pilz.dat <- read.csv(data.in2)

## Design Inputs
sigma <- 1    

k <- 2
n1 <- 85                                         # stage 1 total sample size
n2 <- 170                                        # planned final sample size
n2min <- n2                                      # minimum final sample size after SSR
n2max <- n2*k                                  # maximum final sample size after SSR

w1 <- sqrt(n1/n2)                                # stage 1 weight for CHW combination statistic
w2 <- sqrt((n2 - n1)/n2)                         # stage 2 weight for CHW combination statistic

ce <- 2.96258805                                 # early efficacy boundary 
c2 <- 1.96859564                                 # efficacy boundary at second look
cf <- -1  

delta0 <- 0                                     # delta under null
delta.true <- seq(0.3, 0.8,0.05)                # true/simulation/response generation delta delta
delta.target <- delta.tilde <- 0.35             # target delta used for finding gamma/lambda

gamma<- unique(JT.dat$gamm)
lambda <- unique(Pilz.dat$lambda)
z1vec <- seq(-4, 4, 0.001) #seq(-4, 4, 0.001)

ESS_JT <- JT.dat$ESS_JT
Power_JT <- JT.dat$Power_JT
ESS_Pilz <- Pilz.dat$ESS_Pilz
Power_Pilz <- Pilz.dat$Power_Pilz
##-------------------------------------------------------------------------------------
## SSR and CP Plots

SSR_JT <- vapply(X = z1vec, FUN = n2x_jt, FUN.VALUE = numeric(1),
                 n1 = n1, n2 = n2, n2min = n2min, n2max = n2max,
                 bdry = c2, w1 = w1, w2 = w2,
                 theta_tilde = delta.tilde, sigma = sigma , 
                 ce = ce, cf = cf, gam = gamma,
                 revenue = revenue, cost = cost,
                 objective = objective, stepsize = 1)

# Conditional power of matching JT design
CP_JT <- mapply(FUN = cp,
                n2x = SSR_JT,
                z1 = z1vec,
                n1 = n1, bdry = c2, w1 = w1, w2 = w2,
                theta_cp = delta.tilde, sigma = sigma, ce = ce, cf= cf)


SS.df <- bind_rows(lapply(z1vec, get.n2star.stop, lambda=lambda,
                          delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                          sigma=sigma, c2=c2, ce=ce, cf=cf,
                          n1=n1, n2=n2, k=k))
SSR_Pilz <- SS.df$n2star.final
CP_Pilz <- SS.df$CondProb


ylim.l <- floor(min(c(SSR_Pilz, SSR_JT))/20)*20
ylim.u <- ceiling(max(c(SSR_Pilz, SSR_JT))/20)*20

## Specify limits and text for the Promising Zone
## You will get these from SSR_JT_vs_z1vec.csv
## Promising zone limits
z.pz.l <- -0.199
z.pz.u <-  2.294
pz.mid <- mean(c(z.pz.l, z.pz.u))
pz.lab <- "Promising Zone\n -0.199 <= Z1 <= 2.294"


# Start the Plots
# Operating Characteristics: Sample Size Rule and Conditional Power

png(filename = file.path(out.dir,"SSR_CP_vs_z1.png"),
    width = 1000, height = 600, units = "px")

par(mfrow=c(1,2))
cex.main <- 2
cex.lab  <- 1.6
cex.axis <- 1.6
lwd.curve <- 4
lwd.grid  <- 2

plot(1, type = "n",
     xlim = c(min(z1vec), max(z1vec)),
     ylim = c(ylim.l, ylim.u),
     xlab = "Interim Test Statistic Z1", ylab = "Sample Size Rule",
     main = expression(bold(paste("Sample Size Rule at ", delta[min], "= 7"))),
     font.main=2, font.lab=2, font.axis=1,
     cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis,
     xaxt = "n")

axis(1, at = seq(-4, 4, by = 1), cex.axis = cex.axis)

#  Sample Size Rule
# Get tick positions
yticks <- axTicks(2)
xticks <- axTicks(1)

# Midpoints between existing ticks
y_mid <- (yticks[-1] + yticks[-length(yticks)]) / 2
x_mid <- (xticks[-1] + xticks[-length(xticks)]) / 2

# Extra grid lines (in between)
abline(h = y_mid, v = x_mid, col = "gray90", lwd = 2)
abline(h = yticks,
       v = xticks,
       lwd = 2, lty = 1, col = "gray90")

# Vertical lines for the Promising Zone
abline(v = c(z.pz.l, z.pz.u),
       col = "black", lwd = 1.2, lty = 2)

# Now the SSR curve
lines(z1vec, SSR_JT, type = "l", col = "red", lwd = 3, lty = 1)
#lines(z1vec, SSR_Pilz, type = "l", col = "black", lwd = 3, lty = 2)
#legend("topleft", c("Pilz","JT"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)

# Finally, the text of the Promising Zone
text(x = pz.mid,
     y = ylim.l + 0.25 * (ylim.u - ylim.l),
     labels = pz.lab,
     cex = 1.4)

# Now start plotting
plot(1, type = "n",
     xlim = c(min(z1vec), max(z1vec)),
     ylim = c(0, 1),
     xlab = "Interim Test Statistic Z1", ylab = "Condiltional Power",
     main = expression(bold(paste("Conditional Power at ", delta[min], "= 7"))),
     font.main=2, font.lab=2, font.axis=1,
     cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis,
     xaxt = "n")

axis(1, at = seq(-4, 4, by = 1), cex.axis = cex.axis)

# Expected Sample Size
abline(h=seq(0,1, 0.1),
       v=seq(min(z1vec), max(z1vec), 1),
       lwd = 2, lty = 1, col = "gray90")


# Vertical lines for the Promising Zone
abline(v = c(z.pz.l, z.pz.u),
       col = "black", lwd = 1.2, lty = 2)

# Now the CP curve
lines(z1vec, CP_JT, type = "l", col = "red", lwd = 3, lty = 1)

# Finally the text for the promising zone
text(x = pz.mid+0.75,
     y = 0.18,
     labels = pz.lab,
     cex = 1.4)


#lines(z1vec, CP_Pilz, type = "l", col = "black", lwd = 3, lty = 2)
#legend("topleft", c("Pilz","JT"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)

dev.off()

#-------------------------------------------------------------------------------------
## Plots of ESS and Unconditional Power at delta/sigma between 0.3 and 0.8
xlim.l <- min(JT.dat$delta)
xlim.u <- max(JT.dat$delta)

ylim.l_1 <- floor(min(JT.dat$ESS_JT)/20)*20
ylim.u_1 <- ceiling(max(JT.dat$ESS_JT)/20)*20

ylim.l_2 <- floor(min(JT.dat$Power_JT)*10)/10
ylim.u_2 <- ceiling(max(JT.dat$Power_JT)*10)/10


# Operating Characteristics: Power and Expected Sample Size
png(filename = file.path(out.dir ,"Power_ESS_vs_delta.png"),
    width = 1000, height = 600, units = "px")

cex.main <- 2
cex.lab  <- 1.6
cex.axis <- 1.6
lwd.curve <- 4
lwd.grid  <- 2

par(mfrow=c(1,2))
plot(1, type = "n",
     xlim = c(xlim.l, xlim.u),
     ylim = c(ylim.l_1, ylim.u_1),
     xlab = expression(bold(paste("Effect Size ", delta,"/", sigma))), ylab = "Expected Sample Size",
     main = "Expected Sample Size Plot" ,
     font.main=2, font.lab=2, font.axis=1, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)

# Expected Sample Size
#Existing gridlines
abline(h=seq(ylim.l_1, ylim.u_1, 20),
       v=seq(xlim.l, xlim.u, 0.1),
       lwd = 2, lty = 1, col = "gray90")
# NEW midpoint grid lines
abline(h = seq(ylim.l_1 + 10, ylim.u_1 - 10, 20),
       v = seq(xlim.l + 0.05, xlim.u - 0.05, 0.1),
       lwd = 2, lty = 1, col = "gray90")

lines(delta.true, ESS_JT, type = "l", col = "red", lwd = 3, lty = 1)
#lines(delta.true, ESS_Pilz, type = "l", col = "black", lwd = 3, lty = 2)
#legend("topleft", c("Pilz","JT"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)


plot(1, type = "n",
     xlim = c(xlim.l, xlim.u),
     ylim = c(ylim.l_2, ylim.u_2),
     xlab = expression(bold(paste("Effect Size ", delta,"/", sigma))), ylab = "Power",
     main = "Unconditional Power Plot",
     font.main=2, font.lab=2, font.axis=1, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)

# Expected Sample Size
#Existing grid lines
abline(h=seq(ylim.l_2, ylim.u_2, 0.1),
       v=seq(xlim.l, xlim.u, 0.1),
       lwd = 2, lty = 1, col = "gray90")
# New mid-point grid lines
abline(h = seq(ylim.l_2 + 0.05, ylim.u_2 - 0.05, 0.1),
       v = seq(xlim.l + 0.05, xlim.u - 0.05, 0.1),
       lwd = 2, lty = 1, col = "gray90")

lines(delta.true, Power_JT, type = "l", col = "red", lwd = 3, lty = 1)
#lines(delta.true, Power_Pilz, type = "l", col = "black", lwd = 3, lty = 2)
#legend("topleft", c("Pilz","JT"), col=c("black","red"), lty = c(2,1), lwd = 2, cex=1)

dev.off()

###################################################################################################################
### Create a data frame with SSR versus z1
SSR_JT.df <- data.frame(
  z1vec  = z1vec,
  SSR_JT = SSR_JT
)

write.csv(
  SSR_JT.df,
  file = file.path(out.dir, "SSR_JT_vs_z1vec.csv"),
  row.names = FALSE
)

#####Diagnostics######################################################################################################
## This code will print the integer values of n2x and obj function for a fixed z1 and also print the best line.
z1.fixed <- -0.199

x <- seq(n2, n2max+200, by = 1)

diagnostic.df <- data.frame(
  n2x = x,
  CP = vapply(x, cp, numeric(1),
              n1 = n1,
              z1 = z1.fixed,
              bdry = c2,
              w1 = w1,
              w2 = w2,
              theta_cp = delta.tilde,
              sigma = sigma,
              ce = ce,
              cf = cf)
)

diagnostic.df$cost <- gamma * (diagnostic.df$n2x - n2)
diagnostic.df$objective <- diagnostic.df$CP - diagnostic.df$cost

best.row <- diagnostic.df[which.max(diagnostic.df$objective), ]

print(best.row)


##### Standalone Conditional Power Calculator #####################################################################

## Input values
z1.input  <- 2.294      # change as desired
n2x.input <- 170      # change as desired

## Compute Conditional Power
CP.input <- cp(
  n1       = n1,
  n2x      = n2x.input,
  z1       = z1.input,
  bdry     = c2,
  w1       = w1,
  w2       = w2,
  theta_cp = delta.tilde,
  sigma    = sigma,
  ce       = ce,
  cf       = cf
)

## Print results
cat("\n")
cat("Conditional Power Calculation\n")
cat("-----------------------------\n")
cat("z1        =", z1.input, "\n")
cat("n2x       =", n2x.input, "\n")
cat("theta_cp  =", delta.tilde, "\n")
cat("CP        =", CP.input, "\n")

#################################################################################################################

