
# Simulates optimal promising zone design, whose SSR rule is optimized for delta=delta.min, at any value of delta=delta.true.
#The outputs are the power and sample size both zone wise and overall

rm(list=ls())
#source("D:\\Cyrus-Office\\BookProject\\All-Chapters\\Chapter4\\Software for Chapter 4\\Adaptive with delta-min\\PZdesign.plots.r")

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
  sec.dir = "Chapter 5/Sec 5.2.2/"
)

# root <- getwd()
# base.dir <- file.path(root, "Adaptive with delta-min")
source(file.path(base.dir, "MinDeltaDesignPlots-PD.r"))

#source("Adaptive with delta-min\\MinDeltaDesignPlots-PD.r")

#delta.min and sigma are taken from MinDeltaDesignPlots-PD.R and is used in the get.req.ss function. 
#Simulations are then performed with delta.true to get the zone wise and unconditional power and sample size

delta.true = 10

start_time <- Sys.time()
num.sims <- 10000
mean.1 <- delta.true * sqrt(n1)/(2 * sigma)
sim.data <- as.data.frame(matrix(0, nrow = num.sims, ncol = 7))
colnames(sim.data) <- c("Sim.Ind", "Z.1", "CP", "Zone", "Z.2", "Success", "Sample.Size")


z1 <- rnorm(num.sims, mean = mean.1, sd = 1)
sim.data[, "Sim.Ind"] <- 1:num.sims
sim.data[, "Z.1"] <- z1
sim.data[, "CP"] <-  sapply(z1, cp, delta = delta.true, sigma = sigma, n1 = n1, n2 = n2)  

# split the simulation data, based on interim results
sim.data.prom <- sim.data[sim.data$Z.1 >= z.l &  sim.data$Z.1 <= z.u, ]

# simulation results in promising zone
sim.data.prom$Zone = 2
n2.prom <- sapply(sim.data.prom$Z.1, get.req.ss,  n2,  nmax, delta = delta.min, sigma = sigma, n1 = n1, target.cp = target.cp.max)
n2.prom.incr <- n2.prom - n1
mean.2.prom <- delta.true * sqrt(n2.prom.incr)/(2 * sigma)
z2.prom.incr <- rnorm(nrow(sim.data.prom), mean = mean.2.prom, sd = 1)

#sim.data.prom[, "Z.2"] <- sqrt(n1/n2.prom) * sim.data.prom$Z.1 + sqrt(n2.prom.incr/n2.prom) * z2.prom.incr: FUNDAMENTAL ERROR!
sim.data.prom[, "Z.2"] <- sqrt(n1/n2) * sim.data.prom$Z.1 + sqrt((n2-n1)/n2) * z2.prom.incr

#sim.data.prom[, "Success"] <- ifelse(sim.data.prom[, "Z.2"] > qnorm(0.975), 1, 0) (in case of no early efficacy stopping)
sim.data.prom[, "Success"] <- ifelse(sim.data.prom[, "Z.2"] > 1.9686, 1, 0) #(LDOF critical value)
sim.data.prom[, "Sample.Size"] <- n2.prom

# simulate the rest of the data in the non-promising zones
sim.data.nonprom <- sim.data[sim.data$Z.1 <= z.l | sim.data$Z.1 > z.u, ]
sim.data.nonprom$Zone <- ifelse(sim.data.nonprom$Z.1 < z.l, 1, 3)
n2.nonprom <- n2
n2.nonprom.incr <- n2.nonprom - n1
mean.2.nonprom <- delta.true * sqrt(n2.nonprom.incr)/(2 * sigma)
z2.nonprom.incr <- rnorm(nrow(sim.data.nonprom), mean = mean.2.nonprom, sd = 1)
sim.data.nonprom[, "Z.2"] <- sqrt(n1/n2.nonprom) * sim.data.nonprom$Z.1 + sqrt(n2.nonprom.incr/n2.nonprom) * z2.nonprom.incr
sim.data.nonprom[, "Success"] <- ifelse(sim.data.nonprom[, "Z.2"] > qnorm(0.975), 1, 0)
sim.data.nonprom[, "Sample.Size"] <- n2

# subdivide sim.data.nonprom into four data frames

# 1. Futility zone: Z.1 < -1
sim.data.nonprom.fut <- sim.data.nonprom[sim.data.nonprom$Z.1 < -1, ]
if (nrow(sim.data.nonprom.fut) > 0) {
  sim.data.nonprom.fut$CP <- 0
  sim.data.nonprom.fut$Zone <- 0
  sim.data.nonprom.fut$Success <- 0
  sim.data.nonprom.fut$Sample.Size <- n1   # 85
}

# 2. Unfavorable zone: -1 <= Z.1 < z.l
sim.data.nonprom.unf <- sim.data.nonprom[sim.data.nonprom$Z.1 >= -1 &
                                           sim.data.nonprom$Z.1 < z.l, ]

# 3. Favorable zone: z.u < Z.1 < 2.9626
sim.data.nonprom.fav <- sim.data.nonprom[sim.data.nonprom$Z.1 > z.u &
                                           sim.data.nonprom$Z.1 < 2.9626, ]

# 4. Efficacy zone: Z.1 >= 2.9626
sim.data.nonprom.effi <- sim.data.nonprom[sim.data.nonprom$Z.1 >= 2.9626, ]
if (nrow(sim.data.nonprom.effi) > 0) {
  sim.data.nonprom.effi$Zone <- 4
  sim.data.nonprom.effi$Success <- 1
  sim.data.nonprom.effi$Sample.Size <- n1   # 85
}

#sim.data <- rbind(sim.data.nonprom, sim.data.prom) 

# combine everything
sim.data.big <- rbind(sim.data.nonprom.fut,
                      sim.data.nonprom.unf,
                      sim.data.prom,
                      sim.data.nonprom.fav,
                      sim.data.nonprom.effi)


# Present the Results
fut  <- sim.data.nonprom.fut
unfav <- sim.data.nonprom.unf
prom <- sim.data.prom
fav  <- sim.data.nonprom.fav
effi <- sim.data.nonprom.effi

Zone.Ind <- c("Futility", "Unfavorable", "Promising", "Favorable", "Efficacy", "Unconditional")

Prob.Zone <- c(nrow(fut), nrow(unfav), nrow(prom), nrow(fav), nrow(effi), nrow(sim.data.big)) / num.sims

Power.Zone <- c(sum(fut$Success),
                sum(unfav$Success),
                sum(prom$Success),
                sum(fav$Success),
                sum(effi$Success),
                sum(sim.data.big$Success)) / (num.sims * Prob.Zone)

Samp.Size.Zone <- c(mean(fut$Sample.Size),
                    mean(unfav$Sample.Size),
                    mean(prom$Sample.Size),
                    mean(fav$Sample.Size),
                    mean(effi$Sample.Size),
                    mean(sim.data.big$Sample.Size))

results <- data.frame(Zone.Ind, Prob.Zone, Power.Zone, Samp.Size.Zone,
                      delta.true, sigma, n1, n2, nmax)

results

#Combine futility and unfavorable zones. And combine favorable and efficacy zones.
# ---- Additional combined zone summaries ----


fut_unf <- rbind(fut, unfav)

# Combine Favorable + Efficacy
fav_effi <- rbind(fav, effi)

Zone.Ind.comb <- c("Fut+Unf", "Promising", "Fav+Effi", "Unconditional")

Prob.Zone.comb <- c(nrow(fut_unf),
                    nrow(prom),
                    nrow(fav_effi),
                    nrow(sim.data.big)) / num.sims

Power.Zone.comb <- c(sum(fut_unf$Success),
                     sum(prom$Success),
                     sum(fav_effi$Success),
                     sum(sim.data.big$Success)) / (num.sims * Prob.Zone.comb)

Samp.Size.Zone.comb <- c(mean(fut_unf$Sample.Size),
                         mean(prom$Sample.Size),
                         mean(fav_effi$Sample.Size),
                         mean(sim.data.big$Sample.Size))

results.comb <- data.frame(Zone.Ind = Zone.Ind.comb,
                           Prob.Zone = Prob.Zone.comb,
                           Power.Zone = Power.Zone.comb,
                           Samp.Size.Zone = Samp.Size.Zone.comb,
                           delta.true, sigma, n1, n2, nmax)

results.comb


end_time <- Sys.time()

end_time - start_time

