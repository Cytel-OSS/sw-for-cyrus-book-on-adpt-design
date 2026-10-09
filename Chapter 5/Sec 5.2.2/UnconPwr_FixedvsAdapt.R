######### To compute and compare the unconditional power ####################
## For fixed sample design against adaptive design #########
############################################################################


rm(list=ls())
source("D:\\Cyrus-Office\\BookProject\\All-Chapters\\Chapter4\\Software for this chapter\\Adaptive with delta-min\\PZdesign.plots.R")


## Integrand for unconditional power computation for non-promising zones: evaluated at delta = delta.true
Unc.integrand1 <- Vectorize(function(z1, delta)
{
  f.del.z1 <-  dnorm(z1, mean = delta*sqrt(n1)/(2*sigma), sd = 1)
  term <- cp(z1, delta, sigma, n1, n2)*  f.del.z1
})

## Integrand for unconditional power computation for promising zones: evaluated at delta = delta.true
## Notice that target.cp=target.cp.max is passed as the argument for cp.with.ssr. This is because we want to maximize the conditional power in the promising zone subject to never exceeding target.cp.max. Also delta.min is passed as an argument of cp.with.ssr and will be used when cp.with.ssr calls up get.req.ss
Unc.integrand2 <- Vectorize(function(z1, delta)
{
  f.del.z1 <- dnorm(z1, mean = delta*sqrt(n1)/(2*sigma), sd = 1)
  term <- cp.with.ssr(z1, delta.min, delta, sigma, n1, n2, nmax, target.cp = target.cp.max)* f.del.z1 
})


## Unconditional Power computation for adaptive design: evaluated at delta = delta.true
UnconPwr.adpt <-Vectorize( function(delta)
{
  #delta = del.inp
  pwr.adpt <- integrate( Unc.integrand1,delta, lower = -Inf, upper = z.l )$value + 
              integrate( Unc.integrand1,delta, lower = z.u, upper = Inf)$value + 
              integrate( Unc.integrand2,delta, lower = z.l, upper = z.u)$value
  return(pwr.adpt)
}
)

## Integrand for expected sample size computation outside the promising zone: evaluated at delta = delta.true 
N.del.integrand1 <- function(z1,delta)
{
  f.del.z1 <- dnorm(z1, mean = delta*sqrt(n1)/(2*sigma), sd = 1)
  Ndel.ad <- n2* f.del.z1
}

## Integrand for expected sample size computation inside the promising zone: evaluated at delta = delta.true
N.del.integrand2 <- function(z1,delta)
{
  f.del.z1 <- dnorm(z1, mean = delta*sqrt(n1)/(2*sigma), sd = 1)
  Ndel.ad <- unlist(Vectorize.get.req.ss(z1, n2, nmax, delta.min, sigma, n1, target.cp.max))*f.del.z1 # Vidyadhar - Vectorize was needed for integrate
}

# Computation of expected sample size of adaptive design for a given delta
SS.adapt <- Vectorize(function(delta)
{
  N.delta.adapt <- integrate(N.del.integrand1, delta, lower = -Inf, upper = z.l)$value +
                   integrate(N.del.integrand1, delta, lower = z.u, upper = Inf)$value + 
                   integrate(N.del.integrand2, delta, lower = z.l, upper =z.u)$value
  return(N.delta.adapt)
})

## Unconditional Power computation for fixed sample design
UnconPwr.fixed <- function( delta)
{
  fixed.pwr <- pnorm( (sqrt(SS.adapt(delta) )* delta/ (2* sigma) )- qnorm(1- alpha) )
  return(fixed.pwr)
}

#Create dataset for plotting
delta.range <- seq(from=1.6, to = 2, by = 0.05)
delta.range_final <- rep(seq(from=1.6, to = 2, by = 0.05),2)
pwr.range <- c(UnconPwr.adpt(delta.range),UnconPwr.fixed(delta.range))
type <- rep(c("Adaptive", "Fixed"), each = length(seq(from=1.6, to = 2, by = 0.05)))
pwr.plot.data2 <- data.frame(delta_range = delta.range_final, Unc_Power = type, Unconditional_Power = pwr.range)
Unc.Pwr.Plot <- ggplot(pwr.plot.data2, aes(x = delta_range, y = Unconditional_Power, color = Unc_Power)) +
  geom_line() + 
  #ylim(c(0, 1))+
  scale_y_continuous(breaks = seq(from=0, to=1,by=0.1), limits = c(0, 1))+
  #ggtitle("Title") +
  #xlab("XLAB") +
  #ylab("YLAB") + 
  #theme_minimal()
  theme(text=element_text(size = 15),
        legend.position = c(0.2, 0.85))

Unc.Pwr.Plot



## Code for Creating Conditional Power Plots over a range of delta values
# first get a range of z1 values for each of the zones
z.1 <- seq(0, 3, by = 0.001)
z.1.low <- z.1[z.1 < z.l]
Z1.Low.Data <- as.data.frame(cbind(Z.1 = z.1.low, SS = n2))

z.1.mid <- z.1[z.1 >= z.l & z.1 <= z.u]
req.ss <- rep(n2, length(z.1.mid))
for(i in 1:length(z.1.mid))
{
  req.ss[i] <- get.req.ss(z.1.mid[i],  n2,  nmax, delta = delta.min, sigma = sigma, n1 = n1, target.cp = target.cp.max)
}
Z1.Mid.Data <- as.data.frame(cbind(Z.1 = z.1.mid, SS = req.ss))

z.1.high <- z.1[z.1 > z.u]
Z1.High.Data <- as.data.frame(cbind(Z.1 = z.1.high, SS = n2))


# now get the conditional power in the different zones at delta = delta.true
delta.true = 1.6
CP.L.d1.6 <- sapply(z.1[z.1 < z.l], cp, delta = delta.true, sigma = sigma, n1 = n1, n2 = n2)
CP.H.d1.6 <- sapply(z.1[z.1 > z.u], cp, delta = delta.true, sigma = sigma, n1 = n1, n2 = n2)
z.1.mid <- z.1[z.1 >= z.l & z.1 <= z.u]
CP.M.d1.6 <- sapply(z.1.mid, cp.with.ssr, delta.min, delta.true, sigma = sigma, n1 = n1, n2 = n2, nmax, target.cp = target.cp.max)

delta.true = 1.8
CP.L.d1.8 <- sapply(z.1[z.1 < z.l], cp, delta = delta.true, sigma = sigma, n1 = n1, n2 = n2)
CP.H.d1.8 <- sapply(z.1[z.1 > z.u], cp, delta = delta.true, sigma = sigma, n1 = n1, n2 = n2)
CP.M.d1.8 <- sapply(z.1.mid, cp.with.ssr, delta.min, delta.true, sigma = sigma, n1 = n1, n2 = n2, nmax, target.cp = target.cp.max)

delta.true = 2
CP.L.d2 <- sapply(z.1[z.1 < z.l], cp, delta = delta.true, sigma = sigma, n1 = n1, n2 = n2)
CP.H.d2 <- sapply(z.1[z.1 > z.u], cp, delta = delta.true, sigma = sigma, n1 = n1, n2 = n2)
CP.M.d2 <- sapply(z.1.mid, cp.with.ssr, delta.min, delta.true, sigma = sigma, n1 = n1, n2 = n2, nmax, target.cp = target.cp.max)

delta.vec <- factor(rep(c(1.6, 1.8, 2), each = 3001))

CP <- c(CP.L.d1.6, CP.M.d1.6, CP.H.d1.6,
        CP.L.d1.8, CP.M.d1.8, CP.H.d1.8,
        CP.L.d2,   CP.M.d2,   CP.H.d2)

#plot conditional power versus z1
Z1.CP.Plot.Data <- data.frame(Z.1 = z.1, CP = CP, delta = delta.vec)
CP.Plot <- ggplot(data = Z1.CP.Plot.Data, aes(x = Z.1, y = CP)) + 
  geom_line(aes(group = delta, col = delta)) + 
  geom_vline(linetype="dashed",  color="red", xintercept = z.l )+
  geom_vline(linetype="dashed",  color="red", xintercept= z.u)+
  #ylim(0, 1) +
  scale_x_continuous(breaks = seq(from=0, to=3,by=0.5))+
  scale_y_continuous(breaks = seq(from=0, to=1,by=0.1), limits = c(0, 1))+
  xlab("Z at Interim") + 
  ylab("Conditional Power") +
  annotate("text", x=1.85, y=0.25, label = "Promising", size=5)+
  annotate("text", x=1.85, y=0.2, label = "Zone", size=5)+
  annotate("text", label= TeX("$(1.2 \\leq z_{1} \\leq 2.5) $", output = "character"), x = 1.85, y= 0.15, size=5, parse=TRUE)+
  theme(
    text=element_text(size = 15),
    legend.position = c(0.15, 0.85))
CP.Plot

grid.arrange(Unc.Pwr.Plot, CP.Plot, ncol = 2)
