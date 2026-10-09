# CODE FOR FIG 5.4

##This R code creates the optimal constrained promising zone design with SSR rule optimized for delta=delta.min and 
#CP evaluated at any delta=delta.true.
#The plots display the actual promising zone based on the inputs provided 
#(unlike EstDeltaDesignPlots.r where they are hard coded based on East output).
#Run this program first to get the plots
#After that, run "MinDeltaDesignSimulations-PD.r" to simulate the design and get zone wise results. 
#It will produce the plots again, on the screen

rm(list = ls())
library(ggplot2)
library(gridExtra)
library(latex2exp)
# library(ggforce)

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
      file.path("..", "..", "..", "tools", "helper.R"),
      local = .GlobalEnv
    )
  }
)

# root <- getwd()
# base.dir <- file.path(root, "Adaptive with delta-min")

# Setting up the base directory for executing this code
base.dir <- get.base.dir(
  cur.dir = getwd(),
  sec.dir = "Chapter 5/Sec 5.2/Sec 5.2.2/"
)

#input paramenters
delta.min = 7 # Smallest clinically meaningful delta; used to construct the ssr rule
delta.true = 7 # hypothesized value of delta; used to evaluate the conditional and unconditional power of the ssr rule
#Take this input from PZdesig.simulation.r
sigma = 20
n2 = 170
n1 = 85 #208
nmax = 340#170
f = nmax/n2 #multiplier to increase sample size from n2 to nmax
nmax.prom = 340 #required to create the promising zone when nmax=n2 i.e. NoSSR case
target.cp.min = 0.8 #0.585
target.cp.max = 0.9 #0.8
alpha = 0.025 #one-sided


#CP with n1 at interim and n2 at final
cp <- function(z1, delta, sigma, n1, n2)
{
  #dA <- (qnorm(0.975) * sqrt(n2) - z1 * sqrt(n1))/sqrt(n2 - n1)  use qnorm(0.975) if no early stopping
  dA <- (1.9686 * sqrt(n2) - z1 * sqrt(n1))/sqrt(n2 - n1) #2.9626 is the gsd efficacy boundary
  dMean <- delta * sqrt(n2 - n1)/(2 * sigma)
  pnorm( dA - dMean , lower.tail = FALSE)
}

#cp with n1 at interim and n2star > n2 at final
cp.star <- function(z1, delta, sigma, n1, n2, n2star)
{
  #pnorm((qnorm(0.975) * sqrt(n2) - z1 * sqrt(n1))/sqrt(n2 - n1) - delta * sqrt(n2star - n1)/(2 * sigma) , lower.tail = FALSE)
  pnorm((1.9686 * sqrt(n2) - z1 * sqrt(n1))/sqrt(n2 - n1) - delta * sqrt(n2star - n1)/(2 * sigma) , lower.tail = FALSE)
}

#input z1 and output cp.diff.z1 = cp minus target.cp
cp.diff.z1 <- function(x, delta, sigma, n1, n2, target.cp)
{
  cp(x, delta, sigma,  n1, n2) - target.cp
}

#input z1 and output cp.diff.z1.star = cp with n2star minus target.cp
cp.diff.z1.star <- function(x, delta, sigma, n1, n2, n2star, target.cp)
{
  cp.star(x, delta, sigma,  n1, n2, n2star) - target.cp
}

#input n2 and output cp.diff.n2 = cp minus target.cp 
cp.diff.n2 <- function(x, delta, sigma, n1, z1, target.cp)
{
  cp(z1, delta, sigma, n1, x)  - target.cp
}

#input n2star and output cp.diff.n2star = cp with n2* minus target.cp
cp.diff.n2star <- function(x, delta, sigma, n1, n2, z1, target.cp)
{
  cp.star(z1, delta, sigma, n1, n2, x)  - target.cp
}

#find ss such that at a given z1, cp is equal to target.cp
get.req.ss <- function(z1, n2, nmax, delta, sigma, n1, target.cp)
{
  if (cp.star(z1, delta, sigma, n1, n2, n1+1) >= target.cp)
    ss <- n2
  else
   ss <- uniroot(function(x) cp.diff.n2star(x, delta, sigma, n1, n2, z1, target.cp),c(n1+1,1000000), tol = 1.0E-12)$root
  return(min(max(ss, n2), nmax))
}

Vectorize.get.req.ss <- Vectorize(get.req.ss, SIMPLIFY =FALSE)

#find the cp that corresponds to the get.req.ss for a given z1=x
# Vidyadhar - Vectorize was needed for integrate
cp.with.ssr <- function(x, delta.min, delta.true, sigma, n1, n2, nmax, target.cp)
{
  new.ss <- unlist(Vectorize.get.req.ss(x, n2, nmax, delta.min, sigma, n1, target.cp))
  cp.star(x, delta.true, sigma,n1, n2, new.ss)
}

#Identify the promising zone (z.l, z.u and z.c)
if(n2 == nmax){
  z.l <- uniroot(function(x) cp.diff.z1.star(x,delta.min, sigma, n1, n2, nmax.prom, target.cp = target.cp.min),c(-6,6), tol = 1.0E-12)$root
  z.c <- uniroot(function(x) cp.diff.z1.star(x,delta.min, sigma, n1, n2, nmax.prom, target.cp = target.cp.max),c(-6,6), tol = 1.0E-12)$root
  z.u <- uniroot(function(x) cp.diff.z1(x,delta.min, sigma, n1, n2, target.cp = target.cp.max),c(-6,6), tol = 1.0E-12)$root 
}else{
z.l <- uniroot(function(x) cp.diff.z1.star(x,delta.min, sigma, n1, n2, nmax, target.cp = target.cp.min),c(-6,6), tol = 1.0E-12)$root
z.c <- uniroot(function(x) cp.diff.z1.star(x,delta.min, sigma, n1, n2, nmax, target.cp = target.cp.max),c(-6,6), tol = 1.0E-12)$root
z.u <- uniroot(function(x) cp.diff.z1(x,delta.min, sigma, n1, n2, target.cp = target.cp.max),c(-6,6), tol = 1.0E-12)$root
}

z.1 <- seq(-2.5, 3.5, by = 0.001)

z.1.fut <- z.1[z.1 < -1]
Z1.Fut.Data <- as.data.frame(cbind(Z.1 = z.1.fut, SS = n1))

z.1.low <- z.1[z.1 >= -1 & z.1 < z.l]
Z1.Low.Data <- as.data.frame(cbind(Z.1 = z.1.low, SS = n2))

z.1.mid <- z.1[z.1 >= z.l & z.1 <= z.u]
req.ss <- rep(n2, length(z.1.mid))
for(i in 1:length(z.1.mid))
{
  req.ss[i] <- get.req.ss(z.1.mid[i],  n2,  nmax, delta = delta.min, sigma = sigma, n1 = n1, target.cp = target.cp.max)
}
Z1.Mid.Data <- as.data.frame(cbind(Z.1 = z.1.mid, SS = req.ss))

z.1.high <- z.1[z.1 > z.u & z.1 < 2.9626]
Z1.High.Data <- as.data.frame(cbind(Z.1 = z.1.high, SS = n2))

z.1.sup <- z.1[z.1 >= 2.9626]
Z1.Sup.Data <- as.data.frame(cbind(Z.1 = z.1.sup, SS = n1))

Z1.SS.Plot.Data <- rbind(Z1.Fut.Data, Z1.Low.Data, Z1.Mid.Data, Z1.High.Data, Z1.Sup.Data)

#plot sample size versus z1
SS.Plot <- ggplot(data = Z1.SS.Plot.Data, aes(x = Z.1, y = SS)) + 
  geom_line(linewidth=1.5) + 
  
  coord_cartesian(xlim = c(-2, 3.5), ylim = c(80, 350)) +
  
  scale_x_continuous(
    breaks = seq(-2, 3.5, by = 0.5)
  ) +
  scale_y_continuous(
    breaks = seq(80, 360, by = 20)
  ) +
  
  geom_vline(linetype="dashed",  color="red", linewidth=1.5, xintercept = z.l)+
  geom_vline(linetype="dashed",  color="red", linewidth=1.5, xintercept=z.u)+
  xlab("Interim Test Statistic Z1") + 
  ylab("Sample Size") +
  annotate("text", x=1.65, y=130, label = "Promising Zone", size=10)+
  #annotate("text", x=1.5, y=115, label = "Zone", size=10)+
  #annotate("text", label= TeX("$(1.2 \\leq z_{1} \\leq 2.5) $", output = "character"), x = 1.85, y= 400,size=5, parse=TRUE)+
  annotate("text", label= TeX("$ \\leq z_{1} \\leq  $", output = "character"), x = 1.65, y= 100, size=10, parse=TRUE)+
  annotate("text", label= paste(round(z.l,4)), x = 1.65-0.55, y= 100, size=10, parse=TRUE)+
  annotate("text", label= paste(round(z.u,4)), x = 1.65+0.55, y= 100, size=10, parse=TRUE)+
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x  = element_text(size = 22),
    axis.text.y  = element_text(size = 22)
  )

 
SS.Plot 

#get the conditional power in the different zones
CP.Fut <- rep(0, sum(z.1 < -1))
CP.L <- sapply(z.1[z.1 >= -1 & z.1 < z.l], cp, delta = delta.true, sigma = sigma, n1 = n1, n2 = n2)
CP.H <- sapply(z.1[z.1 > z.u & z.1 < 2.9626], cp, delta = delta.true, sigma = sigma, n1 = n1, n2 = n2)
z.1.mid <- z.1[z.1 >= z.l & z.1 <= z.u]
CP.M <- sapply(z.1.mid, cp.with.ssr, delta.min, delta.true, sigma = sigma, n1 = n1, n2 = n2, nmax, target.cp = target.cp.max)
CP.Sup <- rep(1, sum(z.1 > 2.9626))

Z1.CP.Plot.Data <- as.data.frame(cbind(Z.1 = z.1, CP = c(CP.Fut, CP.L,  CP.M, CP.H, CP.Sup)))

## Dotted red curve: CP at delta.min with NO adaptive increase between z_L and z_U (just for reference)
z.red <- seq(z.l, z.u, by = 0.001)

Red.Curve.Data <- data.frame(
  Z.1 = z.red,
  CP  = sapply(
    z.red,
    cp,
    delta = delta.min,
    sigma = sigma,
    n1 = n1,
    n2 = n2
  )
)

#plot conditional power versus z1
CP.Plot <- ggplot(data = Z1.CP.Plot.Data, aes(x = Z.1, y = CP)) + 
  geom_line(linewidth=1.5) + 
  coord_cartesian(xlim = c(-2, 3.5), ylim = c(0, 1)) +
  scale_x_continuous(
    breaks = seq(-2, 3.5, by = 0.5)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1)
  ) +
  geom_vline(linetype="dashed",  color="red", linewidth=1.5, xintercept = z.l )+
  geom_vline(linetype="dashed",  color="red", linewidth=1.5, xintercept= z.u)+
  #ylim(0, 1) +
  #scale_y_continuous(breaks = seq(from=0, to=1,by=0.2))+
  xlab("Interim Test Statistic Z1") + 
  ylab(TeX("Conditional Power Evaluated at $\\delta = 7$")) +
  annotate("text", x=1.65, y=0.25, label = "Promising Zone", size=10)+
  #annotate("text", x=1.65, y=0.2, label = "Zone", size=8)+
  #annotate("text", label= TeX("$(1.2 \\leq z_{1} \\leq 2.5) $", output = "character"), x = 1.85, y= 0.15, size=5, parse=TRUE)+
  annotate("text", label= TeX("$ \\leq z_{1} \\leq  $", output = "character"), x = 1.65, y= 0.15, size=10, parse=TRUE)+
  annotate("text", label= paste(round(z.l,4)), x = 1.65-0.55, y= 0.15, size=10, parse=TRUE)+
  annotate("text", label= paste(round(z.u,4)), x = 1.65+0.55, y= 0.15, size=10, parse=TRUE)+
  
  geom_line(
    data = Red.Curve.Data,
    aes(x = Z.1, y = CP),
    inherit.aes = FALSE,
    color = "red",
    linetype = "dotted",
    linewidth = 1.5
  ) +
  
  
  
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x  = element_text(size = 22),
    axis.text.y  = element_text(size = 22)
  )

CP.Plot

#arrange the plots in one column

combined_plot <- grid.arrange(SS.Plot, CP.Plot, ncol = 1)
combined_plot

# fig.dir <- "D:/Cyrus-Office/BookProject/All-Chapters/screens"
fig.dir <- base.dir

ggsave(
  filename = file.path(fig.dir, "PD_SampleSizeRuleMinDelta.png"),
  plot     = combined_plot,
  width    = 15,
  height   = 15,
  units    = "in",
  dpi      = 300
)

# ####Also save in base.dir
# #### Also save in base.dir

# ggsave(
#   filename = file.path(base.dir, "PD_SampleSizeRuleMinDelta.png"),
#   plot     = combined_plot,
#   width    = 15,
#   height   = 15,
#   units    = "in",
#   dpi      = 300,
#   bg       = "white"
# )
