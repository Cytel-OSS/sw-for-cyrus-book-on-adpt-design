## Definitions

# sigma: design standard deviation
# n1: sample size at stage 1
# n2: initially planned final sample
# n2star: adapted/new final sample size
# n2x: target sample size to obtain certain level of conditional power
# z1: interim test statistic
# z.alpha: final look critical boundary
# ce: early efficacy stopping boundary
# cf: early futility stopping boundary
# target.cp: target conditional power 


##---------------------------------------------------------------------------------------------------
## Function for the CPZ Method

# CP with n1 at interim and n2 at final
cp <- function(z1, delta, sigma, n1, n2, z.alpha, ce, cf)
{
  if( cf < z1 & z1 < ce) {
    
    dA <- (z.alpha * sqrt(n2) - z1 * sqrt(n1))/sqrt(n2 - n1)
    dMean <- delta * sqrt(n2 - n1)/(2 * sigma)
    pnorm( dA - dMean , lower.tail = FALSE)
    
    } else {
      
    ifelse(z1 > ce, 1.0, 0.0 )
    
  }
  
}

# CP with n1 at interim and n2star > n2 at final
cp.star <- function(z1, delta, sigma, n1, n2, n2star, z.alpha, ce, cf)
{
  if(z1> cf & z1<ce){
    
    pnorm((z.alpha*sqrt(n2) - z1*sqrt(n1))/sqrt(n2 - n1) - delta*sqrt(n2star - n1)/(2*sigma),
          lower.tail = FALSE)
    
  } else {
    
    ifelse(z1 > ce, 1.0, 0.0)
    
  }
  
}

# cp.diff.z1 is used to obtain the interim test statistic, z1 required 
# to attain target sample size and conditional power
cp.diff.z1<- function(x, delta, sigma, n1, n2, n2x, target.cp, ce, cf)
{
   cp.star(x, delta, sigma,  n1, n2, n2star=n2x, z.alpha, ce, cf) - target.cp
}

# obtain test statistic value for given target sample size and target conditional power
get.z1.cp <- function(n1, n2, n2x, z.alpha, delta, sigma, target.cp, ce, cf) {
  
  uniroot(function(x) cp.diff.z1(x, delta, sigma, n1, n2, n2x, target.cp, ce, cf),
          c(-6,6), tol = 1.0E-12)$root

}


# cp.diff.n2star used to obtained required sample size to attain the target CP
cp.diff.n2star <- function(x, delta, sigma, n1, n2, z1, z.alpha, target.cp, ce, cf)
{
  cp.star(z1, delta, sigma, n1, n2, x, z.alpha, ce, cf) - target.cp
}


## obtain required sample size to attain target CP
get.n2star.cp<- function(z1, n2, n2x, delta, sigma, n1, z.alpha, target.cp, ce, cf)
{
  
  ss <- uniroot(function(x) cp.diff.n2star(x, delta, sigma, n1, n2, z1, z.alpha,
                                           target.cp, ce, cf),c(n1+1,10000), tol = 1.0E-12)$root
  return(min(max(ss, n2), n2x))
}


## evaluate sample size function n2(z1) at each values of z1
ss.cpz.function<-function(z1, zl, zc, zu,  ce, cf, z.alpha, target.cp,
                          n1, sigma, delta, n2, n2x){
  
  n2_z1 <-ifelse(z1<= cf | z1>ce, n1,
          ifelse((z1 > cf & z1 < zl )| (z1 < ce & z1 > zu ), n2,
          ifelse((zl < z1 & z1 < zc), n2x, 
                 get.n2star.cp(z1, n2, n2x, delta, sigma, n1, 
                                 z.alpha, target.cp, ce, cf)))) 
  
  
  return(n2_z1)
  
}


## obtain boundaries, sample size rule and conditional power for the CPZ approach
get.ssr.bdry.cpz <- function(z1, n1, n2, k, ce, cf, z.alpha, delta, sigma, cpmin, cpmax){
  
  n2max<- n2*k
  
  z.l <- get.z1.cp(n1, n2, n2x=n2max, z.alpha, delta, sigma, target.cp=cpmin, ce, cf) 
  z.c <- get.z1.cp(n1, n2, n2x=n2max, z.alpha, delta, sigma, target.cp=cpmax, ce, cf)
  z.u <- get.z1.cp(n1, n2, n2x=n2, z.alpha, delta, sigma, target.cp=cpmax, ce, cf)
  
  
  ss <-  sapply(z1, ss.cpz.function, zl=z.l, zc=z.c, zu=z.u,
                ce=ce, cf=cf, z.alpha=z.alpha, target.cp=cpmax,
                n1=n1, sigma=sigma, delta=delta, n2=n2, n2x=n2max)
  
  cps<- mapply(FUN = cp.star, z1,  n2 = n2, n2star=ss, n1=n1, 
               delta=delta, sigma=sigma, z.alpha=z.alpha, ce=ce, cf=cf)
  
  
  zone <-  sapply(z1, get.zone, zl=z.l, zc=z.c, zu=z.u, ce=ce, cf=cf)
  
  promising_bdry <- c(z.l=z.l, z.c=z.c, z.u=z.u)
  early_stop_bdry <- c(cf=cf, ce=ce)
  
  SSR.df <- data.frame(z1, n2_z1=ss, 
                       ConditionalPower=cps,Zone=zone)
  
  return( ssr=list(early_stop_bdry = early_stop_bdry, 
                   promising_bdry = promising_bdry, SSR.df = SSR.df ))
}


# get.zone <- function(z1, zl, zc, zu, ce, cf){
#   
#   Zone <- ifelse(z1< cf , "Futility", 
#                  ifelse((z1 > cf & z1 < zl ), "Unfavorable",
#                         ifelse((zl < z1 & z1 < zu),  "Promising",
#                                ifelse((z1 > zu & z1 < ce ), "Favorable",  "Efficacy"))))
#   
#   
# }

get.zone <- function(z1, zl, zc, zu, ce, cf){
  
  Zone <- ifelse(z1 < zl , "Futility+Unfavorable",
                 ifelse((zl < z1 & z1 < zu),  "Promising", "Favorable + Efficacy"))
  
  
}


get.power <- function(delta.true, z1, zl, zc, zu, ce, cf, z.alpha, n1, n2, n2star, sigma){
  
  
  
  region <- ifelse((z1<cf | z1>ce), "Early Stop", 
                   ifelse((z1 > cf & z1 < zl | z1 < ce & z1 > zu),"Non-Promising",
                          "Promising") )
  
  n2.incr <-n2star-n1
  df <- data.frame(z1,region,n2star, n2.incr)%>%
    mutate(mean.2.incr= delta.true * sqrt(n2.incr)/(2 * sigma))
  
  for(i in 1:length(z1)){
    
    if(region[i]=="Early Stop"){
      
      df$z2.incr[i]<-0
      
    }else{
      
      df$z2.incr[i]<-  rnorm(1, df$mean.2.incr[i], sd=1 )
      
    }
    
  }             
  
  power.df <- df%>% mutate( z2 = ifelse(region=="Early Stop",0,
                                        sqrt(n1/n2) * z1 + sqrt((n2-n1)/n2) * z2.incr),
                            Success= ifelse(z1<cf, 0,
                                            ifelse(z1>ce, 1,
                                                   ifelse(z2 > z.alpha,1,0))))
  
  
  
  return(data.frame(power.df ))
  
}




## obtain boundaries, sample size rule,conditional power and zonewise summaries 
## for the CPZ approach
SSR.Power.Zonewise.CPZ <- function(cpmin, cpmax, n1, n2, k, ce, cf, z.alpha, 
                          delta, delta.true, sigma, n.sim=100000, seed=NULL){
  
  mean.z1 <- delta.true * sqrt(n1)/(2 * sigma)
  
  z1 <- rnorm(n.sim, mean.z1, sd=1)
  
  ssr=get.ssr.bdry.cpz(z1, n1=n1, n2=n2, k=k, ce=ce, cf=cf,
                  z.alpha=z.alpha, delta=delta, sigma=sigma,
                  cpmin=cpmin, cpmax=cpmax)
  

  n2star<-ssr$SSR.df$n2_z1
  
  z.l <- ssr$promising_bdry["z.l"]
  z.c <- ssr$promising_bdry["z.c"]
  z.u <- ssr$promising_bdry["z.u"]
  
  pow <- get.power(delta.true=delta.true, z1,  zl=z.l, zc=z.c, zu=z.u, ce=ce, cf=cf, z.alpha=z.alpha,
                   n1=n1, n2=n2, n2star=n2star, sigma=sigma)
  
  
  ssr.power.df <- cbind(ssr$SSR.df, Success=pow$Success)
  
  zonewise <- ssr.power.df %>% group_by(Zone)%>%
    summarize(Prob_Zone = n()/n.sim*100,
              Average_SS = mean(n2_z1),
              Power = mean(Success)*100,
              Count_simulations=n())
  
  zonewise.df <- data.frame(True_delta = delta.true, zonewise )
  
  all.trials<- data.frame(True_delta = delta.true, Zone="All Trials",
                          Prob_Zone = length(z1)/n.sim*100,
                          Average_SS = mean(ssr.power.df$n2_z1),
                          Power = mean(ssr.power.df$Success)*100,
                          Count_simulations=length(z1))
  
  zonewise_summaries<- rbind(zonewise.df, all.trials )
  
  return(list(ssr.power.df=ssr.power.df, zonewise_summaries=zonewise_summaries,
              boundaries=data.frame(cf=cf, z.l=z.l, z.c=z.c, z.u=z.u, ce=ce)))
  
  # return(list(ssr.power.df=ssr.power.df, zonewise_summaries=zonewise,
  #             boundaries=data.frame(cf=cf, z.l=z.l, z.c=z.c, z.u=z.u, ce=ce)))
  
  # return(list(ssr.power.df=ssr.power.df, all_trials_summaries=all.trials,
  #             boundaries=data.frame(cf=cf, z.l=z.l, z.c=z.c, z.u=z.u, ce=ce)))
}



## pow.ss.diff.deltas.cpz is function used to obtain operating characteristics-power and sample size at
## different values of true delta/ response generation delta for the CPZ approach
pow.ss.diff.deltas.cpz <-function(delta.true, delta, zl, zc, zu, n1, n2, n2max, ce, cf, z.alpha, 
                               cpmin, cpmax, n.sim=10000, seed=NULL){
  
  mean.z1 <- delta.true * sqrt(n1)/(2 * sigma)
  
  z1 <- rnorm(n.sim, mean.z1, sd=1)
  
  ss <-  sapply(z1, ss.cpz.function, zl=zl, zc=zc, zu=zu, 
                ce=ce, cf=cf, z.alpha=z.alpha, target.cp=cpmax,
                n1=n1, n2=n2, n2x=n2max, delta=delta, sigma=sigma)
  

  pow <- get.power(delta.true=delta.true, z1,  zl=zl, zc=zc, zu=zu,  ce=ce, cf=cf, z.alpha=z.alpha,
                   n1=n1, n2=n2, n2star=ss, sigma=sigma)
  
  zone <- get.zone(z1,  zl=zl, zc=zc, zu=zu, ce=ce, cf=cf)
  
  pow.ss.df <- cbind(True_delta=delta.true, pow, zone)
  
  return(pow.ss.df)
}


## operating characteristics- power and sample size at 
## different values of true delta/ response generation delta

OC.Power.SS.deltas.CPZ<- function(delta.vec, delta, n1, n2, k, ce, cf, z.alpha, 
                                   cpmin, cpmax, n.sim=10000, seed=NULL){
  n2max<-n2*k
  
  z.l <- get.z1.cp( n1=n1, n2=n2, n2x=n2max, z.alpha=z.alpha, delta=delta,
               sigma=sigma, target.cp=cpmin, ce=ce, cf=cf)
 
  z.u <- get.z1.cp( n1=n1, n2=n2, n2x=n2, z.alpha=z.alpha, delta=delta,
                sigma=sigma, target.cp=cpmax, ce=ce, cf=cf)
  
  z.c <- get.z1.cp( n1=n1, n2=n2, n2x=n2max, z.alpha=z.alpha, delta=delta,
         sigma=sigma, target.cp=cpmax, ce=ce, cf=cf)
  
  
  out <-bind_rows(lapply(delta.vec, pow.ss.diff.deltas.cpz, delta=delta,
                          zl=z.l, zc=z.c, zu=z.u,
                          n1=n1, n2=n2, n2max=n2max, ce=ce, cf=cf, z.alpha=z.alpha, 
                          cpmin=cpmin, cpmax=cpmax, n.sim=n.sim, seed=seed))
  

  delta_wise_OC <- out %>% group_by(True_delta)%>%
                   summarize(Average_SS = mean(n2star),
                   Power = mean(Success)*100,
                   Total_simulations=n())
  
  
  
  return(list(delta_wise_OC=delta_wise_OC,
              boundaries=data.frame(cf=cf, z.l=z.l, z.c=z.c, z.u=z.u, ce=ce)))
}


##---------------------------------------------------------------------------------------------------
## Function for the PPZ Method


# mean of the prior density 
prior.mean.m0 <- function(delta.l,delta.u, tails){
  
  s0 <- (delta.u-delta.l)/(qnorm(1-tails)-qnorm(tails))
  m0 <- delta.l-s0*qnorm(tails)
  return(m0)
  
}

# standard deviation of the prior density 
prior.sd.s0 <- function(delta.l,delta.u, tails){
  
  s0 <- (delta.u-delta.l)/(qnorm(1-tails)-qnorm(tails))
  return(s0)
  
}

## prior density function truncated between delta.lower and delta.upper
## evaluated at a given value of delta
prior.density.delta<-function(m0, s0, delta.lower, delta.upper, delta){
  
  c <-  1/(sqrt(2*pi)*((pnorm(delta.upper, mean=m0, sd=s0)-pnorm(delta.lower, mean=m0, sd=s0))))
  c/s0*exp(-0.5*((delta-m0)/(s0))^2)
  
}

# prior density 
prior_density<-function(m0, s0, sigma, delta.lower= -Inf, delta.upper= Inf, delta){

  pr_density<- sapply(delta, prior.density.delta, m0=m0, s0=s0,
                      delta.lower=delta.lower, delta.upper=delta.upper)
  
  prior_density.df <- data.frame(delta=delta, prior_density=pr_density)
  
  return(list(prior_density_parameters=list(s0=s0, m0=m0, n0=n0),
              prior_density = prior_density.df))
}

# posterior mean 
post.mean.m1 <- function(m0, s0, n1, sigma, z1){
  
  sigma1 <- 2*sigma/sqrt(n1)
  delta.hat <- 2*sigma*z1/sqrt(n1)
  
  m1 <- (delta.hat*s0^2 + m0*sigma1^2)/(s0^2 + sigma1^2)
  
  return(m1)
}

# posterior standard deviation
post.sd.s1 <- function(s0, n1, sigma){
  
  sigma1 <- 2*sigma/sqrt(n1)
  
  s1<- s0*sigma1/sqrt((s0^2 + sigma1^2))
  
  return(s1)
}

## posterior density function truncated between delta.lower and delta.upper
## evaluated at a given value of delta
post.density.delta <- function(m1, s1, delta.lower= -Inf, delta.upper= Inf, delta){
  
  K <-  1/(sqrt(2*pi)*s1*((pnorm(delta.upper, mean=m1, sd=s1)-pnorm(delta.lower, mean=m1, sd=s1))))
  post <- K*exp(-0.5*((delta-m1)/(s1))^2) 
  
  return (post)
}


# prior density 
post.density <- function(m0, s0, sigma, n1,  z1, delta.lower= -Inf, delta.upper= Inf, delta){
  
  pr_density<- sapply(delta, prior_density_delta, m0=m0, s0=s0,
                      delta.lower=delta.lower, delta.upper=delta.upper)
  
  m1 <- post.mean.m1(m0=m0, s0=s0, n1=n1, sigma=sigma, z1=z1)
  s1 <- post.sd.s1(s0=s0, n1=n1, sigma=sigma)
  
  sigma1 <- 2*sigma/sqrt(n1)
  delta.hat <- 2*sigma*z1/sqrt(n1)

  post <- sapply (delta, post.density.delta, m1=m1, s1=s1,
                  delta.lower=delta.lower, delta.upper=delta.upper)
  
  
  post.density.df <- data.frame( delta = delta, prior = pr_density,
                                 posterior=post)
 
  return(list(prior_params =list(delta.lower=delta.lower, delta.upper=delta.upper, m0=m0, s0=s0),
              posterior_params = list(m1=m1, s1=s1),
              design_params = list( delta=delta, z1=z1, delta=delta.hat, 
                                    sigma1=sigma1), 
              post_density_df = post.density.df))
}


## integrand for finding predictive power (no truncation)
PP_function<-function(m1, s1, n1, sigma, delta,  z1, n2, n2x, z.alpha){
  
  n2star <- n2x
  
  
  K <-  1/(sqrt(2*pi)*s1)
  
  term1 <- K*exp(-0.5*((delta-m1)/(s1))^2)
  term2 <- pnorm((z.alpha*sqrt(n2) - z1*sqrt(n1))/sqrt(n2 - n1) - delta*sqrt(n2star - n1)/(2*sigma),
                 lower.tail = FALSE)
  
  term1*term2
  
}

## predictive power
pp<-function(m0, s0, n1, sigma, z1, n2, n2x, z.alpha, ce,cf){
  
  sigma1 <- 2*sigma/sqrt(n1)
  delta.hat <- 2*sigma*z1/sqrt(n1)
  
  m1 <- (delta.hat*s0^2 + m0*sigma1^2)/(s0^2 + sigma1^2)
  s1<- s0*sigma1/sqrt((s0^2 + sigma1^2))
  
  if(z1>cf & z1<ce){
    
    integrate(PP_function, lower = m1-10*s1, upper= m1+10*s1,
              m1=m1, s1=s1,  n1=n1, sigma=sigma,
              z1=z1, n2=n2, n2x=n2x, z.alpha=z.alpha)$value
  } else{
    
    ifelse(z1 > ce, 1.0, 0.0 )
    
  }
  
}


# pp.diff.z1 is used to obtain the interim test statistic, z1 required 
# to attain target sample size and predictive power
pp.diff.z1 <- function(target.pp, m0, s0, n1, sigma, 
                       z1, n2, n2x, z.alpha, ce,cf){
  
  target.pp - pp(m0, s0, n1, sigma, z1, n2, n2x, z.alpha, ce,cf)
  
}


# obtain the interim test statistic, z1 required 
# to attain target sample size and predictive power
get.z1.pp<-function(target.pp, m0, s0, n1, sigma,
                   n2, n2x, z.alpha, ce,cf){
  
  uniroot(function(x) pp.diff.z1(target.pp, m0, s0, n1, sigma, 
                                 x, n2, n2x=n2x, z.alpha, ce,cf),
          c(-6,6), tol = 1.0E-12)$root 
  
}

## obtain required sample size to attain target PP
get.n2star.pp <- function(z1,target.pp, m0, s0, n1, sigma, 
                          n2, n2x, z.alpha, ce, cf)
{
  ss <-  uniroot(function(x) pp.diff.z1(target.pp, m0, s0, n1, sigma, 
                                        z1, n2, x, z.alpha, ce,cf),
                 c(n1+1, 10000), tol = 1.0E-12)$root 
  
  return(min(max(ss, n2), n2x))
}


## evaluate sample size function n2(z1) at each values of z1
ss.ppz.function<-function(z1, zl, zc, zu,  ce, cf, z.alpha, target.pp, m0, s0, n1,
                      n2, n2x, sigma){
  n2_z1 <-ifelse(z1<= cf | z1>ce, n1,
                 ifelse((z1 > cf & z1 < zl )| (z1 < ce & z1 > zu ), n2,
                        ifelse((zl < z1 & z1 < zc), n2x, 
                               get.n2star.pp(z1, target.pp, m0, s0,
                                             n1, sigma, 
                                             n2, n2x,  z.alpha, ce, cf)))) 
  return(n2_z1)
  
}



## obtain boundaries, sample size rule and conditional power for the PPZ approach

get.ssr.bdry.ppz<-function(z1, ppmin, ppmax, m0, s0, n1, n2, k, ce, cf, z.alpha, sigma){
  # z <- z1
  n2max <- n2 * k
  
  z.l <- get.z1.pp(target.pp=ppmin, m0, s0, n1, sigma, 
                 n2, n2x=n2max, z.alpha, ce, cf)
  z.u <- get.z1.pp(target.pp=ppmax, m0, s0, n1, sigma, 
                  n2, n2x=n2, z.alpha, ce, cf)
  z.c <- get.z1.pp(target.pp=ppmax, m0, s0, n1, sigma, 
                  n2, n2x=n2max, z.alpha, ce, cf)
  
  ss <-  sapply(z1, ss.ppz.function, zl=z.l, zc=z.c, zu=z.u, 
                ce=ce, cf=cf, z.alpha=z.alpha, target.pp=ppmax,
                m0=m0, s0=s0, n1=n1, n2=n2, n2x=n2max, 
                sigma=sigma)
  ## Added
  pps<- mapply(FUN = pp, z1, n2x=ss, m0=m0, s0=s0, 
               n1=n1, sigma=sigma, n2=n2, 
               z.alpha=z.alpha, ce=ce, cf=cf)
  
  zone <-  sapply(z1, get.zone, zl=z.l, zc=z.c, zu=z.u, ce=ce, cf=cf)
  
  promising_bdry<- c(z.l=z.l, z.c=z.c, z.u=z.u)
  early_stop_bdry <- c(cf=cf, ce=ce)
  
  SSR.df <- data.frame(z1, n2_z1=ss, 
                       PredictivePower=pps,Zone=zone)
  
  
  
  return( ssr=list(early_stop_bdry = early_stop_bdry, 
                   promising_bdry = promising_bdry, SSR.df = SSR.df ))
  
}


## obtain boundaries, sample size rule,conditional power and zonewise summaries 
## for the PPZ approach
SSR.Power.Zonewise.PPZ<- function(ppmin, ppmax, m0, s0, n1, n2, k, ce, cf, z.alpha,
                      delta.true,  sigma,
                      n.sim=10000, seed=NULL){
  
  mean.z1 <- delta.true * sqrt(n1)/(2 * sigma)
  
  z1 <- rnorm(n.sim, mean.z1, sd=1)
  
  ssr=get.ssr.bdry.ppz(z1, ppmin=ppmin, ppmax=ppmax, m0=m0, s0=s0,
                      n1=n1, n2=n2, k=k, ce=ce, cf=cf, z.alpha=z.alpha, sigma=sigma)
  
  n2star<-ssr$SSR.df$n2_z1
  
  z.l <- ssr$promising_bdry["z.l"]
  z.c <- ssr$promising_bdry["z.c"]
  z.u <- ssr$promising_bdry["z.u"]
  
  pow <- get.power(delta.true=delta.true, z1,  zl=z.l, zc=z.c, zu=z.u, ce=ce, cf=cf, z.alpha=z.alpha,
                   n1=n1, n2=n2, n2star=n2star, sigma=sigma)
  
  
  ssr.power.df <- cbind(ssr$SSR.df, Success=pow$Success)
  
  zonewise <- ssr.power.df %>% group_by(Zone)%>%
    summarize(Prob_Zone = n()/n.sim*100,
              Average_SS = mean(n2_z1),
              Power = mean(Success)*100,
              Count_simulations=n())
  
  zonewise.df <- data.frame(True_delta = delta.true, zonewise )
  
  all.trials<- data.frame(True_delta = delta.true, Zone="All Trials",
                          Prob_Zone = length(z1)/n.sim*100,
                 Average_SS = mean(ssr.power.df$n2_z1),
                 Power = mean(ssr.power.df$Success)*100,
                 Count_simulations=length(z1))

  zonewise_summaries<- rbind(zonewise.df, all.trials )

  return(list(ssr.power.df=ssr.power.df, zonewise_summaries=zonewise_summaries,
            boundaries=data.frame(cf=cf, z.l=z.l, z.c=z.c, z.u=z.u, ce=ce)))
  
  # return(list(ssr.power.df=ssr.power.df, zonewise_summaries=zonewise,
  #             boundaries=data.frame(cf=cf, z.l=z.l, z.c=z.c, z.u=z.u, ce=ce))) 
}

## pow.ss.diff.deltas.ppz is function used to obtain operating characteristics-power and sample size at
## different values of true delta/ response generation delta for the PPZ approach

pow.ss.diff.deltas.ppz <-function(delta.true, zl, zc, zu, n1, n2, k, ce, cf, z.alpha, 
                           m0, s0, ppmin, ppmax, n.sim=10000, seed=NULL){
  
  
  mean.z1 <- delta.true * sqrt(n1)/(2 * sigma)
  
  z1 <- rnorm(n.sim, mean.z1, sd=1)
  
  ss <-  sapply(z1, ss.ppz.function, zl=zl, zc=zc, zu=zu, 
                ce=ce, cf=cf, z.alpha=z.alpha, target.pp=ppmax,
                m0=m0, s0=s0, n1=n1, n2=n2, n2x=n2max, 
                sigma=sigma)
  
  pow <- get.power(delta.true=delta.true, z1,  zl=zl, zc=zc, zu=zu,  ce=ce, cf=cf, z.alpha=z.alpha,
                   n1=n1, n2=n2, n2star=ss, sigma=sigma)
  
  zone <- get.zone(z1,  zl=zl, zc=zc, zu=zu, ce=ce, cf=cf)
  
  pow.ss.df <- cbind(True_delta=delta.true, pow, zone)
  
  return(pow.ss.df)
}


## operating characteristics- power and sample size at 
## different values of true delta/ response generation delta
OC.Power.SS.deltas.PPZ<- function(ppmin, ppmax, m0, s0, n1, n2, k, ce, cf, z.alpha,
                                   delta.vec, sigma,
                                   n.sim=10000, seed=NULL){
  n2max<-n2*k
  z.l <- get.z1.pp(target.pp=ppmin, m0, s0, n1, sigma, 
                  n2, n2x=n2max, z.alpha, ce, cf)
  z.u <- get.z1.pp(target.pp=ppmax, m0, s0, n1, sigma, 
                  n2, n2x=n2, z.alpha, ce, cf)
  z.c <- get.z1.pp(target.pp=ppmax, m0, s0, n1, sigma, 
                  n2, n2x=n2max, z.alpha, ce, cf)
  
  out <-bind_rows(lapply(delta.vec,  pow.ss.diff.deltas.ppz,  zl=z.l, zc=z.c, zu=z.u, 
                         n1, n2, k, ce, cf, z.alpha, 
                             m0, s0, ppmin, ppmax, 
                             n.sim, seed))
  
  delta_wise_OC <- out %>% group_by(True_delta)%>%
    summarize(Average_SS = mean(n2star),
              Power = mean(Success)*100,
              Total_simulations=n())

  
  
  return(list(delta_wise_OC=delta_wise_OC,
              boundaries=data.frame(cf=cf, z.l=z.l, z.c=z.c, z.u=z.u, ce=ce)))
}


