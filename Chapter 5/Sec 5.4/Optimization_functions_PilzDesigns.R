library(tidyverse)
library(doSNOW)
library(foreach)
library("doFuture")
registerDoFuture()
plan(multicore, workers = 8)

## Definitions:
# n1: stage 1 sample size
# n2: initially planned final sample size
# n2_z1: new/adapted sample size as a function of interim test statistic z1
# k: sample size multiplier 
# ce: early efficacy stopping boundary
# cf: early futility stopping boundary
# c2: critical bondary at the final look
# sigma: design standard deviation
# delta0: delta under the null hypothesis
# delta.target: target delta value (equivalent to delta min in CPZ method) used to find the sample size rule
# delta.true: true delta/response generation delta
# target.power: target power to find lambda at the target delta






##-------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------

# get z statistics
get.z1<-function(delta,delta0, n1, sigma, seed,n.sim){
  set.seed(seed)
  
  delta_st<-(delta-delta0)/(2*sigma)
  z.rs<-rnorm(n.sim,sqrt(n1)*delta_st,1)
  return(z.rs)
}

##-------------------------------------------------------------------------------------------------------------
# get c2 boundary

Az <- function (z, c2, n, n1, n2) 
{
  (1-pnorm((c2*sqrt(n)-z*sqrt(n1))/sqrt(n2))) * dnorm(z)
}

up_prob<-function(c2, cf, ce, n, n1, n2){
  integrate(Az , lower=cf, upper = ce, c2=c2, n=n, n1=n1, n2=n2 )$value
}


get.c2<-function(n,n1,n2,ce,cf, t1e)
{
  
  fn_c2<-function(cf, ce, c2, n, n1, n2, t1e)
    
  {
    (1-pnorm(ce)) + up_prob(c2, cf,ce, n, n1, n2) - t1e
  }
  
  c2.sol<-uniroot(fn_c2,c(0.000001 ,100000),cf=cf, ce=ce, n=n, n1=n1,
                  n2=n2, t1e=t1e)$root
  return(c2.sol)
}



##-------------------------------------------------------------------------------------------------------------
## probability in the continuation region

Cond_Prob<-function(n1, n2, c2, z, n2_z1, delta.true, delta0, sigma){
  del_st<-(delta.true-delta0)/(2*sigma)
  cont.prob<-1-pnorm((sqrt(n2)*c2-sqrt(n1)*z)/sqrt(n2-n1)-sqrt(n2_z1 - n1)*del_st)
  return(cont.prob)
}

##-------------------------------------------------------------------------------------------------------------
## n2(z1) for no efficacy and futility boundaries 
get.n2star<-function(z1,lambda,delta.target,delta0,delta.true,sigma,c2,n1,n2,k)
{                    
  
  n2max<-n2*k
  delta_st<-(delta.target-delta0)/(2*sigma)
  
  
  fn2_star.incr<-function(n2star.incr)
  {
    2*log(sqrt(8*pi)/lambda)-2*log(delta_st)+log(n2star.incr)+
      ((sqrt(n2)*c2-sqrt(n1)*z1)/sqrt(n2-n1)-sqrt(n2star.incr)*delta_st)^2
  }
  
  
  n2star.incr<-tryCatch(
    {
      uniroot(fn2_star.incr,c(0.0001,10000))$root
    },
    
    error=function(e) {
      e
      
      n2-n1                    ## we set the value of incremental sample size to n2 
      ## if root of the equation doesn't exist    
    },
    warning=function(w)
    {
      w
      NA
    }
  )
  
  n2star <- n1 + n2star.incr
  n2star.final <- ifelse(n2star< n2, n2 , ifelse(n2star>n2max, n2max, n2star)) 
  n2star.incr.final<-n2star.final-n1
  CondProb<-Cond_Prob(n1, n2, c2, z=z1, n2_z1=n2star.final,
                              delta.true,delta0,sigma)
  
  n2.ret.df<-data.frame(delta.target, delta.true, lambda, z1, n2star.incr, n2star.incr.final,
                        n2star, n2star.final,CondProb)
  return(n2.ret.df)
}

##-------------------------------------------------------------------------------------------------------------

### n2(z1) for with efficacy and user specified futility boundaries 
get.n2star.stop<-function(z1, lambda, delta.target, delta0, delta.true,
                          sigma, c2, ce, cf, n1, n2, k) 
{
  
  n2max <- n2*k
  delta_st <- (delta.target-delta0)/(2*sigma)
  
   fn2_star.incr<-function(n2star.incr)
     {
          2*log(sqrt(8*pi)/lambda)-2*log(delta_st)+log(n2star.incr)+
          ((sqrt(n2)*c2-sqrt(n1)*z1)/sqrt(n2-n1)-sqrt(n2star.incr)*delta_st)^2
     }
  
  
  if(cf < z1 & z1 <ce){
    
    n2star.incr<-tryCatch(
      {
        uniroot(fn2_star.incr,c(0.001,10000), tol=0.0001)$root
      },
      
      error=function(e) {
        e
        
        n2-n1                    ## we set the value of incremental sample size to n2 
        ## if root of the equation doesn't exist    
      },
      warning=function(w)
      {
        w
        NA
      }
    )
    
    n2star <- n1+ n2star.incr
    n2star.final<- ifelse(n2star< n2, n2 , ifelse(n2star>n2max, n2max, n2star)) 
    #nstar.final<-  ifelse(nstar>nmax, nmax, nstar)
    
    n2star.incr.final<-n2star.final-n1
    
    CondProb<-Cond_Prob(n1=n1, n2=n2, c2=c2, z=z1, n2_z1=n2star.final,
                                delta.true=delta.true, delta0=delta0, sigma=sigma)
    
    
    n2.ret.df<-data.frame(delta.target,delta.true,lambda,z1, n2star.incr, n2star.incr.final, n2star,
                          n2star.final, CondProb)
    
    return(n2.ret.df)
  }else
   {
    n2star.incr.final<-n2star.incr<-0
    n2star.final<-n2star<-n1+n2star.incr
    CondProb<-ifelse(z1<=cf,0,1)
    
    n2.ret.df<-data.frame(delta.target,delta.true,lambda,z1, n2star.incr, n2star.incr.final, n2star,
                          n2star.final, CondProb)
    return(n2.ret.df)
  }
}


##-------------------------------------------------------------------------------------------------------------


## Find lambda for design with no early efficacy and futility boundary without CP constraint
get.lambda.no.stop<-function(target.power, c2, n1, n2, k, delta0, delta.target,
                             sigma,n.sim, seed)
{
  
  n2max <- k*n2
  delta0<-delta0
  
  z <-get.z1(delta=delta.target, delta0, n1, sigma, seed, n.sim)
  
  lambda.vec0<- seq(10000,50000,10000)
  SS.df0<-foreach(i = 1:length(lambda.vec0), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                   bind_rows(lapply(z, get.n2star, lambda=lambda.vec0[i],
                                  delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                                  sigma=sigma, c2=c2, n1=n1, n2=n2, k=k))}

  
  lambda.CP.ESS.df0<-SS.df0%>%group_by(lambda)%>%
                     summarize(ESS=mean(n2star.final),
                               Overall.Power=mean(CondProb))
  
  sel.lam0<-lambda.CP.ESS.df0%>%data.frame()%>%
            slice(which.min(abs(Overall.Power-target.power)))%>%
            select(lambda)%>%as.numeric()
  
  
  lambda.vec1<- seq(sel.lam0-10000, sel.lam0+10000, 1000)[-1]
  SS.df1<-foreach(i = 1:length(lambda.vec1), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z, get.n2star, lambda=lambda.vec1[i],
                                     delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                                     sigma=sigma, c2=c2, n1=n1, n2=n2, k=k))}
  
  
  lambda.CP.ESS.df1<-SS.df1%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final),
              Overall.Power=mean(CondProb))
  
  sel.lam1<-lambda.CP.ESS.df1%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  lambda.vec2<- seq(sel.lam1-1000, sel.lam1+1000, 100)[-1]
  SS.df2<-foreach(i = 1:length(lambda.vec2), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z, get.n2star, lambda=lambda.vec2[i],
                                     delta.target=delta.target, delta0=delta0,
                                     delta.true=delta.target,
                                     sigma=sigma, c2=c2, n1=n1, n2=n2, k=k))}
  
  
  lambda.CP.ESS.df2<-SS.df2%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final),
              Overall.Power=mean(CondProb))
  
  sel.lam2<-lambda.CP.ESS.df2%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  lambda.vec3<- seq(sel.lam2-100, sel.lam1+100, 10)[-1]
  SS.df3<-foreach(i = 1:length(lambda.vec3), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z, get.n2star, lambda=lambda.vec3[i],
                                     delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                                     sigma=sigma, c2=c2, n1=n1, n2=n2, k=k))}
  
  
  lambda.CP.ESS.df3<-SS.df3%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final),
              Overall.Power=mean(CondProb))
  
  sel.lam3<-lambda.CP.ESS.df3%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  lambda.vec4<- seq(sel.lam3-10, sel.lam3+10, 1)[-1]
  SS.df4<-foreach(i = 1:length(lambda.vec4), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z, get.n2star, lambda=lambda.vec4[i],
                                     delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                                     sigma=sigma, c2=c2, n1=n1, n2=n2, k=k))}
  
  
  lambda.CP.ESS.df4<-SS.df4%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final),
              Overall.Power=mean(CondProb))
  
  sel.lam4<-lambda.CP.ESS.df4%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  
  lambda.vec5<- seq(sel.lam4-1, sel.lam4+1, 0.1)[-1]
  SS.df5<-foreach(i = 1:length(lambda.vec5), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z, get.n2star, lambda=lambda.vec5[i],
                                     delta.target=delta.target, delta0=delta0, 
                                     delta.true=delta.target,
                                     sigma=sigma, c2=c2, n1=n1, n2=n2, k=k))}
  
  
  lambda.CP.ESS.df5<-SS.df5%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final),
              Overall.Power=mean(CondProb))
  
  sel.lam5<-lambda.CP.ESS.df5%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  lambda.vec6<- seq(sel.lam4-0.1, sel.lam4+0.1, 0.01)[-1]
  SS.df6<-foreach(i = 1:length(lambda.vec6), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z, get.n2star, lambda=lambda.vec6[i],
                                     delta.target=delta.target, delta0=delta0, 
                                     delta.true=delta.target,
                                     sigma=sigma, n=n, c2=c2, n1=n1, n2=n2, k=k))}
  
  
  lambda.CP.ESS.df6<-SS.df6%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final),
              Overall.Power=mean(CondProb))
  
  sel.lam6<-lambda.CP.ESS.df6%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  
  return(list(opt.lambda=sel.lam6,lambda.Power.ESS=lambda.CP.ESS.df6,
              SS.CP.df=SS.df6))
}

##-------------------------------------------------------------------------------------------------------------

## Find lambda for design with early efficacy and futility boundary without CP constraint
get.lambda.stop<-function(target.power, cf, ce, c2, n1, n2, k,
                          delta0, delta.target,
                          sigma, n.sim, seed) {
  
  
  z1 <- get.z1(delta=delta.target, delta0, n1, sigma, seed, n.sim)
  
  lambda.vec0<- seq(10000,50000,10000)
  
  SS.df0<-foreach(i = 1:length(lambda.vec0), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                  bind_rows(lapply(z1, get.n2star.stop, lambda=lambda.vec0[i],
                            delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                            sigma=sigma, c2=c2, ce=ce, cf=cf, 
                            n1=n1, n2=n2, k=k))}

  
  lambda.CP.ESS.df0<-SS.df0%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final), 
              Overall.Power=mean(CondProb))
  
   sel.lam0<-lambda.CP.ESS.df0%>%data.frame()%>%
             slice(which.min(abs(Overall.Power-target.power)))%>%
             select(lambda)%>%as.numeric()
  
  lambda.vec1<-seq(sel.lam0-10000, sel.lam0+10000, 1000)[-1]
 
  SS.df1<-foreach(i = 1:length(lambda.vec1), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                   bind_rows(lapply(z1, get.n2star.stop,lambda=lambda.vec1[i],
                             delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                             sigma=sigma, c2=c2, ce=ce, cf=cf,
                             n1=n1, n2=n2, k=k))}
  
  
  lambda.CP.ESS.df1<-SS.df1%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final), 
              Overall.Power=mean(CondProb))
  
  sel.lam1<-lambda.CP.ESS.df1%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  lambda.vec2<-seq(sel.lam1-1000, sel.lam1+1000, 100)[-1]
  
  SS.df2<-foreach(i = 1:length(lambda.vec2), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z1, get.n2star.stop,lambda=lambda.vec2[i],
                              delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                              sigma=sigma, c2=c2, ce=ce, cf=cf,
                              n1=n1, n2=n2, k=k))}
  
  lambda.CP.ESS.df2<-SS.df2%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final), 
              Overall.Power=mean(CondProb))
  
  sel.lam2<-lambda.CP.ESS.df2%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  lambda.vec3<-seq(sel.lam2-100, sel.lam2+100, 10)[-1]
  
  SS.df3<-foreach(i = 1:length(lambda.vec3), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z1, get.n2star.stop,lambda=lambda.vec3[i],
                              delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                              sigma=sigma, c2=c2, ce=ce, cf=cf,
                              n1=n1, n2=n2, k=k))}

  lambda.CP.ESS.df3<-SS.df3%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final), 
              Overall.Power=mean(CondProb))
  
  sel.lam3<-lambda.CP.ESS.df3%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  lambda.vec4<-seq(sel.lam3-10, sel.lam3+10, 1)[-1]
  
  SS.df4<-foreach(i = 1:length(lambda.vec4), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z1, get.n2star.stop,lambda=lambda.vec4[i],
                              delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                              sigma=sigma, c2=c2, ce=ce, cf=cf,
                              n1=n1, n2=n2, k=k))}
  
  
  lambda.CP.ESS.df4<-SS.df4%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final), 
              Overall.Power=mean(CondProb))
  
  sel.lam4<-lambda.CP.ESS.df4%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  lambda.vec5<-seq(sel.lam4-1, sel.lam4+1, 0.1)[-1]
  
  SS.df5<-foreach(i = 1:length(lambda.vec5), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z1, get.n2star.stop,lambda=lambda.vec5[i],
                              delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                              sigma=sigma, c2=c2, ce=ce, cf=cf,
                              n1=n1, n2=n2, k=k))}
  

  
  lambda.CP.ESS.df5<-SS.df5%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final), 
              Overall.Power=mean(CondProb))
  
  
  sel.lam5<-lambda.CP.ESS.df5%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  lambda.vec6<-seq(sel.lam5-0.1, sel.lam5+0.1, 0.02)[-1]
  
  SS.df6<-foreach(i = 1:length(lambda.vec6), .combine = "bind_rows" ,
                  .packages = "tidyverse") %dopar% { 
                    bind_rows(lapply(z1, get.n2star.stop,lambda=lambda.vec6[i],
                              delta.target=delta.target, delta0=delta0, delta.true=delta.target,
                              sigma=sigma, c2=c2, ce=ce, cf=cf,
                              n1=n1, n2=n2, k=k))}

  lambda.CP.ESS.df6<-SS.df6%>%group_by(lambda)%>%
    summarize(ESS=mean(n2star.final), 
              Overall.Power=mean(CondProb))
  
  sel.lam6<-lambda.CP.ESS.df6%>%data.frame()%>%
    slice(which.min(abs(Overall.Power-target.power)))%>%
    select(lambda)%>%as.numeric()
  
  
  return(list(opt.lambda=sel.lam6,lambda.Power.ESS=lambda.CP.ESS.df6,
              SS.CP.df=SS.df6))
}

##-------------------------------------------------------------------------------------------------------------
## Find optimal design without CP constraint
## Find optimal lambda at delta_min/target_delta explore operating characteristics over a series of deltas

find.unconstrained.design<-function(delta0, delta.target, delta.true, target.power,
                                    sigma, n1, n2, k, ce=NULL, cf=NULL, c2, n.sim, seed=NULL) {
  
  start_time <- Sys.time()
  
  delta.vec<-delta.true
 
    if(is.null(ce) & is.null(cf)){
    
    lam.opt<-get.lambda.no.stop(target.power, c2, n1, n2, k, delta0, delta.target,
                                 sigma, n.sim, seed)


    lambda<-lam.opt$opt.lambda
    
    SS.df<- foreach(i = 1:length(delta.vec), .combine = "bind_rows" ,
                    .packages = "tidyverse") %dopar% {
                      
                      z=get.z1(delta=delta.vec[i],delta0, n1, sigma, seed,n.sim=10000)
                      
                      bind_rows(lapply(z, get.n2star,lambda=lambda, delta.target=delta.target,
                                       delta0, delta.true=delta.vec[i], sigma, c2, n1, n2, k))  } 
    
    delta.CP.ESS.df<-SS.df%>%group_by(delta.true)%>%
      summarize(ESS=mean(n2star.final),
                Overall.Power=mean(CondProb))
    
    
    names(delta.CP.ESS.df)<-c("delta","ESS","Power")
    
    ## Expected Sample Size vs delta curve
    SS.Plot1<-ggplot()+
      xlab(expression(delta)) + 
      ylab("Expected Sample Size") +
      geom_line(data = delta.CP.ESS.df, aes(x = delta, y = ESS), color="blue", size=1)+
      geom_vline(xintercept =seq(min(delta.vec), max(delta.vec), 0.1), col = "gray")+
      geom_hline(yintercept =seq(roundUp10(n)-50, roundUp10(nmax)+50, 50), col = "gray")+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    
    ## Power vs delta curve
    SS.Plot2<-ggplot()+
      xlab(expression(delta)) + 
      ylab("Power") +
      geom_line(data = delta.CP.ESS.df, aes(x = delta, y = Power), color="blue", size=1)+
      geom_vline(xintercept =seq(min(delta.vec), max(delta.vec), 0.1), col = "gray")+
      geom_hline(yintercept =seq(0,1,0.1), col = "gray")+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    SSR.df<-SS.df%>%filter(delta.true==delta.target)
    
    ## Sample Size Rule Plot
    SS.Plot3<-ggplot()+
      xlab("Z1 at the interim") + 
      ylab("Sample Size Rule")  +
      geom_vline(xintercept =seq(min(SSR.df$z1), max(SSR.df$z1), 0.5), col = "gray")+
      geom_hline(yintercept =seq(roundUp10(n)-10, roundUp10(nmax)+10, 20),
                 col = "gray")+
      geom_line(data = SSR.df, aes(x = z1, y = nstar.final), color="blue", size=1)+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    ## Sample Size Rule Plot
    SS.Plot4<-ggplot()+
      xlab("Z1 at the interim") + 
      ylab("Conditional Power")  +
      geom_vline(xintercept =seq(min(SSR.df$z1), max(SSR.df$z1), 0.5), col = "gray")+
      geom_hline(yintercept =seq(0,1,0.1),
                 col = "gray")+
      geom_line(data = SSR.df, aes(x = z1, y = CondProb), color="blue", size=1)+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    end_time <- Sys.time()
    elapsed.time<-end_time - start_time
    
    return(list(opt.lambda=lam.opt$opt.lambda, 
                lambda.Power.ESS=lam.opt$lambda.Power.ESS,
                SS.CP.all.lambda=lam.opt$SS.CP.df,
                delta.Power.ESS=delta.CP.ESS.df,
                SS.CP.all.delta=SS.df, 
                delta.vs.ESS.plot=SS.Plot1,
                delta.vs.Power.plot=SS.Plot2,
                SSR.plot=SS.Plot3,
                Cond.Power.plot=SS.Plot4,
                elapsed.time=elapsed.time))
    
    
  }
    else{
    
   lam.opt<-get.lambda.stop(target.power, cf, ce, c2, n1, n2, k, delta0,
                             delta.target, sigma, n.sim, seed)
    
   lambda.opt<-lam.opt$opt.lambda
    
    delta.vec<-delta.true
    
    SS.df<- foreach(i = 1:length(delta.vec), .combine = "bind_rows" ,
                    .packages = "tidyverse") %dopar% { 
                      
                      z=get.z1(delta=delta.vec[i],delta0, n1, sigma, seed,n.sim=n.sim)
                      
                      bind_rows(lapply(z, get.n2star.stop, lambda=lambda.opt, delta.target=delta.target,
                                       delta0=delta0, delta.true=delta.vec[i], sigma=sigma,
                                       c2=c2, ce=ce, cf=cf, n1=n1, n2=n2, k=k))}
    
    delta.CP.ESS.df<-SS.df%>%group_by(delta.true)%>%
      summarize(ESS=mean(n2star.final),
                Overall.Power=mean(CondProb))
    
    names(delta.CP.ESS.df)<-c("delta","ESS","Power")
    
    ## Expected Sample Size vs delta curve
    SS.Plot1<-ggplot()+
      xlab(expression(delta)) + 
      ylab("Expected Sample Size") +
      geom_line(data = delta.CP.ESS.df, aes(x = delta, y = ESS), color="blue", size=1)+
      geom_vline(xintercept =seq(min(delta.vec), max(delta.vec), 0.1), col = "gray")+
      geom_hline(yintercept =seq(roundUp10(n)-50, roundUp10(nmax)+50, 50), col = "gray")+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    ## Power vs delta curve
    SS.Plot2<-ggplot()+
      xlab(expression(delta)) + 
      ylab("Power") +
      geom_line(data = delta.CP.ESS.df, aes(x = delta, y = Power), color="blue", size=1)+
      geom_vline(xintercept =seq(min(delta.vec), max(delta.vec), 0.1), col = "gray")+
      geom_hline(yintercept =seq(0,1,0.1), col = "gray")+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    SSR.df<-SS.df%>%filter(delta.true==delta.target)%>%
      mutate(CondProb=ifelse(z1<cf, 1, ifelse(ce>z1,0,CondProb)))
    
    ## Sample Size Rule Plot
    SS.Plot3<-ggplot()+
      xlab("Z1 at the interim") + 
      ylab("Sample Size Rule")  +
      geom_vline(xintercept =seq(min(SSR.df$z1), max(SSR.df$z1), 0.5), col = "gray")+
      geom_hline(yintercept =seq(roundUp10(n1)-10, roundUp10(nmax)+10, 20),
                 col = "gray")+
      geom_line(data = SSR.df, aes(x = z1, y = nstar.final), color="blue", size=1)+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    ## Sample Size Rule Plot
    SS.Plot4<-ggplot()+
      xlab("Z1 at the interim") + 
      ylab("Conditional Power")  +
      geom_vline(xintercept =seq(min(SSR.df$z1), max(SSR.df$z1), 0.5), col = "gray")+
      geom_hline(yintercept =seq(0,1,0.1),
                 col = "gray")+
      geom_line(data = SSR.df, aes(x = z1, y = CondProb), color="blue", size=1)+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    end_time <- Sys.time()
    elapsed.time<-end_time - start_time
    
    return(list(opt.lambda=lam.opt$opt.lambda, 
                lambda.Power.ESS=lam.opt$lambda.Power.ESS,
                SS.CP.all.lambda=lam.opt$SS.CP.df,
                delta.Power.ESS=delta.CP.ESS.df,
                SS.CP.all.delta=SS.df, 
                delta.vs.ESS.plot=SS.Plot1,
                delta.vs.Power.plot=SS.Plot2,
                SSR.plot=SS.Plot3,
                Cond.Power.plot=SS.Plot4,
                elapsed.time=elapsed.time))
    
    
  }
  
}

