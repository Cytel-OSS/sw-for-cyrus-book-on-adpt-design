
# Functions used to generate results in the manuscript:
# Hsiao S, Liu L, Mehta C, Optimal Promising Zone Designs, Biometrical Journal 2018

#rm(list = ls())

# -----------------------------------------------------------------------------
# Definitions of arguments used in the functions
#
#  n1 = stage 1 sample size
#  n2 = planned cumulative final sample size at end of stage 2 without adaptation
#  n2x = Actual cumulative sample size at end of stage 2
#  z1 = stage 1 Z-statistic
#  w1 and w2 = CHW weights
#  theta = treatment effect mean
#  sigma = variance of endpoint
#  theta_cp = reference theta for conditional power calculations
#  * if theta_cp = NULL, calculate CP using estimated theta *
#  bdry = final efficacy boundary
#  c1 = interim efficacy boundary
#  z1min and z1max = lower and upper bound of promising zone
#  targetcp = CP to target for sloped SSR rule
#  cpmin = minimum CP that defines start of promising zone
#  cpmax = maximum CP that defines end of promising zone
#  search_z1_min_jt, search_z1_max_jt = lower and upper limits of search grid for z1
#  search_z1_increment_jt = increments to use in search grid for z1
#  mc = number of sample points to use for Monte Carlo integration
#  revenue_bound = lower bound on revenue function (e.g., CP) for contrained optimization
#  theta_revenue_bound = theta value at which revenue constraint is evaluated (e.g., CP evaluated at this theta value)

# -----------------------------------------------------------------------------
# conditional power (CP) function
cp <- function(n1, n2x, z1, bdry, w1, w2, theta_cp = NULL, sigma , ce = Inf, cf= -Inf) {
  if ( cf < z1 & z1 < ce) {
    if (is.null(theta_cp)) {
      pnorm((w1/w2 + sqrt((n2x - n1)/n1))*z1 - bdry/w2)
    } else {
      pnorm(theta_cp*sqrt(n2x - n1)/(2*sigma) + (w1*z1 - bdry)/w2)
    }
  } else {
    ifelse(z1 > ce, 1.0, 0.0 )
    
  }
}

# convert cp to z1
cp_to_z1 <- function(n1, n2x, cp, bdry, w1, w2, theta_cp = NULL, sigma) {
  if (is.null(theta_cp)) {
    (w2*qnorm(cp) + bdry)/(w1 + w2*sqrt((n2x - n1)/n1))
  } else {
    (w2*qnorm(cp) + bdry - w2*theta_cp*sqrt(n2x - n1)/(2*sigma))/w1
  }
}

# lower promising zone where SSR begins, required to achieve
# .. conditional power = cp_bound at theta_cp
promzone_min_z1 <- function(n1, n2x, bdry, w1, w2, theta_cp, sigma, cp_bound) {
  (w2*qnorm(cp_bound) - w2*theta_cp*sqrt(n2x - n1)/(2*sigma) + bdry)/w1
}

# Expectation of z1
expected_z1 <- function(theta, sigma, n1) {
  theta*sqrt(n1)/(2*sigma)
}


# Density function of z1
density_z1 <- function(z1, theta, sigma, n1) {
  dnorm(z1 - expected_z1(theta = theta, sigma = sigma, n1 = n1))
}


# Probability of being in promising zone
prob_prom <- function(z1min, z1max, theta, sigma, n1) {
  Ez1 <- expected_z1(theta = theta, sigma = sigma, n1 = n1)
  pnorm(z1max, mean = Ez1, sd = 1) - pnorm(z1min, mean = Ez1, sd = 1)
}

# Probability of being in the three zones: unfav, prom, fav
prob_zones <- function(z1min, z1max, ce, cf, theta, sigma, n1) {
  Ez1 <- expected_z1(theta = theta, sigma = sigma, n1 = n1)
  
  if (ce < z1max) {
    stop("Early efficacy boundary less than upper promising zone boundary")
  }else{
    if(cf > z1min){
      stop("Early futility boundary greater than lower promising zone boundary") 
    }
  }
  
  list(fut = pnorm(cf, mean = Ez1, sd = 1),
       unfav = pnorm(z1min, mean = Ez1, sd = 1)-pnorm(cf, mean = Ez1, sd = 1),
       prom = pnorm(z1max, mean = Ez1, sd = 1) - pnorm(z1min, mean = Ez1, sd = 1),
       fav = (1 - pnorm(z1max, mean = Ez1, sd = 1))-(1 - pnorm(ce, mean = Ez1, sd = 1)),
       eff = 1 - pnorm(ce, mean = Ez1, sd = 1))
}


# Expected sample size using JT method
avgn_jt2 <- function(n1, n2, n2min, n2max, bdry, w1, w2, ce, cf, theta, sigma,
                     theta_tilde, gam,
                     revenue, cost, objective, stepsize = 0.1,
                     search_z1_min_jt, search_z1_max_jt,
                     search_z1_increment_jt, mc = 1e5,
                     revenue_bound = NULL, theta_revenue_bound = NULL) {
  # expected values of z1
  Ez1 <- expected_z1(theta = theta, sigma = sigma, n1 = n1)
  # density of z1
  fz1 <- function(z1) {
    dnorm(z1 - Ez1)
  }

  # SSR region boundaries on z scale
  z1min <- find_z1min_jt(n1 = n1, n2 = n2, n2min = n2min, n2max = n2max, bdry = bdry,
                         w1 = w1, w2 = w2, theta_tilde = theta_tilde,
                         sigma = sigma, ce=ce, cf=cf, gam = gam, revenue = revenue, cost = cost,
                         objective = objective, stepsize = stepsize,
                         search_z1_min_jt = search_z1_min_jt,
                         search_z1_max_jt = search_z1_max_jt,
                         search_z1_increment_jt = search_z1_increment_jt, 
                         revenue_bound = revenue_bound,
                         theta_revenue_bound = theta_revenue_bound)
  z1max <- find_z1max_jt(n1 = n1, n2 = n2, n2min = n2min, n2max = n2max, bdry = bdry,
                         w1 = w1, w2 = w2, theta_tilde = theta_tilde,
                         sigma = sigma, ce=ce, cf=cf, gam = gam, revenue = revenue, cost = cost,
                         objective = objective, stepsize = stepsize,
                         search_z1_min_jt = search_z1_min_jt,
                         search_z1_max_jt = search_z1_max_jt,
                         search_z1_increment_jt = search_z1_increment_jt, revenue_bound = revenue_bound,
                         theta_revenue_bound = theta_revenue_bound)
  z1grid <- seq(z1min, z1max, lookup_table_increments)
  n2x_lookup_table <- n2x_jt(n1 = n1, n2 = n2, n2min = n2min, n2max = n2max, z1 = z1grid,
                             bdry = bdry, w1 = w1, w2 = w2,
                             theta_tilde = theta_tilde, sigma = sigma, ce=ce, cf=cf, 
                             gam = gam, revenue = revenue, cost = cost,
                             objective = objective, stepsize = stepsize,
                             revenue_bound = revenue_bound, theta_revenue_bound = theta_revenue_bound)
  # faster computation of n2 star
  n2x_jt_fast <- function(z1) {
    out <- numeric(length(z1))
    ssr_region <- (z1 > z1min & z1 < z1max & z1 < ce & z1 > cf)
    out[z1 > ce | z1 < cf] <- n1
    out[!ssr_region & z1 < ce & z1 > cf] <- n2
    
    out[ssr_region] <- n2x_lookup_table[closest_z1(z1[ssr_region], z1grid,
                                                   lookup_table_increments)]
   
    out
  }

  # monte carlo integration - average n2x over randomly sampled z1 values
  mean(n2x_jt_fast(rnorm(mc, mean = Ez1, sd = 1)))
}


# Expected sample size of CPZ design using sloped SSR rule
avgn_ssrslope <- function(n1, n2, n2min, n2max, bdry, w1, w2, ce, cf,
                         theta, sigma ,
                         theta_tilde = NULL, targetcp,
                         cpmin, cpmax) {
  # SSR region boundaries on z scale
  z1min <- cp_to_z1(n1 = n1, n2x = n2max, cp = cpmin, bdry = bdry,
                    w1 = w1, w2 = w2, theta_cp = theta_tilde, sigma = sigma)
  z1max <- cp_to_z1(n1 = n1, n2x = n2, cp = cpmax, bdry = bdry,
                    w1 = w1, w2 = w2, theta_cp = theta_tilde, sigma = sigma)
  # expected values of z1
  Ez1 <- expected_z1(theta = theta, sigma = sigma, n1 = n1)

  # if (ce < z1max) {
  #   stop("Early efficacy boundary less than upper promising zone boundary")
  # }else{
  #       if(cf > z1min){
  #         stop("Early futility boundary greater than lower promising zone boundary") 
  #       }
  #   
  # }

  # probability that final sample size is n1 (early stop for efficacy)
  prob_n1 <- pnorm(Ez1 - ce)+(1-pnorm(Ez1 - cf))
  
  
  # probability that final sample size is n2
  
  if(z1min<cf & z1max <ce){
    prob_n2 <- pnorm(ce, mean = Ez1, sd = 1) - pnorm(z1max, mean = Ez1, sd = 1)
  }else{
    if(z1min > cf & z1max > ce){
      prob_n2 <- pnorm(z1min, mean = Ez1, sd = 1) - pnorm(cf, mean = Ez1, sd = 1)
      
    }else{
      if(z1min<cf & z1max >ce){
        prob_n2 <- 0
      }else{
      prob_n2 <- pnorm(z1min, mean = Ez1, sd = 1) + pnorm(z1max, mean = Ez1, sd = 1, lower.tail = FALSE) - prob_n1
    }
 }
} 
  # probability that final sample size is n2max
  if (is.null(theta_tilde)) {
    z1turn <- (w2*qnorm(targetcp) + bdry)/(w1 + w2*sqrt((n2max - n1)/n1))
  } else {
    z1turn <- (w2*(qnorm(targetcp) - theta_tilde/(2*sigma)*sqrt(n2max - n1)) + bdry)/w1
  }
  z1turn <- min(z1turn, z1max)
  z1turn <- max(z1turn, z1min)
  
  if(z1min<cf){
    prob_n2max <- pnorm(z1turn, mean = Ez1, sd = 1) - pnorm(cf, mean = Ez1, sd = 1) 
  }else{
    prob_n2max <- pnorm(z1turn, mean = Ez1, sd = 1) - pnorm(z1min, mean = Ez1, sd = 1)
  }
  

  # average sample size in sloped region
  integrand <- function(z1) {
    n2x_slope(z1 = z1, z1min = z1min, z1max = z1max, n2 = n2, n2max = n2max,
              theta_tilde = theta_tilde, sigma = sigma, targetcp = targetcp,
              bdry = bdry, w1 = w1, w2 = w2, ce=ce, cf=cf)*dnorm(z1 - Ez1)
  }
  avgnslope <- integrate(integrand, lower = z1turn, upper = z1max)$value
  n1*prob_n1 + n2*prob_n2 + n2max*prob_n2max + avgnslope
}


# Expected sample size of CPZ design using sloped SSR rule, when maximum sample
#   size is never achieved (is no "flat" part of SSR rule)
avgn_ssrslope_noflat <- function(n1, n2, n2min, n2max, bdry, w1, w2, ce, cf, 
                          theta, sigma,
                          theta_tilde = NULL, targetcp,
                          z1min, z1max) {
  # expected values of z1
  Ez1 <- expected_z1(theta = theta, sigma = sigma, n1 = n1)
  if (ce < z1max) {
    stop("Early efficacy boundary less than upper promising zone boundary")
  }else{
    if(cf > z1min){
      stop("Early futility boundary greater than lower promising zone boundary") 
    }
    
  }

  # probability that final sample size is n1 (early stop for efficacy)
  prob_n1 <- pnorm(Ez1 - ce)+(1-pnorm(Ez1 - cf))
  # probability that final sample size is n2
  prob_n2 <- pnorm(z1min, mean = Ez1, sd = 1) + pnorm(z1max, mean = Ez1, sd = 1, lower.tail = FALSE) - prob_n1

  # average sample size in sloped region
  integrand <- function(z1) {
    n2x_slope(z1 = z1, z1min = z1min, z1max = z1max, n2 = n2, n2max = n2max,
              theta_tilde = theta_tilde, sigma = sigma, targetcp = targetcp,
              bdry = bdry, w1 = w1, w2 = w2, ce=ce, cf=cf)*dnorm(z1 - Ez1)
  }
  avgnslope <- integrate(integrand, lower = z1min, upper = z1max)$value
  n1*prob_n1 + n2*prob_n2 + avgnslope
}


# Find gamma for JT design with specified average sample size
find_gam <- function(targetn, gam_min, gam_max, nsteps = 10,
                     n1, n2, n2min, n2max, bdry, w1, w2, ce=ce, cf=cf, theta, sigma,
                     theta_tilde,
                     revenue, cost, objective, stepsize,
                     search_z1_min_jt, search_z1_max_jt,
                     search_z1_increment_jt, mc = 1e5,
                     revenue_bound = NULL, theta_revenue_bound = NULL) {
  # local copy of avgn_jt2
  avgn_jt2_local <- function(gam) {
    avgn_jt2(
      n1 = n1, n2 = n2, n2min = n2min, n2max = n2max, bdry = bdry,
      w1 = w1, w2 = w2, ce=ce, cf=cf, theta = theta, sigma = sigma,
      theta_tilde = theta_tilde, gam = gam,
      revenue = revenue, cost = cost, objective = objective, stepsize = stepsize,
      search_z1_min_jt = search_z1_min_jt, search_z1_max_jt = search_z1_max_jt,
      search_z1_increment_jt = search_z1_increment_jt, revenue_bound = revenue_bound,
      theta_revenue_bound = theta_revenue_bound)
  }
  # check gamma range is sufficient
  if(avgn_jt2_local(gam_min) < targetn) {
    stop("gamma is out of range: use smaller gam_min value")
  }
  if (avgn_jt2_local(gam_max) > targetn) {
    stop("gamma is out of range: use larger gam_max value")
  }

  # search for gamma using bisection
  counter <- 0
  while (counter < nsteps) {
    gam_curr <- mean(c(gam_min, gam_max))
    if (avgn_jt2_local(gam_curr) > targetn) {
      gam_min <- gam_curr
    } else {
      gam_max <- gam_curr
    }
    counter <- counter + 1
  }

  mean(c(gam_min, gam_max))
}



# unconditional power JT method at given true theta, using monte carlo integration
uncond_pwr_jt2 <- function(n1, n2, n2min, n2max, bdry, w1, w2, ce=ce, cf=cf,  theta, sigma,
                          theta_tilde, gam,
                          revenue, cost, objective, stepsize = 0.1,
                          search_z1_min_jt, search_z1_max_jt,
                          search_z1_increment_jt, mc = 1e5,
                          revenue_bound = NULL, theta_revenue_bound = NULL) {
  z1min <- find_z1min_jt(n1 = n1, n2 = n2, n2min = n2min, n2max = n2max, bdry = bdry,
                         w1 = w1, w2 = w2, theta_tilde = theta_tilde,
                         sigma = sigma, ce=ce, cf=cf, gam = gam, revenue = revenue, cost = cost,
                         objective = objective, stepsize = stepsize,
                         search_z1_min_jt = search_z1_min_jt,
                         search_z1_max_jt = search_z1_max_jt,
                         search_z1_increment_jt = search_z1_increment_jt,
                         revenue_bound = revenue_bound, theta_revenue_bound = theta_revenue_bound)
  z1max <- find_z1max_jt(n1 = n1, n2 = n2, n2min = n2min, n2max = n2max, bdry = bdry,
                         w1 = w1, w2 = w2, theta_tilde = theta_tilde,
                         sigma = sigma, ce=ce, cf=cf,  gam = gam, revenue = revenue, cost = cost,
                         objective = objective, stepsize = stepsize,
                         search_z1_min_jt = search_z1_min_jt,
                         search_z1_max_jt = search_z1_max_jt,
                         search_z1_increment_jt = search_z1_increment_jt,
                         revenue_bound = revenue_bound, theta_revenue_bound = theta_revenue_bound)
  # set up lookup table for n2x to speed up computation
  z1grid <- seq(z1min, z1max, lookup_table_increments)
  n2x_lookup_table <- n2x_jt(n1 = n1, n2 = n2, n2min = n2min, n2max = n2max, z1 = z1grid,
                             bdry = bdry, w1 = w1, w2 = w2,
                             theta_tilde = theta_tilde, sigma = sigma, ce=ce, cf=cf, 
                             gam = gam, revenue = revenue, cost = cost,
                             objective = objective, stepsize = stepsize,
                             revenue_bound = revenue_bound, theta_revenue_bound = theta_revenue_bound)
  # expected values of z1
  Ez1 <- expected_z1(theta = theta, sigma = sigma, n1 = n1)
  # density of z1
  fz1 <- function(z1) {
    dnorm(z1 - Ez1)
  }

  # faster computation of n2 star
  n2x_jt_fast <- function(z1) {
    out <- numeric(length(z1))
    ssr_region <- (z1 > z1min & z1 < z1max & z1 < ce & z1 > cf)
    
    early_stop_region <- (z1 > ce | z1 < cf)
    
    out[early_stop_region] <- n1
    
    #out[!ssr_region & z1 > ce & z1 < cf] <- n1
    out[!ssr_region & !early_stop_region] <- n2
    
    out[ssr_region] <- n2x_lookup_table[closest_z1(z1[ssr_region], z1grid,
                                                   lookup_table_increments)]
    out
  }

  # fast computation of conditional power using lookup table
  cp_fast <- function(z1) {
    out <- numeric(length(z1))
    z1nostop <- (z1 < ce & z1 > cf)
    out[z1 > ce] <- 1 # cp=1 if early efficacy stop
    out[z1 < cf] <- 0 # cp=0 if early futility stop
    out[z1nostop] <- mapply(FUN = cp,
                           n2x = n2x_jt_fast(z1[z1nostop]),
                           z1 = z1[z1nostop],
                           n1 = n1, bdry = bdry, w1 = w1, w2 = w2, theta_cp = theta,
                           sigma = sigma, ce = ce, cf=cf)
    out
  }

  mean(cp_fast(rnorm(mc, mean = Ez1, sd = 1)))
}



# Unconditonal power for CPZ design with sloped SSR rule
uncond_pwr_ssrslope <- function(n1, n2, n2min, n2max, bdry, w1, w2, ce, cf,
                                theta, sigma,
                                theta_tilde = NULL, targetcp,
                                cpmin, cpmax) {
  # true unconditional power for single true theta
  ucp <- function(theta0) {
    # SSR region boundaries on z scale
    z1min <- cp_to_z1(n1 = n1, n2x = n2max, cp = cpmin, bdry = bdry,
                      w1 = w1, w2 = w2, theta_cp = theta_tilde, sigma = sigma)
    z1max <- cp_to_z1(n1 = n1, n2x = n2, cp = cpmax, bdry = bdry,
                      w1 = w1, w2 = w2, theta_cp = theta_tilde, sigma = sigma)
    # expected values of z1
    Ez1 <- expected_z1(theta = theta0, sigma = sigma, n1 = n1)
    # density of z1
    fz1 <- function(z1) {
      dnorm(z1 - Ez1)
    }
    # probability of early stopping (where cp = 1)
    prob_stop <- pnorm(Ez1 - ce)
    
    # set power as sum of integrals involving conditional power
    integrand_n2 <- function(z1) {
      fz1(z1)*vapply(X = z1, FUN = cp, FUN.VALUE = numeric(1),
                     n2x = n2,
                     n1 = n1, bdry = bdry, w1 = w1, w2 = w2, theta_cp = theta0, sigma = sigma,
                     ce = ce, cf= cf)
    }
    integrand_n2max <- function(z1) {
      fz1(z1)*vapply(X = z1, FUN = cp, FUN.VALUE = numeric(1),
                     n2x = n2max,
                     n1 = n1, bdry = bdry, w1 = w1, w2 = w2, theta_cp = theta0, sigma = sigma,
                     ce = ce, cf= cf)
    }
    integrand_n2slope <- function(z1) {
      fz1(z1)*mapply(FUN = cp,
                     n2x = n2x_slope(z1, z1min = z1min, z1max = z1max, n2 = n2,
                                     n2max = n2max,
                                     theta_tilde = theta_tilde, sigma = sigma,
                                     targetcp = targetcp, bdry = bdry, w1 = w1, w2 = w2, 
                                     ce=ce, cf=cf),
                     z1 = z1,
                     n1 = n1, bdry = bdry, w1 = w1, w2 = w2, theta_cp = theta0, sigma = sigma,
                     ce = ce, cf= cf)
    }
    # turning point where flat increase becomes sloped
    if (is.null(theta_tilde)) {
      z1turn <- (w2*qnorm(targetcp) + bdry)/(w1 + w2*sqrt((n2max - n1)/n1))
    } else {
      z1turn <- (w2*(qnorm(targetcp) - theta_tilde/(2*sigma)*sqrt(n2max - n1)) + bdry)/w1
    }
    z1turn <- min(z1turn, z1max)
    z1turn <- max(z1turn, z1min)

    if(z1min<cf & z1max<ce){
      prob_stop + integrate(f = integrand_n2, lower = cf, upper = z1min)$value +
        integrate(f = integrand_n2max, lower = cf, upper = z1turn)$value +
        integrate(f = integrand_n2slope, lower = z1turn, upper = z1max)$value +
        integrate(f = integrand_n2, lower = z1max, upper = ce)$value
    }else{
      if(z1min>cf & z1max>ce){
        prob_stop + integrate(f = integrand_n2, lower = cf, upper = z1min)$value +
          integrate(f = integrand_n2max, lower = z1min, upper = z1turn)$value +
          integrate(f = integrand_n2slope, lower = z1turn, upper = ce)$value +
          integrate(f = integrand_n2, lower = z1max, upper = ce)$value
      }else{
        if(z1min<cf & z1max>ce){
          prob_stop + integrate(f = integrand_n2, lower = cf, upper = z1min)$value +
            integrate(f = integrand_n2max, lower = cf, upper = z1turn)$value +
            integrate(f = integrand_n2slope, lower = z1turn, upper = ce)$value +
            integrate(f = integrand_n2, lower = z1max, upper = ce)$value
        }else{
          prob_stop + integrate(f = integrand_n2, lower = cf, upper = z1min)$value +
            integrate(f = integrand_n2max, lower = z1min, upper = z1turn)$value +
            integrate(f = integrand_n2slope, lower = z1turn, upper = z1max)$value +
            integrate(f = integrand_n2, lower = z1max, upper = ce)$value
        }
          
      }
    }
    
    
    
  }
  vapply(X = theta, FUN = ucp, FUN.VALUE = numeric(1))
}


# Unconditonal power for CPZ design with sloped SSR rule and with no "flat part"
#  i.e., maximum sample size is never reached
uncond_pwr_ssrslope_noflat <- function(n1, n2, n2min, n2max, bdry, w1, w2, ce, cf,
                                theta, sigma,
                                theta_tilde = NULL, targetcp,
                                z1min, z1max) {
  # true unconditional power for single true theta
  ucp <- function(theta0) {
    # expected values of z1
    Ez1 <- expected_z1(theta = theta0, sigma = sigma, n1 = n1)
    # density of z1
    fz1 <- function(z1) {
      dnorm(z1 - Ez1)
    }
    
    # probability of early stopping (where cp = 1)
    prob_stop <- pnorm(Ez1 - ce) + (1-pnorm(Ez1 - cf))
    
    # set power as sum of integrals involving conditional power
    integrand_n2 <- function(z1) {
      fz1(z1)*vapply(X = z1, FUN = cp, FUN.VALUE = numeric(1),
                     n2x = n2,
                     n1 = n1, bdry = bdry, w1 = w1, w2 = w2, 
                     theta_cp = theta0, sigma = sigma, ce=ce, cf=cf)
    }
    integrand_n2slope <- function(z1) {
      fz1(z1)*mapply(FUN = cp,
                     n2x = n2x_slope(z1, z1min = z1min, z1max = z1max, n2 = n2,
                                     n2max = n2max,
                                     theta_tilde = theta_tilde, sigma = sigma,
                                     targetcp = targetcp, bdry = bdry,
                                     w1 = w1, w2 = w2, ce=ce, cf=cf),
                     z1 = z1,
                     n1 = n1, bdry = bdry, w1 = w1, w2 = w2,
                     theta_cp = theta0, sigma = sigma, ce=ce, cf=cf)
    }

    prob_stop + integrate(f = integrand_n2, lower = cf, upper = z1min)$value +
      integrate(f = integrand_n2max, lower = z1min, upper = z1turn)$value +
      integrate(f = integrand_n2slope, lower = z1min, upper = z1max)$value +
      integrate(f = integrand_n2, lower = z1max, upper = ce)$value
  }
  vapply(X = theta, FUN = ucp, FUN.VALUE = numeric(1))
}

# Sample size re-assessment rule of Jennison-Turnbull
# Gives final re-estimated sample size given the interim statistic is z1
n2x_jt <- function(n1, n2, n2min, n2max, z1, bdry, w1, w2,
                   theta_tilde, sigma , ce, cf, gam, revenue, cost, objective, stepsize = 1,
                   revenue_bound = NULL, theta_revenue_bound = NULL) 
  {
  optimize_objective(n1 = n1, n2 = n2, n2min = n2, n2max = n2max, z1 = z1, bdry = bdry,
                     w1 = w1, w2 = w2, theta_tilde = theta_tilde, sigma = sigma ,
                     ce=ce, cf=cf, gam = gam,
                     revenue = revenue, cost = cost, objective, stepsize = stepsize,
                     revenue_bound = revenue_bound, theta_revenue_bound = theta_revenue_bound)
}


# Turning point in z1 statistic for SSR design, where flat increase to n2max turns to
#  a sloped increase to a sample size smaller than n2max
z1_turn <- function(w1, w2, targetcp, theta_tilde, sigma, n1, n2max, bdry) {
  (w2*(qnorm(targetcp) - theta_tilde/(2*sigma)*sqrt(n2max - n1)) + bdry)/w1
}

# Final re-estimated sample size using sloped SSR rule targeting a fixed CP at theta_tilde
n2x_slope <- function(z1, z1min, z1max, n2, n2max,
                      theta_tilde, sigma, targetcp, bdry, w1, w2, ce, cf) {
  out <- numeric(length(z1))
  
  ssr_region <- (z1 >cf & z1 > z1min & z1 < z1max & z1 < ce)
  early_stop_region <- (z1 > ce | z1 < cf)
  
  out[!ssr_region & early_stop_region] <- n1
  out[!ssr_region & !early_stop_region] <- n2
  
  if (is.null(theta_tilde)) {
    z1turn <- (w2*qnorm(targetcp) + bdry)/(w1 + w2*sqrt((n2max - n1)/n1))
    flat_region <- ssr_region & z1 < z1turn
    slope_region <- ssr_region & z1 >= z1turn
    out[flat_region] <- n2max
    out[slope_region] <- n1 + ((w2*qnorm(targetcp) + bdry - z1[slope_region]*w1))^2*n1/(w2^2*z1[slope_region]^2)
  } else {
    z1turn <- (w2*(qnorm(targetcp) - theta_tilde/(2*sigma)*sqrt(n2max - n1)) + bdry)/w1
    flat_region <- ssr_region & z1 < z1turn
    slope_region <- ssr_region & z1 >= z1turn
    out[flat_region] <- n2max
    out[slope_region] <- n1 + ((2*sigma)/theta_tilde*(qnorm(targetcp) + (bdry - w1*z1[slope_region])/w2))^2
  }
  out
}


# Find lower limit of "promising zone" in JT method; z1 value where sample size
#  starts increasing
find_z1min_jt <- function(n1, n2, n2min, n2max, bdry, w1, w2, theta_tilde,
                          sigma , ce, cf, gam, revenue, cost, objective, stepsize = 1,
                          search_z1_min_jt = 0, search_z1_max_jt = 4.0,
                          search_z1_increment_jt = search_z1_increment_jt,
                          revenue_bound = NULL, theta_revenue_bound = NULL) {
  z1current <- search_z1_min_jt
  while (n2x_jt(n1 = n1, n2 = n2, n2min = n2min, n2max = n2max, z1 = z1current,
                bdry = bdry, w1 = w1, w2 = w2, theta_tilde = theta_tilde,
                sigma = sigma , ce=ce, cf=cf, gam = gam, revenue = revenue, cost = cost,
                objective = objective, stepsize = stepsize, revenue_bound = revenue_bound,
                theta_revenue_bound = theta_revenue_bound) <= n2 &&
         z1current < search_z1_max_jt) {
    z1current <- z1current + search_z1_increment_jt
  }
  if (z1current >= search_z1_max_jt) {
    stop("did not locate z1min: try increasing search_z1_max_jt")
  }
  if (z1current <= search_z1_min_jt) {
    stop("did not locate z1min: decrease search_z1_min_jt")
  }
  z1current - search_z1_increment_jt
}



# Find upper limit of "promising zone" in JT method; the z1 value where
#  sample size stops increasing
find_z1max_jt <- function(n1, n2, n2min, n2max, bdry, w1, w2, theta_tilde,
                          sigma , ce, cf, gam, revenue, cost, objective, stepsize = 1,
                          search_z1_min_jt = 0, search_z1_max_jt = 4.0,
                          search_z1_increment_jt = search_z1_increment_jt,
                          revenue_bound = NULL, theta_revenue_bound = NULL) {
  z1current <- search_z1_max_jt
  while (n2x_jt(n1 = n1, n2 = n2, n2min = n2min, n2max = n2max, z1 = z1current,
                bdry = bdry, w1 = w1, w2 = w2, theta_tilde = theta_tilde,
                sigma = sigma, ce=ce, cf=cf, gam = gam, revenue = revenue, cost = cost,
                objective = objective, stepsize = stepsize, revenue_bound = revenue_bound,
                theta_revenue_bound = theta_revenue_bound) <= n2 &&
         z1current > search_z1_min_jt) {
    z1current <- z1current - search_z1_increment_jt
  }
  if (z1current >= search_z1_max_jt) {
    stop("did not locate z1max: try increasing search_z1_max_jt")
  }
  if (z1current <= search_z1_min_jt) {
    stop("did not locate z1max: decrease search_z1_min_jt")
  }
  z1current + search_z1_increment_jt
}

# Find closest grid points to z1
closest_z1 <- function(z1, z1grid, lookup_table_increments) {
  grid_midpts <- c(-Inf, z1grid[-length(z1grid)] + 0.5*lookup_table_increments)
  vapply(X = z1, FUN = function(x) {sum(x > grid_midpts)}, FUN.VALUE = numeric(1))
}

# Faster computation of n2 star using a lookup table
n2x_jt_fast <- function(z1, z1min, z1max, ce, cf, z1grid, n2x_lookup_table) {
  out <- numeric(length(z1))
  ssr_region <- (z1 >cf & z1 > z1min & z1 < z1max & z1 < ce)
  early_stop_region<-(z1> ce | z1< cf)
  
  out[early_stop_region] <- n1
  out[!ssr_region & !early_stop_region] <- n2
  
  out[ssr_region] <- n2x_lookup_table[closest_z1(z1[ssr_region], z1grid,
                                                 lookup_table_increments)]
  out
}


# Final sample size optimizing the objective function in the sense of Jennison-Turnbull,
#  allow specification of lower bound on revenue function (e.g., CP) for constrained optimization
optimize_objective <- function(n1, n2, n2min, n2max, z1, bdry, w1, w2,
                               theta_tilde, sigma , ce, cf, gam, 
                               revenue, cost, objective,
                               stepsize = 1,
                               revenue_bound = NULL, theta_revenue_bound = NULL) {
  x <- seq(n2min, n2max, stepsize)
  # z1 might be a vector, so need function that optimizes at each z1 value
  optimize_at_each_z1 <- function(z1_singleton) {
    if(z1_singleton > ce | z1_singleton < cf ){
      
      n1
    }else{
      
      revenue_bound_satisfied <- rep(TRUE, length(x))
      if (!is.null(revenue_bound)) {
        if (!is.null(theta_revenue_bound)) {
          r <- vapply(X = x, FUN = revenue, FUN.VALUE = numeric(1),
                      n1 = n1, n2 = n2, z1 = z1_singleton, bdry = bdry, w1 = w1, w2 = w2,
                      theta_tilde = theta_revenue_bound, sigma = sigma,
                      ce = ce, cf= cf)
        } else {
          r <- vapply(X = x, FUN = revenue, FUN.VALUE = numeric(1),
                      n1 = n1, n2 = n2, z1 = z1_singleton, bdry = bdry, w1 = w1, w2 = w2,
                      theta_tilde = theta_tilde, sigma = sigma,
                      ce = ce, cf= cf)
        }
        revenue_bound_satisfied <- r > revenue_bound
      }
      if (any(revenue_bound_satisfied)) {
        x <- x[revenue_bound_satisfied]
        y <- vapply(X = x, FUN = objective, FUN.VALUE = numeric(1),
                    n1 = n1, n2 = n2, z1 = z1_singleton, bdry = bdry, w1 = w1, w2 = w2,
                    theta_tilde = theta_tilde, sigma = sigma, ce=ce, cf=cf,
                    gam = gam,
                    revenue = revenue, cost = cost)
        
        x[which.max(y)]
      } else {
        # minimum revenue requirements not met, return smallest possible sample size
        n2min
      }
      
    }
    
  }
  # vector of optimized values across all z1 inputs
  vapply(X = z1, FUN = optimize_at_each_z1, FUN.VALUE = numeric(1))
}



# objective function defined abstractly as revenue minus cost
objective <- function(n2x, n1, n2, z1, bdry, w1, w2, theta_tilde, sigma, 
                      ce, cf, gam, revenue, cost) {
  revenue(n1 = n1, n2 = n2, n2x = n2x, z1 = z1, bdry = bdry,
          w1 = w1, w2 = w2, theta_tilde = theta_tilde, sigma,
          ce = ce, cf= cf) -
    cost(n1 = n1, n2 = n2, n2x = n2x, z1 = z1, bdry = bdry,
         theta_tilde = theta_tilde, sigma = sigma, gam = gam)
}

# "revenue" defined in terms of conditional power
revenue <- function(n1, n2, n2x, z1, bdry, w1, w2, theta_tilde, sigma,
                    ce, cf) {
  cp(n1 = n1, n2x = n2x, z1 = z1, bdry = bdry,
                 w1 = w1, w2 = w2, theta_cp = theta_tilde, sigma = sigma,
                 ce = ce, cf= cf) 
}

# "cost" is linearly related to sample size
cost <- function(n1, n2, n2x, z1, bdry, w1, w2, theta_tilde, sigma, gam) {
  gam*(n2x - n2)
}


## find_gam_match_pow finds gamma to match desired power for the JT design
## at desired target delta / smallest clinically meaningful delta

find_gam_match_pow<- function(target.pow, gam_min, gam_max, nsteps,
                              n1, n2, n2min, n2max, bdry, w1, w2, ce, cf, theta, sigma,
                              theta_tilde,
                              revenue, cost, objective, stepsize,
                              search_z1_min_jt, search_z1_max_jt,
                              search_z1_increment_jt, mc = 1e5,
                              revenue_bound = NULL, theta_revenue_bound = NULL) {
  
  uncond_pwr_jt2_local<-function(gam)
  {
    uncond_pwr_jt2 (n1, n2, n2min, n2max, bdry, w1, w2, ce=ce, cf=cf, theta, sigma,
                    theta_tilde, gam, revenue, cost, objective,
                    stepsize = 0.1, search_z1_min_jt=-4, search_z1_max_jt,
                    search_z1_increment_jt, mc = 1e5,
                    revenue_bound = NULL, theta_revenue_bound = NULL)
    
  }
  
  # check gamma range is sufficient
  if(uncond_pwr_jt2_local(gam_min) < target.pow ) {
    stop("gamma is out of range: use smaller gam_min value")
  }
  if (uncond_pwr_jt2_local(gam_max) > target.pow) {
    stop("gamma is out of range: use larger gam_max value")
  }
  
  # search for gamma using bisection
  counter <- 0
  while (counter < nsteps) {
    gam_curr <- mean(c(gam_min, gam_max))
    if (uncond_pwr_jt2_local(gam_curr) >target.pow ) {
      gam_min <- gam_curr
    } else {
      gam_max <- gam_curr
    }
    counter <- counter + 1
  }
  
  mean(c(gam_min, gam_max))
}



