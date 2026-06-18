# Simulation code for assessing necessary minimum sample sizes to detect a treatment
#   effect on annual survival of northern bobwhite (Colinus virginianus). From 
#   Lewis, W. B., and J. A. Martin. 2026. Estimating statistical power for detecting 
#   treatment effects on northern bobwhite (Colinus virginianus) vital rates.

# We use data simulations to assess statistical power for detecting a significant treatment
#   effect on productivity (chicks/female) between a control and treatment group. Simulations
#   vary in terms of the sample size (number of females with breeding data ranging from 40 - 400),
#   baseline productivity in the control group (2 or 4 chicks/female), and productivity
#   in the treatment group being 10%, 25%, 50%, 75%, or 100% higher relative to the
#   control group.
# We simulate per-capita productivity via a zero-inflated Poisson model based on 1) the 
#   probability of successfully hatching at least one brood during the breeding season,
#   and 2) the number of hatched chicks per successful nest. First, we calculate the 
#   probability of a female fledging at least one nest based on the experimental group; 
#   then use this probability in a Bernoulli distribution to simulate an indicator 
#   variable P denoting whether or not each bird successfully hatched a brood. We then
#   model the expected number of chicks from successful nests based on the experimental
#   group. We simulate the number of chicks produced by each female from a Poisson process
#   based on the product of P and the expected number of chicks. We specify parameters for
#   each submodel such that variation in productivity is ~ 90% driven by variation in the
#   probability of successfully hatching a brood. Lastly, we analyze simulated datasets as
#   the opposite of the data simulation process.


require(jagsUI)
require(coda)




# Simulation Parameters ---------------------------------------------------------
H.baseline <- c(2, 4)
H.baseline.length <- length(H.baseline)
Nbirds <- seq(40, 400, by=40)
Nbirds.len <- length(Nbirds)
Effect <- c(0.1,0.25,0.5,0.75,1)
Effect.len <- length(Effect)
Nsims <- 1000

# Intercept for calculating probability of hatching at least one brood (theta) by baseline productivity in the control group (2/4 chicks/female).
b.0 <- c(-1.442, -0.560)
# Treatment effect on theta. This varies by effect size (columns) and baseline productivity (rows).
b.1 <- matrix(c(0.107, 0.251, 0.463, 0.648, 0.813,
                0.138, 0.332, 0.637, 0.932, 1.229), nrow = 2, byrow = T)
# Intercept for calculating number of chicks hatched from successful nests (lambda) by baseline productivity in the control group (2/4 chicks/female)
g.0 <- c(2.347, 2.398)
# Treatment effect on lambda. This varies by effect size; however, most variation in productivity is attributed to variation in theta rather than
#   lambda. As such, the treatment effect varies little by baseline productivity.
g.1 <- c(0.01, 0.025, 0.049, 0.072, 0.095)





# Simulating and analyzing datasets --------------------------------------------

for(h in 1:H.baseline.length){
  for(b in 1:Nbirds.len){
    for(e in 1:Effect.len){
      
      sim.data <- mcmc.samps <- GR <- vector(mode = "list", length = Nsims)
      
      # Looping through simulations
      for(x in 1:Nsims){
        
        
        ## Simulating datasets -------------------------------------------------
        
        sim.data.info <- data.frame(Group = rep(NA,times = Nbirds[b]), theta = NA, 
                                    P = NA, lambda = NA, H = NA)
        
        ## Assigning experimental group (0=Control, 1=Treatment)
        sim.data.info$Group <- rep(0:1, each = Nbirds[b]/2)
        
        ## Simulating productivity for each female
        for(i in 1:nrow(sim.data.info)){
          
          ## Zero-inflation process (P), based on theta
          sim.data.info$theta[i] <- plogis(b.0[h] + b.1[h,e] * sim.data.info$Group[i])
          sim.data.info$P[i] <- rbinom(1, 1, sim.data.info$theta[i])
          
          ## Conditional productivity, based on P and lambda
          sim.data.info$lambda[i] <- exp(g.0[h] + g.1[e] * sim.data.info$Group[i])
          sim.data.info$H[i] <- rpois(1, sim.data.info$P[i] * sim.data.info$lambda[i])
          
        } # end i
        
        
        
        ## Analyzing simulated productivity datasets in JAGS -------------------
        
        sink("productivitymodel.jags")
        cat("
            model{
              
            # Priors
            b.0 ~ dnorm(0, 0.001)
            b.1 ~ dnorm(0, 0.001)
            g.0 ~ dnorm(0, 0.001)
            g.1 ~ dnorm(0, 0.001)
        
            # Model
            for(f in 1:N){
              logit(theta[f]) <- b.0 + b.1 * Group[f]
              p[f] ~ dbern(theta[f])
              lambda[f] <- exp(g.0 + g.1 * Group[f])
              H[f] ~ dpois(p[f] * lambda[f] + 0.000001) # JAGS can have problems with 0 expected value in Poisson
            }
        
            # Derived Parameters. Calculating expected productivity in each experimental group, then calculating
            #   the absolute (H.Effect) and relative (Relative.H.Effect) difference in productivity between
            #   the control and treatment groups.
            for(g in 1:Ngroups){
              logit(p.pred[g]) <- b.0 + b.1 * Group.pred[g]
              lambda.pred[g] <- exp(g.0 + g.1 * Group.pred[g])
              H.pred[g] <- p.pred[g] * lambda.pred[g]
            }
            H.Effect <- H.pred[2] - H.pred[1]
            Relative.H.Effect <- H.Effect/(H.pred[1])
             
            }
          ",fill=TRUE)
        sink()
        
        params <- c("b.0","b.1","g.0","g.1","H.Effect","Relative.H.Effect")
        
        mod.data <- list(N = Nbirds[b],
                         Group = sim.data.info$Group,
                         H = sim.data.info$H,
                         Ngroups = length(unique(sim.data.info$Group)),
                         Group.pred = unique(sim.data.info$Group))
        
        inits.fun <- function() list(p = ifelse(mod.data$H>0, 1, 0),
                                     b.0 = runif(1, -1, 1),
                                     g.0 = runif(1, 2, 3),
                                     b.1 = runif(1, -0.5, 0.5),
                                     g.1 = runif(1, -0.5, 0.5))
        
        jags.post <- jags.basic(data = mod.data, 
                                model.file = "productivitymodel.jags",
                                parameters.to.save = params,
                                n.chains = 3, n.adapt = 3000, n.burnin = 3000,
                                n.iter = 30000, parallel = TRUE, n.thin = 10)
        
        sim.data[[x]] <- sim.data.info
        mcmc.samps[[x]] <- rbind(jags.post[[1]], jags.post[[2]], jags.post[[3]])
        GR[[x]] <- gelman.diag(jags.post, multivariate=F)$psrf[,2]
        
        rm(list = c("sim.data.info", "params", "mod.data", "inits.fun", "jags.post"))
        
      } # end x
      
      sim.params <- list(H.baseline = H.baseline[h],
                         Nbirds = Nbirds[b],
                         Effect = Effect[e],
                         b.0 = b.0[h],
                         b.1 = b.1[h, e],
                         g.0 = g.0[h],
                         g.1 = g.1[e],
                         Nsims = Nsims)
      
      sim.out <- list(sim.params = sim.params,
                      sim.data = sim.data,
                      mcmc.samps = mcmc.samps,
                      GR = GR)
      
      filename <- paste0("Bobwhite_productivity_simulationoutput_", H.baseline[h], "BaselineProductivity_", Nbirds[b], "Birds_", Effect[e], "TreatmentEffect.gzip")
      save(sim.out, file=filename)
      
      rm(list = c("sim.params", "sim.data", "mcmc.samps", "GR", "sim.out", "filename"))
      
    } # end e
  } # end b
} # end h

