# Simulation code for assessing necessary minimum sample sizes to detect a treatment
#   effect on annual survival of northern bobwhite (Colinus virginianus). From 
#   Lewis, W. B., and J. A. Martin. 2026. Estimating necessary minimum sample sizes
#   for detecting treatment effects on annual survival and per-capita productivity 
#   of northern bobwhite (Colinus virginianus).

# We use data simulations to assess statistical power for detecting a significant treatment
#   effect on survival between a control and treatment group. Simulations vary in terms
#   of the sample size (number of radtiotransmittered birds from 40 - 400),
#   baseline annual survival in the control group (0.12 or 0.22), and survival in 
#   the treatment group being 10%, 25%, 50%, 75%, or 100% higher relative to the
#   control group.
# We simulate weekly capture histories for birds over a 104 week period (52 weeks
#   per year and 2 years). We simulate birds as being captured and fitted with 
#   radiotransmitters during the first week of March and the first week of October.
#   We divide start times of birds evenly between experimental groups (control and
#   treatment), deployment seasons (March and October), and years (first and second).
#   For example, a total sample size of 400 corresponds to 50 transmitters deployed
#   in each experimental group in each deployment season and year. We specify birds
#   as being alive in the tracking period in which they transmitter was deployed;
#   we simulate the alive/dead state (z) of each individual in subsequent weeks 
#   based on a two-step process incorporating survival and radio failure. First, 
#   we model the weekly probability of retaining an active transmitter, i.e., the
#   radiotrasnmitter is still attached to the bird and actively transmitting. Next,
#   we calculate the 26-week seasonal survival for each individual in each week 
#   based on an intercept, experimental group of the bird, and season of the 
#   tracking week. Season is defined as breeding (April - September) and 
#   non-breeding (October - March). We specify non-breeding survival being 
#   ~ 10 points higher on the percentage scale than breeding survival (e.g., 
#   50% non-breeding and 40% breeding survival). We then convert period survival
#   to weekly survival and simulate weekly capture histories based on a Bernoulli
#   proces and the weekly survival probability, conditional on retaining an active
#   transmitter and being alive in the previous period. We assume perfect detection
#   of birds during telemetry as in other studies. We subset capture histories to
#   the period between deployment and either mortality or censoring, the latter of
#   which occurs on radio failure or if a bird is alive at the end of the study.
#   Lastly, we analyze simulated datasets with a Bayesian known-fate model. The 
#   modeling process is the opposite of the data simulation process; however, we
#   do not directly model tag retention because capture histories are constrained
#   to periods with an active transmitter.



require(jagsUI)
require(coda)




# Simulation Parameters ---------------------------------------------------------
SA.baseline <- c(0.12, 0.22)
SA.baseline.length <- length(SA.baseline)
Nbirds <- seq(40, 400, by=40)
Nbirds.len <- length(Nbirds)
Effect <- c(0.1,0.25,0.5,0.75,1)
Effect.len <- length(Effect)
Nyears <- 2
Nsims <- 1000

# Survival intercept of baseline annual survival in the control group (0.12/0.22)
a.0 <- c(-0.405, 0.088)
# Treatment effect on survival. This varies by effect size (columns) and baseline survival (rows).
a.1 <- matrix(c(c(0.074, 0.177, 0.332, 0.470, 0.598),
                c(0.092, 0.223, 0.426, 0.619, 0.807)), byrow = T, nrow = 2)
# Difference in survival between breeding and non-breeding season (intercept corresponds to non-breeding). 
#   This varies by baseline survival.
a.2 <- c(-0.442, -0.403)

# Weekly probability of retaining an active radiotransmitter. This corresponds to a median
#   time to failure of ~ 33 weeks.
delta <- 0.9792

# Number of tracking periods in each season (breeding vs. non-breeding)
Nweeks <- 26
# Total number of tracking periods across the two years
Nperiods <- Nweeks*2*Nyears
# Season of each tracking period, 0 is non-breeding, 1 is breeding. First week corresponds to the start of
#   the non-breeding season (first week of October) in the first year.
season.sim <- rep(rep(0:1, each = Nweeks), times = Nyears)
# Capture/deployment weeks in each year. Rows correspond to years while capture periods correspond to columns (1st: early October, 2nd: early March)
Deploy.sim <- matrix(c(1 + (1:Nyears-1) * 22, Nweeks * 2 + 1 + (1:Nyears-1) * 22), ncol = 2, byrow = T)





# Simulating and analyzing datasets --------------------------------------------

for(sa in 1:SA.baseline.length){
  for(b in 1:Nbirds.len){
    for(e in 1:Effect.len){
      
      sim.data <- CH <- mcmc.samps <- GR <- vector(mode = "list", length = Nsims)
      
      # Looping through simulations
      for(x in 1:Nsims){
        
        
        ## Simulating datasets -------------------------------------------------
        
        sim.data.info <- data.frame(Group = rep(NA, times = Nbirds[b]), Year = NA, Period = NA,
                                    First = NA, Last = NA)
        
        ## Assigning experimental group (0=Control, 1=Treatment), year (1 or 2),
        ##  and deployment period (1 is October, 2 is March).
        sim.data.info$Group <- rep(0:1, each = Nbirds[b]/2)
        sim.data.info$Year <- rep(rep(1:2, each = Nbirds[b]/4), times = 2)
        sim.data.info$Period <- rep(rep(1:2, each = Nbirds[b]/8), times = 4)
        for(k in 1:nrow(sim.data.info)){
          sim.data.info$First[k] <- Deploy.sim[sim.data.info$Year[k], sim.data.info$Period[k]]
        } # end k
        
        ## Creating capture history matrix
        
        CH.sim <-  matrix(NA, nrow=Nbirds[b], ncol=Nperiods)
        
        for(i in 1:nrow(CH.sim)){
          
          CH.sim[i, sim.data.info$First[i]] <- 1 # Set alive on capture day

          # Filling out weekly until bird either is censored or dies. Birds generally
          #   tracked every few days during studies, so assuming that will be tracked
          #   every week
          
          for(t in (sim.data.info$First[i]+1):Nperiods){
            
            # Simulating indicator variable TA denoting if tag still active (TA=1) on a given day
            # If tag fails (TA=0), censoring capture histories
            TA <- rbinom(1, 1, delta)
            if(TA == 0){
              sim.data.info$Last[i] <- t-1
              break
            } # end if
            
            # If transmitter is active, estimating season-specific survival then converting to
            #   weekly survival
            S <- plogis(a.0[sa] + a.1[sa, e] * sim.data.info$Group[i] + a.2[sa] * season.sim[t-1])
            phi <- S^(1/Nweeks)
            CH.sim[i, t] <- rbinom(1, 1, phi * CH.sim[i, t-1])
            
            # Censoring in event of mortality of if still alive at the end of the survey period
            if(CH.sim[i, t]==0){
              sim.data.info$Last[i] <- t
              break
            } # end if
            if(t==Nperiods & CH.sim[i,t]==1){
              sim.data.info$Last[i] <- t
            } # end if
          } # end t
        } # end i
        
        
        # Can occasionally have instances where tag failed immediately after release before
        #   getting any tracking data. These have first and last cap dates the same. Removing.
        CH.sim <- CH.sim[sim.data.info$First!=sim.data.info$Last,]
        sim.data.info <- sim.data.info[sim.data.info$First!=sim.data.info$Last,]
        
        
        
        ## Analyzing simulated survival datasets in JAGS -----------------------
        
        sink("survivalmodel.jags")
        cat("
          model{
              
          # Priors
          a.0 ~ dnorm(0, 0.001)
          a.1 ~ dnorm(0, 0.001)
          a.2 ~ dnorm(0, 0.001)
        
          # Model
          for(i in 1:N){
            for(w in (first[i]+1):last[i]){
              logit(S[i,w]) <- a.0 + a.1 * Group[i] + a.2 * Season[w-1]
              phi[i,w] <- pow(S[i,w], 1/Nweeks)
              CH[i,w] ~ dbern(phi[i,w])
            }
          }
        
          # Derived Parameters. Calculating expected seasonal survival in each experimental group and season,
          #   then calculating expected annual survival in each experimental group. Finally, calculating
          #   the absolute (S.Effect) and relative (Relative.S.Effect) difference in survival between
          #   the control and treatment groups.
          for(g in 1:Ngroups){
            for(s in 1:Nseasons){
              logit(S.pred[g,s]) <- a.0 + a.1 * Group.pred[g] + a.2 * Season.pred[s]
            }
            Annual.S.pred[g] <- prod(S.pred[g,1:Nseasons])
          }
          S.Effect <- Annual.S.pred[2] - Annual.S.pred[1]
          Prop.S.Effect <- S.Effect/(Annual.S.pred[1])
             
            }
          ",fill=TRUE)
        sink()
        
        params <- c("a.0","a.1","a.2","S.Effect","Prop.S.Effect")
        
        mod.data <- list(CH = CH.sim,
                         first = sim.data.info$First,
                         last = sim.data.info$Last,
                         N = nrow(sim.data.info),
                         Season = season.sim,
                         Nseasons = length(unique(season.sim)),
                         Season.pred = unique(season.sim),
                         Group = sim.data.info$Group,
                         Ngroups = length(unique(sim.data.info$Group)),
                         Group.pred = unique(sim.data.info$Group),
                         Nweeks = Nweeks)
        
        inits.fun <- function() list(a.0 = runif(1, qlogis(0.95), qlogis(0.999)),
                                     a.1 = runif(1, -0.05, 0.05),
                                     a.2 = runif(1, -0.05, 0.05))
        
        jags.post <- jags.basic(data = mod.data, 
                                model.file = "survivalmodel.jags",
                                parameters.to.save = params,
                                n.chains = 3, n.adapt = 3000, n.burnin = 3000,
                                n.iter = 30000, parallel = TRUE, n.thin = 10)
        
        sim.data[[x]] <- sim.data.info
        CH[[x]] <- CH.sim
        mcmc.samps[[x]] <- rbind(jags.post[[1]], jags.post[[2]], jags.post[[3]])
        GR[[x]] <- gelman.diag(jags.post, multivariate=F)$psrf[,2]
        
        rm(list = c("sim.data.info", "CH.sim", "TA", "S", "phi", "params", "mod.data", "inits.fun", "jags.post"))
        
      } # end x
      
      sim.params <- list(SA.baseline = SA.baseline[sa],
                         Nbirds = Nbirds[b],
                         Effect = Effect[e],
                         a.0 = a.0[sa],
                         a.1 = a.1[sa, e],
                         a.2 = a.2[sa],
                         Nsims = Nsims,
                         delta = delta,
                         Nweeks = Nweeks,
                         season = season.sim)
      
      sim.out <- list(sim.params = sim.params,
                      sim.data = sim.data,
                      CH = CH,
                      mcmc.samps = mcmc.samps,
                      GR = GR)
      
      filename <- paste0("Bobwhite_annualsurvival_simulationoutput_", SA.baseline[sa], "BaselineAnnualSurvival_", Nbirds[b], "Birds_", Effect[e], "TreatmentEffect.gzip")
      save(sim.out, file=filename)
      
      rm(list = c("sim.params", "sim.data", "CH", "mcmc.samps", "GR", "sim.out", "filename"))
      
      
    } # end e
  } # end b
} # end h
