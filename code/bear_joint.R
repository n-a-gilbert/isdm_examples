library(here)
library(nimble)
library(coda)

setwd(here::here("data"))

load("isdm_data_constants.RData")

str(data)
str(constants)

inits <- function() {
  base::list(z = rep(1, constants$ncell),
             alpha = runif(1,-1,1), 
             b_forest = runif(1, -1, 1),
             b_imperv = runif(1, -1, 1),
             tau = rgamma(1, 1, 1), 
             s = rep(0, base::length(data$num)),
             p = runif(1, 0, 1),
             eps_psi = rnorm(constants$ncell, 0, runif(constants$ncell, 0, 1)),
             sd_psi = runif(1, 0, 2))
}

joint <- nimble::nimbleCode( {
  
  # CAR prior for spatial random effects
  s[1:ncell] ~ dcar_normal(adj[1:k],
                           weights[1:k],
                           num[1:ncell],
                           tau)
  
  tau ~ dgamma(1, 1)
  
  # Other priors
  alpha ~ dlogis(0, 1)
  b_forest ~ dlogis(1, 1)
  b_imperv ~ dlogis(-1, 1)
  p ~ dunif(0, 1)
  
  for(i in 1:ncell){
    eps_psi[i] ~ dnorm(0, sd = sd_psi)
  }
  
  sd_psi ~ dunif(0, 2)
  
  # prediction-resolution grid
  for(i in 1:ncell){
    log(lambda[i]) <- min(s[i]
                          + alpha
                          + b_forest*forest[i]
                          + b_imperv*imperv[i]
                          + eps_psi[i], 10)
    psi[i] <- 1 - exp(-lambda[i])
    z[i] ~ dbern(psi[i])
  }
  
  # cameras within prediction-resolution grid
  for(j in 1:nsite){
    muy[j] <- z[cell[j]]*p
    y[j] ~ dbin(muy[j], nsurveys[j])
  }
  
  for(i in 1:ncounty){
    # e: effort coefficient
    # low & high: indices to account for number of cells within each county
    mu_p[i] <- e[i]*inprod(lambda[low[i]:high[i]], z[low[i]:high[i]])
    w[i] ~ dpois(mu_p[i])
  }
  
} )

# parameters to monitor
keepers <- c('psi', 'p', "alpha", 'b_forest', "b_imperv")

nc <- 3 # number of chains
nb <- 1000 # number of initial MCMC iterations to discard
ni <- nb + 1000 # total number  of iterations
nt <- 5 # thinning interval

model <- nimble::nimbleModel(code = joint, 
                             data = data, 
                             constants = constants, 
                             inits = inits())

# check to see if everything is initialized
model$initializeInfo()

# compile the model
c_model <- nimble::compileNimble(model)

model_conf <- nimble::configureMCMC(model)

model_conf$addMonitors(keepers)

model_mcmc <- nimble::buildMCMC(model_conf)

c_model_mcmc <- nimble::compileNimble(model_mcmc, project = model)

samples <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc,
                           thin = nt)

samples_mcmc <- as.mcmc(lapply(samples, mcmc))

coda::traceplot(samples_mcmc[, "b_forest"])
coda::gelman.diag(samples_mcmc[,1:5])
                  