# .......................................................................
# REQUIRED PACKAGES
# .......................................................................
# install.packages(c("here", "nimble", "coda", "stringr", "sf"))
library(here)
library(nimble)
library(coda)
library(stringr)
library(sf)

setwd(here::here("data"))
load("isdm_data_constants.RData")
cov <- st_read("covariates.shp")

# inspect data and constants used in modeling
# refer to README file for description of what they are
str(data)
str(constants)

# inspect environmental predictors
plot(cov[1], border = NA)
plot(cov[2], border = NA)

# .......................................................................
# MODEL CODE
# .......................................................................

harvest <- nimble::nimbleCode( {
  
  # .............................................................
  # PRIORS
  # .............................................................
  
  ## Priors for SDM ##
  
  # CAR prior for spatial random effect
  s[1:ncell] ~ dcar_normal(adj[1:neigh], weights[1:neigh], num[1:ncell], tau)
  # precision of CAR prior
  tau ~ dgamma(1, 1)
  
  # regression coefficient for canopy cover
  b_forest ~ dnorm(0, 2)
  # regression coefficient for impervious cover
  b_imperv ~ dnorm(0, 2)
  
  # .............................................................
  # LIKELIHOOD
  # .............................................................
  
  # SDM - model for the latent state
  for(i in 1:ncell){
    log(lambda[i]) <- s[i] + b_forest*forest[i] + b_imperv*imperv[i]
    n[i] ~ dpois(lambda[i])
  }
  
  # Harvest submodel
  for(c in 1:ncounty){
    # e: effort coefficient
    # low & high: indices to account for number of cells within each county
    mu_p[c] <- e[c]*sum(lambda[low[c]:high[c]])
    w[c] ~ dpois(mu_p[c])
  }
  
} )

# .......................................................................
# PREPARE MODEL TO RUN
# .......................................................................

# function to provide random initial values for parameters
inits <- function() {
  base::list(z = rep(1, constants$ncell),
             b_forest = runif(1, -1, 1),
             b_imperv = runif(1, -1, 1),
             tau = rgamma(1, 1, 1),
             s = rep(0, base::length(data$num)))}

# parameters to monitor
keepers <- c("lambda", 'b_forest', "b_imperv")

data_harvest_only <- list(
  e = data$e, 
  w = data$w, 
  num = data$num, 
  adj = data$adj, 
  weights = data$weights, 
  forest = data$forest, 
  imperv = data$imperv)

constants_harvest_only <- list(
  ncell = constants$ncell, 
  ncounty = constants$ncounty,
  low = constants$low, 
  high = constants$high, 
  neigh = constants$neigh
)

# Will have to run chains for much longer to approach convergence
# running with 200 iterations took about 10 minutes on my laptop with 4 cores
# to speed things up, particularly for longer chains, you can run chains in parallel
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
nc <- 3 # number of chains
nb <- 100 # number of initial MCMC iterations to discard
ni <- nb + 100 # total number  of iterations

# .......................................................................
# RUN MODEL
# .......................................................................

# create model
model <- nimble::nimbleModel(code = harvest, 
                             data = data_harvest_only, 
                             constants = constants_harvest_only, 
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
                           nchains = nc)

# .......................................................................
# INSPECT RESULTS
# .......................................................................

# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc(lapply(samples, coda::mcmc))

# Look at traceplots of the first five parameters
par(mfrow=c(1,2))
coda::traceplot(samples_mcmc[, 1:2])
# calculate Rhat convergence diagnostic for first five parameters
coda::gelman.diag(samples_mcmc[,1:2])

# extract mean and SD occurrence probability of each grid cell
samplesdf <- do.call(rbind, samples_mcmc)
lambda <- samplesdf[, which(stringr::str_detect(string = colnames(samplesdf), pattern = 'lambda\\['))]
lambda_mean <- apply(lambda, 2, mean)
lambda_sd <- apply(lambda, 2, sd)

# map mean and standard deviation of occurrence probability
par(mfrow=c(1,1))
cov$lambda_mean <- lambda_mean
plot(cov["lambda_mean"], border = NA)
cov$lambda_sd <- lambda_sd
plot(cov["lambda_sd"], border = NA)