# .......................................................................
# REQUIRED PACKAGES
# .......................................................................
# install.packages(c("here", "nimble", "coda", "pROC", "stringr", "sf"))
library(here)
library(nimble)
library(coda)
library(pROC)
library(stringr)
library(sf)

setwd(here::here("data"))
load("isdm_data_constants.RData")
cov <- st_read("covariates.shp")

# inspect data and constants used in modeling
# refer to README file for description of what they are
str(data)
str(constants)

# .......................................................................
# MODEL CODE
# .......................................................................

joint <- nimble::nimbleCode( {
  
  # .............................................................
  # PRIORS
  # .............................................................
  
  # CAR prior for spatial random effect
  s[1:ncell] ~ dcar_normal(adj[1:k], weights[1:k], num[1:ncell], tau)
  # precision of CAR prior
  tau ~ dgamma(1, 1)
  
  # intercept of SDM
  alpha ~ dlogis(0, 1)
  # regression coefficient for canopy cover
  b_forest ~ dlogis(1, 1)
  # regression coefficient for impervious cover
  b_imperv ~ dlogis(-1, 1)
  # detection probability
  p ~ dunif(0, 1)
  
  # .............................................................
  # LIKELIHOOD
  # .............................................................
  
  # SDM - model for the latent state
  for(i in 1:ncell){
    log(lambda[i]) <- min(s[i] + alpha + b_forest*forest[i] + b_imperv*imperv[i], 10)
    psi[i] <- 1 - exp(-lambda[i])
    z[i] ~ dbern(psi[i])
  }
  
  # Camera submodel
  for(j in 1:nsite){
    muy[j] <- z[cell[j]]*p
    y[j] ~ dbin(muy[j], nsurveys[j])
  }
  
  # Harvest submodel
  for(i in 1:ncounty){
    # e: effort coefficient
    # low & high: indices to account for number of cells within each county
    mu_p[i] <- e[i]*inprod(lambda[low[i]:high[i]], z[low[i]:high[i]])
    w[i] ~ dpois(mu_p[i])
  }
  
} )

# .......................................................................
# PREPARE MODEL TO RUN
# .......................................................................

# function to provide random initial values for parameters
inits <- function() {
  base::list(z = rep(1, constants$ncell),
             alpha = runif(1,-1,1), 
             b_forest = runif(1, -1, 1),
             b_imperv = runif(1, -1, 1),
             tau = rgamma(1, 1, 1), 
             s = rep(0, base::length(data$num)),
             p = runif(1, 0, 1)) }

# parameters to monitor
keepers <- c('psi', 'p', "alpha", 'b_forest', "b_imperv")

nc <- 3 # number of chains
nb <- 1000 # number of initial MCMC iterations to discard
ni <- nb + 1000 # total number  of iterations
nt <- 5 # thinning interval

# .......................................................................
# PREPARE MODEL TO RUN
# .......................................................................

# create model
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

# .......................................................................
# INSPECT RESULTS
# .......................................................................

# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc(lapply(samples, coda::mcmc))

# Look at traceplots of the first four parameters
par(mfrow=c(2,2))
coda::traceplot(samples_mcmc[, 1:4])
# calculate Rhat convergence diagnostic for first four parameters
coda::gelman.diag(samples_mcmc[,1:4])

# extract mean and SD occurrence probability of each grid cell
samples <- do.call(rbind, samples_mcmc)
psi <- samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'psi\\['))]
psi_mean <- apply(psi, 2, mean)
psi_sd <- apply(psi, 2, sd)

# mean occurrence probability of the cells that contain testing data
test_psi <- psi_mean[test$INDEX]

# compute AUC to evaluate predictive performance
pROC::roc(response = test$Y, predictor = test_psi, ci = TRUE)

# map mean and standard deviation of occurrence probability
par(mfrow=c(1,1))
cov$psi_mean <- psi_mean
plot(cov["psi_mean"])
cov$psi_sd <- psi_sd
plot(cov["psi_sd"])                  
