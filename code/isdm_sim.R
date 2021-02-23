# simulation study - what are the consequences of poorly-measured effort data for iSDMs?

library(raster)
library(tidyverse)
library(sf)
library(nimble)
library(coda)
library(parallel)
library(patchwork)

# set seed for entire session
addTaskCallback(function(...) {set.seed(1);TRUE})

# set resolution of the fine grid
nrow_f <- 36
ncol_f <- 36

# set resolution of the coarse grid
nrow_c <- sqrt(nrow_f)
ncol_c <- sqrt(ncol_f)

# per-individual harvest probability (ph) will vary randomly by coarse cell ("county")
# each cell has probability of harvest somewhere in the 0-0.5 range
# this is intended to mimic reality, since all individuals in a population are not harvested 
ph_true <- runif(n = nrow_c*ncol_c, min = 0, max = 0.5)

# true probability is of course unknown
# instead, we rely on measures of harvest effort to inform probability of harvest
# assuming that increased effort will translate to higher probability of harvest
# we will simulate four effort scenarios:
# 1 - mediocre effort, which is modestly correlated with ph_true;
# 2 - moderate effort, which is fairly strongly correlated with ph_true;
# 3 - quality effort, which is strongly correlated with ph_true
# 4 - perfect effort, which is exactly equal to ph_true (unrealistic, but humor me)
noise1 <- ph_true + rnorm(n = length(ph_true), mean = 0, sd = 0.35)
effort1 <- (1 - 0.01)*((noise1 - min(noise1))/(max(noise1) - min(noise1))) + 0.01
cor(ph_true, effort1) # correlation = 0.52

noise2 <- ph_true + rnorm(n = length(ph_true), mean = 0, sd = 0.15)
effort2 <- (1 - 0.01)*((noise2 - min(noise2))/(max(noise2) - min(noise2))) + 0.01
cor(ph_true, effort2) # correlation = 0.77

noise3 <- ph_true + rnorm(n = length(ph_true), mean = 0, sd = 0.05)
effort3 <- (1 - 0.01)*((noise3 - min(noise3))/(max(noise3) - min(noise3))) + 0.01
cor(ph_true, effort3) # correlation = 0.95

effort4 <- ph_true
cor(ph_true, effort4) # correlation = 1 :)

# key for the simulation scenarios
# 2 abundance levels
# abun_int is to be used in function to generate true abundance in fine cells
# camera-only, harvest-only, and four integrated models for each abundancec level
# different effort quality for the four integrated models
scenarios <- tibble(
  scenario    = paste0("Scenario", 1:12),
  abundance   = c(rep("Low", 6), rep("High", 6)),
  abun_int    = c(rep(0, 6), rep(1, 6)),
  data        = c(rep(c("Camera", "Harvest", rep("Both", 4)), 2)),
  effort_name = c(rep(c(NA, "Moderate", "Mediocre", "Moderate", "Quality", "Perfect"), 2)),
  effort      = rep(c(list(tibble(id_c = 1:(nrow_c*ncol_c), ph = rep(NA, nrow_c*ncol_c), effort = rep(1, nrow_c*ncol_c))),
                      list(tibble(id_c = 1:(nrow_c*ncol_c), ph = ph_true, effort = effort2)),
                      list(tibble(id_c = 1:(nrow_c*ncol_c), ph = ph_true, effort = effort1)),
                      list(tibble(id_c = 1:(nrow_c*ncol_c), ph = ph_true, effort = effort2)),
                      list(tibble(id_c = 1:(nrow_c*ncol_c), ph = ph_true, effort = effort3)),
                      list(tibble(id_c = 1:(nrow_c*ncol_c), ph = ph_true, effort = effort4))), 2))

# per-individual detection probability by cameras within fine-resolution cells
pf <- 0.1

# what percentage of fine-resolution cells are NOT surveyed?
percent_na <- 0.85

# slope of linear relationship between covariate and occurrence
b1 <- -0.5
b2 <- 1.25

# MCMC settings - burnin and total number of iterations
nb <- 45000
ni <- nb + 5000

# fine-resolution grid
fine <- raster(nrows = nrow_f, ncols = ncol_f, xmn = -1, xmx = 1, ymn = -1, ymx = 1)

# coarse-resolution grid
coarse <- raster(nrows = nrow_c, ncols = ncol_c, xmn = -1, xmx = 1, ymn = -1, ymx = 1)

# coarse-resolution grid in polygon form
coarse_poly <- rasterToPolygons(coarse) %>% 
  sf::st_as_sf(.,  crs = 4326, agr = "constant") %>% 
  mutate(id_c = row_number()) %>% # county id
  dplyr::select(id_c, geometry)

# identify rows that should NOT be surveyed
na_rows <- tibble(
  id_f    = 1:ncell(fine),
  cell_na = rbinom(ncell(fine), 1, percent_na))

lambda_results <- list()
p_beta_results <- list()
for(i in 1:nrow(scenarios)){

    (  final <- rasterToPolygons(fine) %>% 
       sf::st_as_sf(., crs = 4326, agr = "constant") %>% 
       mutate(id_f       = 1:nrow(.),
              xcoord     = st_coordinates(st_centroid(.))[,1],
              ycoord     = st_coordinates(st_centroid(.))[,2],
              noise      = rnorm(n = nrow(.), mean = 0, 0.1),
              log_lambda = scenarios[[i, "abun_int"]] + b1*xcoord + b2*ycoord + noise, 
              lambda_f   = exp(log_lambda)) %>% 
       group_by(xcoord, ycoord) %>% 
       mutate(n_f     = rpois(n = 1, lambda = lambda_f),
              survey1 = rbinom(n = 1, size = n_f, prob = pf), 
              survey2 = rbinom(n = 1, size = n_f, prob = pf), 
              survey3 = rbinom(n = 1, size = n_f, prob = pf),
              survey4 = rbinom(n = 1, size = n_f, prob = pf), 
              survey5 = rbinom(n = 1, size = n_f, prob = pf), 
              survey6 = rbinom(n = 1, size = n_f, prob = pf),
              survey7 = rbinom(n = 1, size = n_f, prob = pf), 
              survey8 = rbinom(n = 1, size = n_f, prob = pf), 
              survey9 = rbinom(n = 1, size = n_f, prob = pf),
              survey10 = rbinom(n = 1, size = n_f, prob = pf), 
              survey11 = rbinom(n = 1, size = n_f, prob = pf), 
              survey12 = rbinom(n = 1, size = n_f, prob = pf),
              survey13 = rbinom(n = 1, size = n_f, prob = pf), 
              survey14 = rbinom(n = 1, size = n_f, prob = pf), 
              survey15 = rbinom(n = 1, size = n_f, prob = pf)) %>% 
       pivot_longer(survey1:survey15, names_to = "survey", values_to = "count") %>% 
       # "degrade" counts to detection/nondetection
       mutate(count = ifelse(count > 0, 1, 0)) %>% 
       full_join(na_rows) %>% 
       # convert not-surveyed cells to NA
       mutate(count = ifelse(cell_na == 1, NA, count)) %>% 
       pivot_wider(names_from = survey, values_from = count) %>%
       ungroup(.) %>%   
       sf::st_as_sf(., crs = 4326, agr = "constant") %>%   
       sf::st_join(coarse_poly, join = st_within) %>% 
       arrange(id_c) %>%
       group_by(id_c) %>%
       mutate(n_c = sum(n_f)) %>%  # coarse-cell abundance
       full_join(scenarios[[i, "effort"]][[1]]) %>% # add in the effort information
       mutate(harvest = rbinom(n = 1, size = n_c, ph)) %>% 
       ungroup(.) %>%
       mutate(cell_id = row_number()) %>%
       dplyr::select(cell_id, id_f, id_c, xcoord, ycoord, lambda_f, n_f,
                     survey1:survey15, n_c, ph, effort, harvest, geometry) )
  
  # function to calculate neighbor structure of grid (for spatial smoothing)
  st_queen <- function(a, b = a) st_relate(a, b, pattern = "F***T****")
  
  # for each cell, who are the neighboring grid cells?
  adj <- st_queen(dplyr::select(final, geometry))
  
  # for each cell in the fine grid - how many neighboring grid cells does it have?
  num <- sapply(adj, length)
  adj_final <- unlist(adj)
  
  # prepare the coarse-cell data
  # low and high are indices to tell the model which county each fine cell falls inside
  ( coarse_indices <- final %>%
      dplyr::select(id_f, id_c, effort, harvest) %>%
      st_drop_geometry(.) %>%
      ungroup(.) %>%
      mutate(row = row_number()) %>%
      group_by(id_c) %>%
      mutate(low = min(row), 
             high = max(row)) %>%
      dplyr::select(-id_f, -row) %>%
      distinct(.) )
  
  # prepare the fine-resolution data
  surveyed_cells <- final %>%
    filter(!is.na(survey1)) %>%
    ungroup(.) %>%
    dplyr::select(id_f, survey1:survey15, cell_id) %>%
    sf::st_drop_geometry(.)
  
  # camera-only model
  if(scenarios[[i, "data"]] == "Camera") {
    # 
    # package up the data for the model
    data <- list(
      xcoord  = as.vector(final$xcoord), # this is the x-coordinate of each cell, one of our covariates
      ycoord  = as.vector(final$ycoord), # and the y-coord
      num     = num, # number of adjacencies
      adj     = adj_final, # who's who - cell neighbors
      weights = rep(1, length(adj_final)), # weights for spatial random effect, set to 1
      y       = unname(as.matrix(dplyr::select(surveyed_cells, starts_with("survey")))))
    
    constants <- list(
      ncellf         = length(data$num), # number of fine-res cells
      nadj           = length(adj_final), # number of cell adjacencies
      nsurveyedcells = nrow(data$y), # number of cells with survey data
      nsurveys       = ncol(data$y), # number of surveys
      cell           = surveyed_cells$cell_id) # cell ID of the surveyed cells
    
    code <- nimble::nimbleCode( {
      
      # Priors
      # Spatial smoothing
      s[1:ncellf] ~ dcar_normal(adj[1:nadj], weights[1:nadj], num[1:ncellf], tau)
      tau ~ dgamma(1, 1) # scalar precision for SRE
      b1  ~ dnorm(0, sd = 2) # coefficient for x-coord
      b2  ~ dnorm(0, sd = 2) # coefficient for y-coord
      p   ~ dbeta(1, 1)  # detection probability
      
      # Likelihood
      # SDM - all fine-res cells
      for(i in 1:ncellf){
        log(lambda[i]) <- s[i] + b1*xcoord[i] + b2*ycoord[i]
        n[i] ~ dpois(lambda[i]) # latent true abundance
      }
      
      # camera submodel - detection prob is held constant
      for(j in 1:nsurveyedcells){
        for(k in 1:nsurveys){
          muy[j, k] <- 1 - pow(1 - p, n[cell[j]])
          y[j, k] ~ dbern(muy[j, k])
        }}
      
    })
    
    # funtion for initial values
    car_inits <- function() {
      base::list(n   = rep(1, constants$ncellf),
                 tau = rgamma(1, 1, 1),
                 b1  = runif(1, -1, 1),
                 b2  = runif(1, -1, 1),
                 s   = rep(0, base::length(data$num)),
                 p   = runif(1, 0, 1))}
    
    # parameters that we want to get results for
    keepers <- c("lambda", "p", "b1", "b2")
    
    # run the model
    out <- nimbleMCMC(code      = code,
                      constants = constants, 
                      data      = data, 
                      inits     = car_inits(), 
                      monitors  = keepers, 
                      nb        = nb, 
                      ni        = ni)
    
    lambda_result <- as_tibble(out) %>%
      dplyr::select(starts_with("lambda")) %>%
      pivot_longer(1:constants$ncellf, names_to = "cell", values_to = "lambda_est") %>%
      mutate(cell_id = extract_numeric(cell)) %>%
      group_by(cell_id) %>%
      summarise(lambda_est = median(lambda_est)) %>%
      full_join(final) %>%
      add_column(data        = scenarios[[i, "data"]],
                 effort_name = scenarios[[i, "effort_name"]],
                 abundance   = scenarios[[i, "abundance"]]) %>%
      dplyr::select(abundance, data, effort_name, cell_id, id_f, id_c, lambda_est, lambda_actual = lambda_f, geometry)
    
    p_beta_result <- as_tibble(out) %>%
      dplyr::select("p", "b1", "b2") %>%
      pivot_longer(p:b2, names_to = "parameter", values_to = "value") %>%
      group_by(parameter) %>%
      summarise(lower = quantile(value, 0.025),
                mean  = mean(value),
                upper = quantile(value, 0.975)) %>%
      add_column(data        = scenarios[[i, "data"]],
                 effort_name = scenarios[[i, "effort_name"]],
                 abundance   = scenarios[[i, "abundance"]]) %>% 
      dplyr::select(abundance, data, effort_name, parameter, mean, lower, upper)
    
    lambda_results[[i]] <- lambda_result
    p_beta_results[[i]] <- p_beta_result
    
    print(paste("Finished", i, "of", nrow(scenarios)))
    
  } else if(scenarios[[i, "data"]] == "Harvest"){
    
    # package up the data for the model
    data <- list(
      xcoord  = as.vector(final$xcoord), # this is the x-coordinate of each cell, one of our covariates
      ycoord  = as.vector(final$ycoord), # and the y-coord
      num     = num, # number of adjacencies
      adj     = adj_final, # who's who - cell neighbors
      weights = rep(1, length(adj_final)), # weights for spatial random effect, set to 1
      harvest = coarse_indices$harvest,
      effort  = coarse_indices$effort)
    
    constants <- list(
      ncellf = length(data$num), # number of fine-res cells
      nadj   = length(adj_final), # number of cell adjacencies
      ncellc = nrow(coarse_indices), 
      low    = coarse_indices$low, 
      high   = coarse_indices$high)
    
    code <- nimble::nimbleCode( {
      
      # Priors
      # Spatial smoothing
      s[1:ncellf] ~ dcar_normal(adj[1:nadj], weights[1:nadj], num[1:ncellf], tau)
      tau ~ dgamma(1, 1) # scalar precision for SRE
      b1  ~ dnorm(0, sd = 2) # coefficient for x-coord
      b2  ~ dnorm(0, sd = 2) # coefficient for y-coord
      p   ~ dbeta(1, 1)  # detection probability
      
      gamma0 ~ dnorm(0, sd = 2)
      gamma1 ~ dnorm(0, sd = 2)
      
      # Likelihood
      # SDM - all fine-res cells
      for(i in 1:ncellf){
        log(lambda[i]) <- s[i] + b1*xcoord[i] + b2*ycoord[i]
        n[i] ~ dpois(lambda[i]) # latent true abundance
      }
      
      # harvest submodel
      for(q in 1:ncellc){
        # function to scale fine-scale to county-scale expected abundance
        log(lambda_county[q]) <- gamma0 + gamma1*log(sum(lambda[low[q]:high[q]])) 
        harvest[q] ~ dpois(effort[q]*lambda_county[q])
      }
    })
    
    # funtion for initial values
    car_inits <- function() {
      base::list(n      = rep(1, constants$ncellf),
                 tau    = rgamma(1, 1, 1),
                 b1     = runif(1, -1, 1),
                 b2     = runif(1, -1, 1),
                 gamma0 = runif(1, -1, 1), 
                 gamma1 = runif(1, -1, 1),
                 s      = rep(0, base::length(data$num)),
                 p      = runif(1, 0, 1))}
    
    # parameters that we want to get results for
    keepers <- c("lambda", "p", "b1", "b2")
    
    # run the model
    out <- nimbleMCMC(code      = code,
                      constants = constants, 
                      data      = data, 
                      inits     = car_inits(), 
                      monitors  = keepers, 
                      nb        = nb, 
                      ni        = ni)
    
    lambda_result <- as_tibble(out) %>%
      dplyr::select(starts_with("lambda")) %>%
      pivot_longer(1:constants$ncellf, names_to = "cell", values_to = "lambda_est") %>%
      mutate(cell_id = extract_numeric(cell)) %>%
      group_by(cell_id) %>%
      summarise(lambda_est = median(lambda_est)) %>%
      full_join(final) %>%
      add_column(data        = scenarios[[i, "data"]],
                 effort_name = scenarios[[i, "effort_name"]],
                 abundance   = scenarios[[i, "abundance"]]) %>%
      dplyr::select(abundance, data, effort_name, cell_id, id_f, id_c, lambda_est, lambda_actual = lambda_f, geometry)
    
    p_beta_result <- as_tibble(out) %>%
      dplyr::select("p", "b1", "b2") %>%
      pivot_longer(p:b2, names_to = "parameter", values_to = "value") %>%
      group_by(parameter) %>%
      summarise(lower = quantile(value, 0.025),
                mean  = mean(value),
                upper = quantile(value, 0.975)) %>%
      add_column(data        = scenarios[[i, "data"]],
                 effort_name = scenarios[[i, "effort_name"]],
                 abundance   = scenarios[[i, "abundance"]]) %>% 
      dplyr::select(abundance, data, effort_name, parameter, mean, lower, upper)
    
    lambda_results[[i]] <- lambda_result
    p_beta_results[[i]] <- p_beta_result
    
    print(paste("Finished", i, "of", nrow(scenarios)))
    
  } else {
    
    # package up the data for the model
    data <- list(
      xcoord  = as.vector(final$xcoord), # this is the x-coordinate of each cell, one of our covariates
      ycoord  = as.vector(final$ycoord), # and the y-coord
      num     = num, # number of adjacencies
      adj     = adj_final, # who's who - cell neighbors
      weights = rep(1, length(adj_final)), # weights for spatial random effect, set to 1
      y       = unname(as.matrix(dplyr::select(surveyed_cells, starts_with("survey")))),
      harvest = coarse_indices$harvest,
      effort  = coarse_indices$effort)
    
    constants <- list(
      ncellf         = length(data$num), # number of fine-res cells
      nadj           = length(adj_final), # number of cell adjacencies
      nsurveyedcells = nrow(data$y), # number of cells with survey data
      nsurveys       = ncol(data$y), # number of surveys
      cell           = surveyed_cells$cell_id,
      ncellc         = nrow(coarse_indices), 
      low            = coarse_indices$low, 
      high           = coarse_indices$high) # cell ID of the surveyed cells
    
    code <- nimble::nimbleCode( {
      
      # Priors
      # Spatial smoothing
      s[1:ncellf] ~ dcar_normal(adj[1:nadj], weights[1:nadj], num[1:ncellf], tau)
      tau ~ dgamma(1, 1) # scalar precision for SRE
      b1  ~ dnorm(0, sd = 2) # coefficient for x-coord
      b2  ~ dnorm(0, sd = 2) # coefficient for y-coord
      p   ~ dbeta(1, 1)  # detection probability
      gamma0 ~ dnorm(0, sd = 2)
      gamma1 ~ dnorm(0, sd = 2)
      
      # Likelihood
      # SDM - all fine-res cells
      for(i in 1:ncellf){
        log(lambda[i]) <- s[i] + b1*xcoord[i] + b2*ycoord[i]
        n[i] ~ dpois(lambda[i]) # latent true abundance
      }
      
      # camera submodel - detection prob is held constant
      for(j in 1:nsurveyedcells){
        for(k in 1:nsurveys){
          muy[j, k] <- 1 - pow(1 - p, n[cell[j]])
          y[j, k] ~ dbern(muy[j, k])
        }}
      
      for(q in 1:ncellc){
        # function to scale fine-scale to county-scale expected abundance
        log(lambda_county[q]) <- gamma0 + gamma1*log(sum(lambda[low[q]:high[q]])) 
        harvest[q] ~ dpois(effort[q]*lambda_county[q])
      }
    })
    
    # funtion for initial values
    car_inits <- function() {
      base::list(n      = rep(1, constants$ncellf),
                 tau    = rgamma(1, 1, 1),
                 b1     = runif(1, -1, 1),
                 b2     = runif(1, -1, 1),
                 gamma0 = runif(1, -1, 1), 
                 gamma1 = runif(1, -1, 1),
                 s      = rep(0, base::length(data$num)),
                 p      = runif(1, 0, 1))}
    
    # parameters that we want to get results for
    keepers <- c("lambda", "p", "b1", "b2")
    
    # run the model
    out <- nimbleMCMC(code      = code,
                      constants = constants, 
                      data      = data, 
                      inits     = car_inits(), 
                      monitors  = keepers, 
                      nb        = nb, 
                      ni        = ni)
    
    lambda_result <- as_tibble(out) %>%
      dplyr::select(starts_with("lambda")) %>%
      pivot_longer(1:constants$ncellf, names_to = "cell", values_to = "lambda_est") %>%
      mutate(cell_id = extract_numeric(cell)) %>%
      group_by(cell_id) %>%
      summarise(lambda_est = median(lambda_est)) %>%
      full_join(final) %>%
      add_column(data        = scenarios[[i, "data"]],
                 effort_name = scenarios[[i, "effort_name"]],
                 abundance   = scenarios[[i, "abundance"]]) %>%
      dplyr::select(abundance, data, effort_name, cell_id, id_f, id_c, lambda_est, lambda_actual = lambda_f, geometry)
    
    p_beta_result <- as_tibble(out) %>%
      dplyr::select("p", "b1", "b2") %>%
      pivot_longer(p:b2, names_to = "parameter", values_to = "value") %>%
      group_by(parameter) %>%
      summarise(lower = quantile(value, 0.025),
                mean  = mean(value),
                upper = quantile(value, 0.975)) %>%
      add_column(data        = scenarios[[i, "data"]],
                 effort_name = scenarios[[i, "effort_name"]],
                 abundance   = scenarios[[i, "abundance"]]) %>% 
      dplyr::select(abundance, data, effort_name, parameter, mean, lower, upper)
    
    lambda_results[[i]] <- lambda_result
    p_beta_results[[i]] <- p_beta_result
    
    print(paste("Finished", i, "of", nrow(scenarios)))
  }
}

lambs <- bind_rows(lambda_results) %>% 
  mutate(data = factor(data, levels = c("Camera", "Harvest", "Both"))) %>% 
  mutate(name = paste(data, effort_name))

betas <- bind_rows(p_beta_results)

( labDat <- lambs %>% 
    mutate(data = as.character(data)) %>% 
    mutate(data = ifelse(data == "Both", "Int.", data)) %>% 
    mutate(name = paste(data, effort_name)) %>% 
    mutate(name = ifelse(name == "Camera NA", "Camera", name)) %>% 
    mutate(name = ifelse(name == "Harvest Moderate", "Harvest", name)) %>% 
    mutate(name = factor(name, levels = c("Camera", "Harvest", "Int. Mediocre", 
                                          "Int. Moderate", "Int. Quality", "Int. Perfect"))) %>%
    mutate(abundance = ifelse(abundance == "Low", "Low abundance", "High abundance")) %>% 
    mutate(abundance = factor(abundance, levels = c("Low abundance", "High abundance"))) %>% 
    group_by(abundance, name) %>%
    mutate(lambCor = cor(lambda_est, lambda_actual)) %>%
    mutate(bias = lambda_est - lambda_actual) %>% 
    dplyr::select(abundance, name, bias, lambCor) %>%
    summarise(bias = mean(bias), 
              lambCor = mean(lambCor)) %>% 
    add_column(x = 8, 
               y = 32) %>% 
    mutate(lambCor = as.character(round(lambCor, 2)),
           bias = as.character(round(bias, 2))) %>% 
    dplyr::select(x, y, abundance, name, lambCor, bias) )

( lamb_plot <- lambs %>% 
    mutate(data = as.character(data)) %>% 
    mutate(data = ifelse(data == "Both", "Int.", data)) %>% 
    mutate(name = paste(data, effort_name)) %>% 
    mutate(name = ifelse(name == "Camera NA", "Camera", name)) %>% 
    mutate(name = ifelse(name == "Harvest Moderate", "Harvest", name)) %>% 
    mutate(name = factor(name, levels = c("Camera", "Harvest", "Int. Mediocre", 
                                          "Int. Moderate", "Int. Quality", "Int. Perfect"))) %>%
    mutate(abundance = ifelse(abundance == "Low", "Low abundance", "High abundance")) %>% 
    mutate(abundance = factor(abundance, levels = c("Low abundance", "High abundance"))) %>% 
    ggplot(aes(x = lambda_actual, y = lambda_est)) + 
    geom_point(alpha = 0.05) + 
    facet_grid(abundance ~ name) +
    geom_text(aes(x = x, y = y, label = paste("Bias:", bias)),
              data = labDat, size = 2.5) +
    geom_text(aes(x = x, y = 38, label = paste("Corr:", lambCor)), 
              data = labDat, size = 2.5) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    theme_classic() + 
    labs(x = expression("Actual  "*lambda[i]*"\n"),
         y = expression("Estimated  "*lambda[i]),
         title = "(a)") +
    theme(strip.background = element_rect(fill = "gray90", color = NA),
          plot.title = element_text(size = 8), 
          strip.text = element_text(size = 7), 
          axis.title = element_text(size = 8), 
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.1),
          axis.ticks = element_line(size = 0.1)) )

( betas2 <- betas %>% 
    filter(!parameter == "p") %>% 
    mutate(label = ifelse(parameter == "b1", "beta[1]", "beta[2]")) )

( betas_plot <- betas %>% 
  filter(!parameter == "p") %>% 
  mutate(label = ifelse(parameter == "b1", "beta[1]", "beta[2]")) %>% 
  mutate(data = ifelse(data == "Both", "Int.", data)) %>% 
  mutate(name = paste(data, effort_name)) %>% 
  mutate(name = ifelse(name == "Camera NA", "Camera", name)) %>% 
  mutate(name = ifelse(name == "Harvest Moderate", "Harvest", name)) %>% 
  mutate(name = factor(name, levels = c("Int. Perfect",
                                        "Int. Quality", 
                                        "Int. Moderate", 
                                        "Int. Mediocre", 
                                        "Harvest", 
                                        "Camera"))) %>%
  ggplot(aes(x = mean, y = name)) + 
  geom_errorbar(aes(xmin = lower, xmax = upper), size = 1, width = 0) +
  geom_point(size = 2) +
  facet_grid(factor(abundance, levels = c("Low", "High")) ~ label,
             scales = "free", labeller = label_parsed) +
  geom_vline(xintercept = 0, color = "gray60", linetype = 3) +
  geom_vline(data = filter(betas2, label == "beta[1]"), aes(xintercept = b1), color = "red") +
  geom_vline(data = filter(betas2, label == "beta[2]"), aes(xintercept = b2), color = "red") + 
  theme_classic() + 
  labs(x = "Effect",
       title = "(b)") +
  theme(strip.background = element_rect(fill = "gray90", color = NA),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1),
        plot.title = element_text(size = 8), 
        strip.text = element_text(size = 7), 
        axis.title = element_text(size = 8), 
        axis.text = element_text(size = 6)) )

# Fig. 7
( sim_results <- lamb_plot + betas_plot + plot_layout(widths = c(3, 1)) )
