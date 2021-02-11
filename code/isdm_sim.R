# simulation study - what are the consequences of poorly-measured effort data for iSDMs?

library(raster)
library(tidyverse)
library(sf)
library(nimble)
library(coda)
library(parallel)

# set resolution of the fine grid
nrowF <- 36
ncolF <- 36

# set resolution of the coarse grid
nrowC <- sqrt(nrowF)
ncolC <- sqrt(ncolF)

# harvest detection probabilities for scenarios in which probability of harvest varies among coarse cells
set.seed(1)
# each cell has probability of harvest somewhere in the 0-0.75 range
pC_het <- runif(nrowC*ncolC, 0, 0.5)
set.seed(1)

# now, for the case in which there is harvest heterogeneity and we have imperfect, but correlated,
# measures of harvest effort
# add some noise and scale to fall between 0.01 and 1 (pC = 0 causes problems for model)
pC_het_noise1 <- pC_het + rnorm(length(pC_het), 0, 0.5)
pC_cor1 <- (1 - 0.01)*((pC_het_noise1 - min(pC_het_noise1))/(max(pC_het_noise1) - min(pC_het_noise1))) + 0.01

set.seed(1)
pC_het_noise2 <- pC_het + rnorm(length(pC_het), 0, 0.1)
pC_cor2 <- (1 - 0.01)*((pC_het_noise2 - min(pC_het_noise2))/(max(pC_het_noise2) - min(pC_het_noise2))) + 0.01

# key for the simulation scenarios
scenarios <- tibble(
  Scenario = paste0("Scenario", 1:6),
  Description = c("Camera only", 
                  "Harvest only, mediocre effort",
                  "Perfect harvest",
                  "Mediocre effort",
                  "Quality effort",
                  "Perfect effort"),
  pC = c(list(tibble(idC = 1:(nrowC*ncolC),
                     pC = rep(NA,nrowC*ncolC),
                     pCMod = rep(1, nrowC*ncolC))),
         list(tibble(idC = 1:(nrowC*ncolC),
                     pC = pC_het,
                     pCMod = pC_cor1)),
         list(tibble(idC = 1:(nrowC*ncolC),
                     pC = pC_het,
                     pCMod = rep(1, nrowC*ncolC))),
         list(tibble(idC = 1:(nrowC*ncolC),
                     pC = pC_het,
                     pCMod = pC_cor1)),
         list(tibble(idC = 1:(nrowC*ncolC),
                     pC = pC_het,
                     pCMod = pC_cor2)),
         list(tibble(idC = 1:(nrowC*ncolC),
                     pC = pC_het,
                     pCMod = pC_het))))

# fine-resolution detection probability
pF <- 0.25

# what percentage of fine-resolution cells are NOT surveyed?
percent_NA <- 0.75

# slope of linear relationship between covariate and occurrence
slope <- 1.5

# number of MCMC iterations
nb <- 5000
ni <- nb + 1000
nc <- 3

# fine-resolution grid
fine <- raster(nrows = nrowF, ncols = ncolF, xmn = -1, xmx = 1, ymn = -1, ymx = 1)

coarse <- raster(nrows = nrowC, ncols = ncolC, xmn = -1, xmx = 1, ymn = -1, ymx = 1)

set.seed(1)
fine_df <- rasterToPoints(fine) %>% 
  as_tibble() %>% 
  mutate(idF = row_number(),
         noise = rnorm(nrow(.), 0, 0.1),
         log_lambda = x*slope + noise,
         lambda = exp(log_lambda)) %>% 
  group_by(x, y) %>% 
  mutate(nF = rpois(1, lambda),
         yF1 = rbinom(1, nF, pF), 
         yF2 = rbinom(1, nF, pF), 
         yF3 = rbinom(1, nF, pF)) %>% 
  pivot_longer(yF1:yF3, names_to = "survey", values_to = "count") %>% 
  mutate(count = ifelse(count > 0, 1, 0)) %>% 
  pivot_wider(names_from = survey, values_from = count)

surveyStack <- stack(setValues(fine, fine_df$idF), 
                     setValues(fine, fine_df$x),
                     setValues(fine, fine_df$lambda),
                     setValues(fine, fine_df$nF), 
                     setValues(fine, fine_df$yF1), 
                     setValues(fine, fine_df$yF2), 
                     setValues(fine, fine_df$yF3)) 

names(surveyStack) <- c("idF", "canopy", "lambda", "nF", "yF1", "yF2", "yF3")

fine_poly <- rasterToPolygons(surveyStack) %>% 
  sf::st_as_sf(.)

coarse_poly <- rasterToPolygons(coarse) %>% 
  sf::st_as_sf(.) %>% 
  mutate(idC = row_number()) %>% 
  dplyr::select(idC, geometry)

na_rows <- tibble(
  idF = 1:nrow(fine_poly),
  cellNA = rbinom(nrow(fine_poly), 1, percent_NA))

lambda_results <- list()
p_beta_results <- list()
for(i in 1:nrow(scenarios)){
  
  final <- fine_poly %>%
    st_join(coarse_poly, join = st_within) %>%
    arrange(idC) %>%
    group_by(idC) %>%
    mutate(nC = sum(nF)) %>%
    full_join(scenarios[[i, "pC"]][[1]]) %>%
    mutate(yC = rbinom(1, nC, pC)) %>%
    full_join(na_rows) %>%
    pivot_longer(yF1:yF3, names_to = "survey", values_to = "value") %>%
    mutate(value = ifelse(cellNA == 1, NA, value)) %>%
    pivot_wider(names_from = survey, values_from = value) %>%
    ungroup(.) %>%
    mutate(cellID = row_number()) %>%
    dplyr::select(cellID, idF, idC, canopy, lambda, nF, yF1:yF3, nC, pC, pCMod, yC, geometry) %>%
    sf::st_as_sf(.)

  # function to calculate neighbor structure of grid (for spatial smoothing)
  st_queen <- function(a, b = a) st_relate(a, b, pattern = "F***T****")
  
  # for each cell, who are the neighboring grid cells?
  adj <- st_queen(dplyr::select(final, geometry))
  
  # for each cell in the fine grid - how many neighboring grid cells does it have?
  num <- sapply(adj, length)
  adj_final <- unlist(adj)
  
  ( coarse_indices <- final %>%
      dplyr::select(idF, idC, pCMod, yC) %>%
      st_drop_geometry(.) %>%
      ungroup(.) %>%
      mutate(row = row_number()) %>%
      group_by(idC) %>%
      mutate(low = min(row),
             high = max(row)) %>%
      dplyr::select(-idF, -row) %>%
      distinct(.) )
  
  surveyed_cells <- final %>%
    filter(!is.na(yF1)) %>%
    ungroup(.) %>%
    dplyr::select(idF, yF1:yF3, cellID) %>%
    sf::st_drop_geometry(.)

  if(scenarios[[i, "Scenario"]] == "Scenario1") {

    data <- list(
      canopy = final$canopy,
      num = num,
      adj = adj_final,
      weights = rep(1, length(adj_final)),
      y = unname(as.matrix(dplyr::select(surveyed_cells, starts_with("yF")))))
    
    constants <- list(
      nCellF = length(data$num),
      nAdj = length(adj_final),
      nSurveyedCells = nrow(data$y),
      nSurveys = ncol(data$y),
      cell = surveyed_cells$cellID)
    
    code <- nimble::nimbleCode( {
      
      s[1:nCellF] ~ dcar_normal(adj[1:nAdj], weights[1:nAdj], num[1:nCellF], tau)
      tau ~ dgamma(1, 1)
      b1 ~ dnorm(0, 2)
      p ~ dbeta(1, 1)
      
      for(i in 1:nCellF){
        log(lambda[i]) <- s[i] + b1*canopy[i]
        N[i] ~ dpois(lambda[i])
      }
      
      for(j in 1:nSurveyedCells){
        for(k in 1:nSurveys){
          muy[j, k] <- 1 - pow(1 - p, N[cell[j]])
          y[j, k] ~ dbern(muy[j, k])
        }}

    })
    
    car_inits <- function() {
      base::list(N = rep(1, constants$nCellF),
                 tau = rgamma(1, 1, 1),
                 b1 = runif(1, 0, 1),
                 s = rep(0, base::length(data$num)),
                 p = runif(1, 0, 1))}
    
    # parameters that we want to get results for
    keepers <- c("lambda", "p", "b1")
    
    cl <- makeCluster(nc)
    parallel::clusterExport(cl, c("code",
                                  "car_inits",
                                  "data",
                                  "constants",
                                  "keepers",
                                  "nb", 
                                  "ni"))
    
    for(j in seq_along(cl)) {
      set.seed(j)
      init <- car_inits()
      clusterExport(cl[j], "init")
    }
    
    out <- clusterEvalQ(cl, {
      library(nimble)
      library(coda)
      
      model <- nimbleModel(code = code,
                           name = "code", 
                           constants = constants, 
                           data = data, 
                           inits = init)
      
      Cmodel <- compileNimble(model)
      modelConf <- configureMCMC(model)
      modelConf$addMonitors(keepers)
      modelMCMC <- buildMCMC(modelConf)
      CmodelMCMC <- compileNimble(modelMCMC, project = model)
      out1 <- runMCMC(CmodelMCMC,
                      nburnin = nb,
                      niter = ni)
      return(as.mcmc(out1))
    })
    
    stopCluster(cl)
    
    out.mcmc <- as.mcmc(out)
    
    lambda_result <- do.call(rbind, out) %>%
      as_tibble() %>% 
      dplyr::select(starts_with("lambda")) %>% 
      pivot_longer(1:constants$nCellF, names_to = "cell", values_to = "lambdaEst") %>% 
      mutate(cellID = extract_numeric(cell)) %>% 
      group_by(cellID) %>% 
      summarise(lambdaEst = median(lambdaEst)) %>% 
      full_join(final) %>% 
      add_column(Scenario = scenarios[[i, "Scenario"]]) %>% 
      dplyr::select(Scenario, cellID, idF, idC, lambdaEst, lambdaActual = lambda, geometry)
    
    p_beta_result <- do.call(rbind, out) %>%
      as_tibble() %>% 
      dplyr::select("p", "b1") %>% 
      pivot_longer(p:b1, names_to = "parameter", values_to = "value") %>% 
      group_by(parameter) %>% 
      summarise(lower = quantile(value, 0.025), 
                mean = mean(value), 
                upper = quantile(value, 0.975)) %>% 
      add_column(Scenario = scenarios[[i, "Scenario"]])
    
    lambda_results[[i]] <- lambda_result
    p_beta_results[[i]] <- p_beta_result
    
  } else if(scenarios[[i, "Scenario"]] == "Scenario2"){
    
    data <- list(
      canopy = final$canopy,
      num = num, 
      adj = adj_final, 
      weights = rep(1, length(adj_final)),
      w = coarse_indices$yC,
      pC = coarse_indices$pCMod)
    
    constants <- list(
      nCellF = length(data$num), 
      nAdj = length(adj_final), 
      nCellC = nrow(coarse_indices),
      low = coarse_indices$low,
      high = coarse_indices$high)
    
    code <- nimble::nimbleCode( {
      
      s[1:nCellF] ~ dcar_normal(adj[1:nAdj], weights[1:nAdj], num[1:nCellF], tau)
      tau ~ dgamma(1, 1)
      b1 ~ dnorm(0, 2)
      p ~ dbeta(1, 1)
      
      for(i in 1:nCellF){
        log(lambda[i]) <- s[i] + b1*canopy[i]
        N[i] ~ dpois(lambda[i])
      }

      for(q in 1:nCellC){
        muP[q] <- pC[q]*sum(lambda[low[q]:high[q]])
        w[q] ~ dpois(muP[q])
      }
    })
    
    car_inits <- function() {
      base::list(N = rep(1, constants$nCellF),
                 tau = rgamma(1, 1, 1), 
                 b1 = runif(1, 0, 1),
                 s = rep(0, base::length(data$num)),
                 p = runif(1, 0, 1))}
    
    # parameters that we want to get results for
    keepers <- c("lambda", "p", "b1")
    
    cl <- makeCluster(nc)
    parallel::clusterExport(cl, c("code",
                                  "car_inits",
                                  "data",
                                  "constants",
                                  "keepers",
                                  "nb", 
                                  "ni"))
    
    for(j in seq_along(cl)) {
      set.seed(j)
      init <- car_inits()
      clusterExport(cl[j], "init")
    }
    
    out <- clusterEvalQ(cl, {
      library(nimble)
      library(coda)
      
      model <- nimbleModel(code = code,
                           name = "code", 
                           constants = constants, 
                           data = data, 
                           inits = init)
      
      Cmodel <- compileNimble(model)
      modelConf <- configureMCMC(model)
      modelConf$addMonitors(keepers)
      modelMCMC <- buildMCMC(modelConf)
      CmodelMCMC <- compileNimble(modelMCMC, project = model)
      out1 <- runMCMC(CmodelMCMC,
                      nburnin = nb,
                      niter = ni)
      return(as.mcmc(out1))
    })
    
    stopCluster(cl)
    
    out.mcmc <- as.mcmc(out)
    
    lambda_result <- do.call(rbind, out) %>%
      as_tibble() %>% 
      dplyr::select(starts_with("lambda")) %>% 
      pivot_longer(1:constants$nCellF, names_to = "cell", values_to = "lambdaEst") %>% 
      mutate(cellID = extract_numeric(cell)) %>% 
      group_by(cellID) %>% 
      summarise(lambdaEst = median(lambdaEst)) %>% 
      full_join(final) %>% 
      add_column(Scenario = scenarios[[i, "Scenario"]]) %>% 
      dplyr::select(Scenario, cellID, idF, idC, lambdaEst, lambdaActual = lambda, geometry)
    
    p_beta_result <- do.call(rbind, out) %>%
      as_tibble() %>% 
      dplyr::select("p", "b1") %>% 
      pivot_longer(p:b1, names_to = "parameter", values_to = "value") %>% 
      group_by(parameter) %>% 
      summarise(lower = quantile(value, 0.025), 
                mean = mean(value), 
                upper = quantile(value, 0.975)) %>% 
      add_column(Scenario = scenarios[[i, "Scenario"]])
    
    lambda_results[[i]] <- lambda_result
    p_beta_results[[i]] <- p_beta_result
    
  } else {
    
    data <- list(
      canopy = final$canopy,
      num = num,
      adj = adj_final,
      weights = rep(1, length(adj_final)),
      y = unname(as.matrix(dplyr::select(surveyed_cells, starts_with("yF")))),
      w = coarse_indices$yC,
      pC = coarse_indices$pCMod)
    
    constants <- list(
      nCellF = length(data$num),
      nAdj = length(adj_final),
      nSurveyedCells = nrow(data$y),
      nSurveys = ncol(data$y),
      cell = surveyed_cells$cellID,
      nCellC = nrow(coarse_indices),
      low = coarse_indices$low,
      high = coarse_indices$high)
    
    code <- nimble::nimbleCode( {
      
      s[1:nCellF] ~ dcar_normal(adj[1:nAdj], weights[1:nAdj], num[1:nCellF], tau)
      tau ~ dgamma(1, 1)
      b1 ~ dnorm(0, 2)
      p ~ dbeta(1, 1)
      
      for(i in 1:nCellF){
        log(lambda[i]) <- s[i] + b1*canopy[i]
        N[i] ~ dpois(lambda[i])
      }
      
      for(j in 1:nSurveyedCells){
        for(k in 1:nSurveys){
          muy[j, k] <- 1 - pow(1 - p, N[cell[j]])
          y[j, k] ~ dbern(muy[j, k])
        }}
      
      for(q in 1:nCellC){
        muP[q] <- pC[q]*sum(lambda[low[q]:high[q]])
        w[q] ~ dpois(muP[q])
      }
    })
    
    car_inits <- function() {
      base::list(N = rep(1, constants$nCellF),
                 tau = rgamma(1, 1, 1),
                 b1 = runif(1, 0, 1),
                 s = rep(0, base::length(data$num)),
                 p = runif(1, 0, 1))}
    
    # parameters that we want to get results for
    keepers <- c("lambda", "p", "b1")
    
    cl <- makeCluster(nc)
    parallel::clusterExport(cl, c("code",
                                  "car_inits",
                                  "data",
                                  "constants",
                                  "keepers",
                                  "nb", 
                                  "ni"))
    
    for(j in seq_along(cl)) {
      set.seed(j)
      init <- car_inits()
      clusterExport(cl[j], "init")
    }
    
    out <- clusterEvalQ(cl, {
      library(nimble)
      library(coda)
      
      model <- nimbleModel(code = code,
                           name = "code", 
                           constants = constants, 
                           data = data, 
                           inits = init)
      
      Cmodel <- compileNimble(model)
      modelConf <- configureMCMC(model)
      modelConf$addMonitors(keepers)
      modelMCMC <- buildMCMC(modelConf)
      CmodelMCMC <- compileNimble(modelMCMC, project = model)
      out1 <- runMCMC(CmodelMCMC,
                      nburnin = nb,
                      niter = ni)
      return(as.mcmc(out1))
    })
    
    stopCluster(cl)
    
    out.mcmc <- as.mcmc(out)
    
    lambda_result <- do.call(rbind, out) %>%
      as_tibble() %>% 
      dplyr::select(starts_with("lambda")) %>% 
      pivot_longer(1:constants$nCellF, names_to = "cell", values_to = "lambdaEst") %>% 
      mutate(cellID = extract_numeric(cell)) %>% 
      group_by(cellID) %>% 
      summarise(lambdaEst = median(lambdaEst)) %>% 
      full_join(final) %>% 
      add_column(Scenario = scenarios[[i, "Scenario"]]) %>% 
      dplyr::select(Scenario, cellID, idF, idC, lambdaEst, lambdaActual = lambda, geometry)
    
    p_beta_result <- do.call(rbind, out) %>%
      as_tibble() %>% 
      dplyr::select("p", "b1") %>% 
      pivot_longer(p:b1, names_to = "parameter", values_to = "value") %>% 
      group_by(parameter) %>% 
      summarise(lower = quantile(value, 0.025), 
                mean = mean(value), 
                upper = quantile(value, 0.975)) %>% 
      add_column(Scenario = scenarios[[i, "Scenario"]])
    
    lambda_results[[i]] <- lambda_result
    p_beta_results[[i]] <- p_beta_result
  }
}
  
lambs <- bind_rows(lambda_results) %>% 
  full_join(scenarios)

ggplot(lambs, aes(x = lambdaActual, y = lambdaEst)) + 
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1) + 
  facet_wrap(~Description) + 
  theme_classic()

pb <- bind_rows(p_beta_results) %>% 
  full_join(scenarios)

ggplot(filter(pb, parameter == "b1"), aes(x = mean, y = Description)) + 
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0) + 
  geom_point() +
  geom_vline(xintercept = 1.5, color = "red", size = 2) 