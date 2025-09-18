# install packages
install.packages(c("devtools", "tidyverse", "unmarked"))
devtools::install_github("inbo/camtraptor")
devtools::install_github("MarcusRowcliffe/camtrapDensity")

install.packages("lubridate")

# load packages
library(camtrapDensity)
library(dplyr)
library(purrr)
library(lubridate)

# Read and check data
pkg <- read_camtrapDP("./juliana-martins-20250409154914/datapackage.json")
plot_deployment_schedule(pkg)

locNameElements <- pkg$data$deployments$locationName %>%
  strsplit(" ")

pkg$data$deployments <- pkg$data$deployments %>%
  mutate(distance = map_chr(locNameElements, \(x) x[1]),
         transect = map_chr(locNameElements, \(x) x[3]),
         zone = map_chr(locNameElements, \(x) x[4]),
         phase = map_chr(locNameElements, \(x) x[6]))


library(dplyr)

#' Make a sequence of occasion time cuts
#' 
#' INPUT
#' start / end: matching vectors of POSIX deployment start / end date-times
#' interval: occasion intervals in days
#' start_hour: numeric hour of day at which to start the sequence
#' 
#' OUTPUT
#' A regular sequence of POSIX data-times spanning the range of start / end
#' 
make_cutSeq <- function(start, end, interval=7, start_hour=0){
  mn <- min(start) %>%
    as.POSIXlt()
  mx <- max(end)
  time <- mn$hour + mn$min/60 + mn$sec/3600
  mn <- if(time >= start_hour) 
    mn - 3600 * (time + start_hour) else
      mn - 3600 * (time + 24 - start_hour)
  seq(mn, mx+interval*86400, interval*86400)
}

#' Make an effort matrix
#' 
#' INPUT
#' deployments: a dataframe of deployment data with columns
#'    locationName: location identifiers, typically globally unique
#'    locationName: alternative location identifiers, typically shorter, locally unique
#'    start / end: POSIX date-times at which deployments start and end
#' cuts: a sequence of POSIX date-times defining detection occasions
#'    Generated internally using make_cutSeq if NULL
#' interval / start_hour: passed to make_cutSeq if cuts is NULL
#' 
#' OUTPUT
#' A list with elements:
#'    effort: a sites x occasions numeric matrix of effort in days
#'    cuts: a POSIX vector of the occasion date-time cutpoints
#'    
make_emat <- function(deployments, cuts = NULL,
                      interval = 7, start_hour = 0){
  if(is.null(cuts)) cuts <- make_cutSeq(deployments$start, 
                                        deployments$end, 
                                        interval = interval,
                                        start_hour = start_hour)
  intervals <- as.numeric(diff(cuts))
  if(any(intervals <= 0))
    stop("cuts are not continually increasing")
  if(min(deployments$start) < min(cuts) | max(deployments$end) > max(cuts))
    stop("cuts do not span start/end times")
  
  nocc <- length(intervals)
  ndep <- nrow(deployments)
  emat <- data.frame(loc = rep(deployments$locationName, nocc),
                     occ = rep(1:nocc, each = ndep),
                     s = rep(deployments$start, nocc),
                     e = rep(deployments$end, nocc),
                     c1 = rep(head(cuts, -1), each=ndep),
                     c2 = rep(tail(cuts, -1), each=ndep),
                     i = rep(intervals, each=ndep),
                     z = 0) %>%
    dplyr::mutate(e_c1 = as.numeric(difftime(e, c1, units="day")),
                  c2_s = as.numeric(difftime(c2, s, units="day")),
                  e_s = as.numeric(difftime(e, s, units="day")),
                  sums = (c1<s) + (c2<=s) + (c1<e) + (c2<=e),
                  sel = case_when(sums==1 ~ "e_c1",
                                  sums==2 & c2>e ~ "e_s",
                                  sums==2 & c2<=e ~ "i",
                                  sums==3 ~ "c2_s",
                                  .default = "z")) %>%
    rowwise() %>%
    mutate(eff = get(sel)) %>%
    select(loc, occ, eff) %>%
    dplyr::group_by(loc, occ) %>%
    dplyr::summarise(eff = sum(eff), .groups = "drop") %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = occ, values_from = eff) %>%
    tibble::column_to_rownames("loc") %>%
    as.matrix()
  list(effort = emat, cuts = cuts)
}

#' Make a detection matrix from dataframes
#' 
#' INPUT
#' deployments: dataframe of deployment data with columns as for make_emat plus:
#'    deploymentID: deployment identifier matched with observations$deploymentID
#' observations: dataframe of observation data with columns:
#'    deploymentID: deployment identifiers
#'    timestamp: POSIX date-times of observation occurence
#' cuts: vector sequence of POSIX values (see make_emat)
#' trim: logical, whether incomplete cells (effort<interval) should be set to NA
#' interval / start_hour: passed to make_emat then make_cutSeq if effort is NULL
#' 
#' OUTPUT
#' A sites x occasions detection / non-detection matrix
#' 
#' DETAILS
#' Cuts must span deployment start/end times but not observations times - any
#' observations outside cuts range will simply be ignored without a warning.
make_dmat <- function(deployments, observations, 
                      cuts = NULL,
                      trim=FALSE, 
                      interval=7, 
                      start_hour=0){
  
  effort <- make_emat(deployments, 
                      cuts=cuts, 
                      interval=interval, 
                      start_hour=start_hour)
  if("locationName" %in% names(observations))
    observations <- dplyr::select(observations, -locationName)
  observations <- deployments %>%
    dplyr::select(deploymentID, locationName) %>%
    dplyr::right_join(observations, by=join_by(deploymentID))
  locs <- sort(unique(deployments$locationName))
  nobs <- nrow(observations)
  nloc <- length(locs)
  nocc <- length(effort$cuts) - 1
  ijk <- expand.grid(loc=1:nloc, occ=1:nocc, obs=1:nobs)
  deploc <- locs[ijk$loc]
  obsloc <- observations$locationName[ijk$obs]
  ts <- observations$timestamp[ijk$obs]
  cut1 <- effort$cuts[ijk$occ]
  cut2 <- effort$cuts[ijk$occ+1]
  isin <- ts >= cut1 & ts < cut2 & obsloc == deploc %>%
    array(c(nloc, nocc, nobs))
  mat <- apply(isin, 1:2, any) %>%
    {+.}
  emult <- if(trim) ifelse(effort$effort<interval, NA, 1) else 
    ifelse(effort$effort==0, NA, 1)
  mat * emult
}

# Make a detection matrix from a camtrapDP datapackage
# 
# INPUT
# pkg: a camtrapDP-like list of camera trap data containing
#    package$data$deployments and package$data$observations, dataframes with
#    required columns as for make_dmat, plus scientificName required in observations
# species: a character vector giving one or more species to create matrices for
# trim / interval / start_hour: arguments passed to make_dmat
# 
# OUTPUT
# A list with elements:
#    effort: a sites x occasions effort matrix
#    cuts: POSIX occasion date-time cutpoints
#    matrix: named a list of detection-nondetection matrices, one for each 
#      species named in the species argument
#      
make_detection_matrix <- function(pkg,
                                  species,
                                  trim=FALSE,
                                  interval=7,
                                  start_hour=0){
  
  obsReq <- c("deploymentID", "scientificName", "timestamp")
  depReq <- c("deploymentID", "locationName", "start", "end")
  fieldsOK <- all(obsReq %in% names(pkg$data$observations),
                  depReq %in% names(pkg$data$deployments))
  if(!fieldsOK) 
    stop("Can't find the necessary data: 
         obsdat must contain columns named timestamp and locationName; 
         depdat must contain columns named start, end and locationName")
  
  if(!all(species %in% pkg$data$observations$scientificName))
    stop("Can't find any observations for that/those species")
  
  depdat <- pkg$data$deployments %>%
    dplyr::mutate(deploymentID = as.character(deploymentID),
                  locationName = as.character(locationName))
  obsdat <- pkg$data$observations %>%
    dplyr::mutate(deploymentID = as.character(deploymentID))
  
  missingDeps <- unique(obsdat$deploymentID[!obsdat$deploymentID %in% depdat$deploymentID])
  if(length(missingDeps)>0)
    stop(paste("These deploymentID values in obsdat are missing from depdat:",
               paste(missingDeps, collapse = " ")))
  
  checkdat <- dplyr::left_join(obsdat, 
                               dplyr::select(depdat, deploymentID, start, end),
                               by="deploymentID")
  bad <- with(checkdat, timestamp<start | timestamp>end)
  if(any(bad)){
    message("Error: some observations occur outside their deployment time: 
            returning problematic observations")
    return(checkdat[bad, ])
  } else{
    effort <- make_emat(depdat)
    dmats <- lapply(species, function(sp) 
      make_dmat(depdat, 
                subset(obsdat, scientificName==sp),
                trim = trim, 
                interval = interval, 
                start_hour=start_hour))
    names(dmats) <- species
    return(c(effort, matrix=list(dmats)))
  }
  }

# Result for Tapir

result_tapir <- make_detection_matrix(pkg, species = "Tapirus terrestris")

matriz_tapir <- result_tapir$matrix[["Tapirus terrestris"]]  

# View table

View(matriz_tapir)  


#
#
#
# Result for Puma

result_puma <- make_detection_matrix(pkg, species = "Puma concolor")

matriz_puma <- result_puma$matrix[["Puma concolor"]]  

View(matriz_puma)  


#
#
#
# Result for Agouti Dasyprocta

result <- make_detection_matrix(pkg, species = "Dasyprocta punctata")

matriz_dasyprocta <- result$matrix[["Dasyprocta punctata"]]  

View(matriz_dasyprocta)  


# See how are the scientific names on the table
unique(pkg$data$observations$scientificName)
table(pkg$data$observations$scientificName)

# Install unmarked package

install.packages("unmarked")
library(unmarked)

# To verify if the matriz is in the correct format (sites x occasions)
dim(matriz_tapir) 

# Create unmarked object
umf_tapir <- unmarkedFramePCount(y = matriz_tapir)

# Fit null N-mixture model (constant abundance and detection) (no covariates yet)
modelo_rn <- pcount(~1 ~1, data = umf_tapir)
summary(modelo_rn)

# The average abundance is on the log scale (log(λ) = 0.355), so convert to natural scale:
exp(0.355) # the result is individuals per site


# Confirm parametres

ls()

# View table

View(pkg$data$deployments)

# Confirm names

names(pkg$data$deployments)

# Create site covariate

library(dplyr)

site_covs <- pkg$data$deployments %>%
  distinct(locationName, .keep_all = TRUE) %>%  
  select(
    deploymentID = locationName,   
    distance,                
    transect,                
    zone,                    
    phase                    
  )

# Criar a coluna season corretamente como fator
site_covs$season <- factor(
  ifelse(site_covs$phase == "1", "Dry", "Beggining of rain"),
  levels = c("Dry", "Beggining of rain")
)

View(site_covs)


## Até aqui está certo


# Create object at unmarked

library(unmarked)

umf_tapir <- unmarkedFramePCount(
  y = matriz_tapir,
  siteCovs = site_covs
)


# Simple Model (Null model)

modelo_simples <- pcount(~1 ~1, data = umf_tapir)
summary(modelo_simples)


# Include distance covariate

modelo_dist <- pcount(~1 ~ distance, data = umf_tapir)
summary(modelo_dist)


# Include zone covariate

modelo_zone <- pcount(~1 ~ zone, data = umf_tapir)
summary(modelo_zone)


# Graph distance from road for tapir
install.packages("ggplot2")
library(ggplot2)


# Correcting order of the factors (distance)

site_covs$distance <- factor(site_covs$distance, levels = c("10m", "100m", "250m", "500m", "1000m"))


# Create a new data frame with distance levels
new_data <- data.frame(distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                                         levels = c("10m", "100m", "250m", "500m", "1000m")))


# Predict abundance estimates from the model
pred <- predict(modelo_dist, type = "state", newdata = new_data)

# Combine distance levels with prediction results
plot_data <- cbind(new_data, pred)


# Plot estimated abundance by distance

ggplot(plot_data, aes(x = distance, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(
    title = "Estimated Abundance of Tapirs at Different Distances from the Road",
    x = "Distance from road (m)",
    y = "Estimated Abundance (λ)"
  ) +
  theme_minimal()



#### Plot graph per zone

# Create a new data frame with zone levels
new_data_zone <- data.frame(zone = factor(c("INPA", "TIWA"),
                                          levels = c("INPA", "TIWA")))

# Predict abundance from modelo_zone
pred_zone <- predict(modelo_zone, type = "state", newdata = new_data_zone)

# Combine with prediction results
plot_data_zone <- cbind(new_data_zone, pred_zone)

# Plot estimated abundance by zone
ggplot(plot_data_zone, aes(x = zone, y = Predicted)) +
  geom_point(size = 3, color = "sienna4") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(
    title = "Estimated Abundance of Tapir by Zone",
    x = "Zone",
    y = "Estimated Abundance (λ)"
  ) +
  theme_minimal()

modelo_zone <- pcount(~1 ~ zone, data = umf)
summary(modelo_zone)




####### Combine model distance + zone

modelo_comb <- pcount(~1 ~ distance + zone, data = umf_tapir)
summary(modelo_comb)

# Put combined model in a graph
new_data <- expand.grid(
  distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                    levels = c("10m", "100m", "250m", "500m", "1000m")),
  zone = factor(c("INPA", "TIWA"),
                levels = c("INPA", "TIWA"))
)

pred_comb <- predict(modelo_comb, type = "state", newdata = new_data)
plot_data_comb <- cbind(new_data, pred_comb)
library(ggplot2)

ggplot(plot_data_comb, aes(x = distance, y = Predicted, color = zone, group = zone)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  labs(
    title = "Estimated abundance by distance from road and zone",
    x = "Distance from road (m)",
    y = "Estimated abundance (λ)",
    color = "Zone"
  ) +
  theme_minimal()



#### Coeficients

# Estimativas do modelo zone para a anta (log-scale)
intercept_log <- 0.297  # INPA (referência)
zone_effect_log <- 0.406  # Efeito de estar na TIWA

# Converter para escala real (lambda: abundância esperada)
lambda_inpa <- exp(intercept_log)
lambda_tiwa <- exp(intercept_log + zone_effect_log)

# Calcular a diferença percentual
percent_diff <- (lambda_tiwa - lambda_inpa) / lambda_inpa * 100

# Mostrar o resultado
cat("Estimated abundance INPA:", round(lambda_inpa, 3), "\n")
cat("Estimated abundance TIWA:", round(lambda_tiwa, 3), "\n")
cat("Percentual difference between TIWA and INPA:", round(percent_diff, 1), "%\n")



#### Include season as covariate

modelo_season <- pcount(~1 ~ season, data = umf_tapir)
summary(modelo_season)

### Graph for season

# Prever abundância para cada estação
new_data <- data.frame(season = c("Dry", "Beggining of rain"))

# Previsão com intervalo de confiança
pred_season <- predict(modelo_season, type = "state", newdata = new_data)

# Juntar previsões com os dados
plot_data <- cbind(new_data, pred_season)

# Carregar o pacote ggplot2
library(ggplot2)

ggplot(plot_data, aes(x = season, y = Predicted)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(
    title = "Estimated Abundance of Tapir by Season",
    x = "Season",
    y = "Estimated Abundance (λ)"
  ) +
  theme_minimal()




########
# Doing the same thing but for Dasyprocta punctata

result_cotia <- make_detection_matrix(pkg, species = "Dasyprocta punctata")
matriz_cotia <- result$matrix[["Dasyprocta punctata"]]
View(matriz_cotia)

# Create object at unmarked

library(unmarked)

umf_cotia <- unmarkedFramePCount(
  y = matriz_cotia,
  siteCovs = site_covs
)

# Correcting order of the factors

site_covs$distance <- factor(site_covs$distance, levels = c("10m", "100m", "250m", "500m", "1000m"))

# Model considering distances

modelo_dist_cotia <- pcount(~1 ~ distance, data = umf_cotia)
summary(modelo_dist_cotia)


# Criar tabela com os níveis de distância
new_data <- data.frame(distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                                         levels = c("10m", "100m", "250m", "500m", "1000m")))


# Prever abundância estimada pelo modelo
pred <- predict(modelo_dist, type = "state", newdata = new_data)

# Juntar dados para plotagem
plot_data <- cbind(new_data, pred)

library(ggplot2)

new_data <- data.frame(distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                                         levels = c("10m", "100m", "250m", "500m", "1000m")))

pred_cotia <- predict(modelo_dist_cotia, type = "state", newdata = new_data)

plot_data_cotia <- cbind(new_data, pred_cotia)

ggplot(plot_data_cotia, aes(x = distance, y = Predicted)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(
    title = "Estimated abundance of Agouti at different distances",
    x = "Distance from road (m)",
    y = "Estimated abundance (λ)"
  ) +
  theme_minimal()


## Agouti per zone

modelo_zone_cotia <- pcount(~1 ~ zone, data = umf_cotia)
summary(modelo_zone_cotia)

# Criar novo data.frame com os níveis da variável "zone"
new_data_zone_cotia <- data.frame(zone = factor(c("INPA", "TIWA"),
                                                levels = c("INPA", "TIWA")))

# Prever abundância estimada com base no modelo
pred_zone_cotia <- predict(modelo_zone_cotia, type = "state", newdata = new_data_zone_cotia)

# Combinar os dados com as previsões
plot_data_zone_cotia <- cbind(new_data_zone_cotia, pred_zone_cotia)

# Plotar gráfico
library(ggplot2)

ggplot(plot_data_zone_cotia, aes(x = zone, y = Predicted)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(
    title = "Estimated Abundance of Agouti (Dasyprocta punctata) by Zone",
    x = "Zone",
    y = "Estimated Abundance (λ)"
  ) +
  theme_minimal()


####### Combine model distance + zone for Agouti

modelo_comb_cotia <- pcount(~1 ~ distance + zone, data = umf_cotia)
summary(modelo_comb_cotia)

# Put combined model in a graph
new_data_cotia <- expand.grid(
  distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                    levels = c("10m", "100m", "250m", "500m", "1000m")),
  zone = factor(c("INPA", "TIWA"),
                levels = c("INPA", "TIWA"))
)

pred_comb_cotia <- predict(modelo_comb_cotia, type = "state", newdata = new_data_cotia)
plot_data_comb_cotia <- cbind(new_data_cotia, pred_comb_cotia)
library(ggplot2)

ggplot(plot_data_comb_cotia, aes(x = distance, y = Predicted, color = zone, group = zone)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  labs(
    title = "Estimated abundance by distance from road and zone (Agouti)",
    x = "Distance from road (m)",
    y = "Estimated abundance (λ)",
    color = "Zone"
  ) +
  theme_minimal()


#### Coeficients

# Estimativas do modelo zone para a cotia (log-scale)
intercept_log <- -0.2753  # INPA (referência)
zone_effect_log <- -0.0355  # Efeito de estar na TIWA

# Converter para escala real (lambda: abundância esperada)
lambda_inpa <- exp(intercept_log)
lambda_tiwa <- exp(intercept_log + zone_effect_log)

# Calcular a diferença percentual
percent_diff <- (lambda_tiwa - lambda_inpa) / lambda_inpa * 100

# Mostrar o resultado
cat("Estimated abundance INPA:", round(lambda_inpa, 3), "\n")
cat("Estimated abundance TIWA:", round(lambda_tiwa, 3), "\n")
cat("Percentual difference between TIWA and INPA:", round(percent_diff, 1), "%\n")




########
# Doing the same thing but for queixada 

result_tayassu <- make_detection_matrix(pkg, species = "Tayassu pecari")
matriz_tayassu <- result_tayassu$matrix[["Tayassu pecari"]]
View(matriz_tayassu)

# Create object at unmarked

library(unmarked)

umf_tayassu <- unmarkedFramePCount(
  y = matriz_tayassu,
  siteCovs = site_covs
)

# Correcting order of the factors

site_covs$distance <- factor(site_covs$distance, levels = c("10m", "100m", "250m", "500m", "1000m"))

# Model considering distances

modelo_dist_tayassu <- pcount(~1 ~ distance, data = umf_tayassu)
summary(modelo_dist_tayassu)


# Criar tabela com os níveis de distância
new_data <- data.frame(distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                                         levels = c("10m", "100m", "250m", "500m", "1000m")))


# Prever abundância estimada pelo modelo
pred <- predict(modelo_dist, type = "state", newdata = new_data)

# Juntar dados para plotagem
plot_data <- cbind(new_data, pred)

library(ggplot2)

new_data <- data.frame(distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                                         levels = c("10m", "100m", "250m", "500m", "1000m")))

pred_tayassu <- predict(modelo_dist_tayassu, type = "state", newdata = new_data)

plot_data_tayassu <- cbind(new_data, pred_tayassu)

ggplot(plot_data_tayassu, aes(x = distance, y = Predicted)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(
    title = "Estimated abundance of Tayassu at different distances",
    x = "Distance from road (m)",
    y = "Estimated abundance (λ)"
  ) +
  theme_minimal()


## Tayassu per zone

modelo_zone_tayassu <- pcount(~1 ~ zone, data = umf_tayassu)
summary(modelo_zone_tayassu)

# Criar novo data.frame com os níveis da variável "zone"
new_data_zone_tayassu <- data.frame(zone = factor(c("INPA", "TIWA"),
                                                levels = c("INPA", "TIWA")))

# Prever abundância estimada com base no modelo
pred_zone_tayassu <- predict(modelo_zone_tayassu, type = "state", newdata = new_data_zone_tayassu)

# Combinar os dados com as previsões
plot_data_zone_tayassu <- cbind(new_data_zone_tayassu, pred_zone_tayassu)

# Plotar gráfico
library(ggplot2)

ggplot(plot_data_zone_tayassu, aes(x = zone, y = Predicted)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(
    title = "Estimated Abundance of Tayassu by Zone",
    x = "Zone",
    y = "Estimated Abundance (λ)"
  ) +
  theme_minimal()


####### Combine model distance + zone for Tayassu

modelo_comb_tayassu <- pcount(~1 ~ distance + zone, data = umf_tayassu)
summary(modelo_comb_tayassu)

# Put combined model in a graph
new_data_tayassu <- expand.grid(
  distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                    levels = c("10m", "100m", "250m", "500m", "1000m")),
  zone = factor(c("INPA", "TIWA"),
                levels = c("INPA", "TIWA"))
)

pred_comb_tayassu <- predict(modelo_comb_tayassu, type = "state", newdata = new_data_tayassu)
plot_data_comb_tayassu <- cbind(new_data_tayassu, pred_comb_tayassu)
library(ggplot2)

ggplot(plot_data_comb_tayassu, aes(x = distance, y = Predicted, color = zone, group = zone)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  labs(
    title = "Estimated abundance by distance from road and zone (Tayassu)",
    x = "Distance from road (m)",
    y = "Estimated abundance (λ)",
    color = "Zone"
  ) +
  theme_minimal()


#### Coeficients

# Estimativas do modelo zone para a tayassu (log-scale)
intercept_log <- -1.602  # INPA (referência)
zone_effect_log <- 0.465  # Efeito de estar na TIWA

# Converter para escala real (lambda: abundância esperada)
lambda_inpa <- exp(intercept_log)
lambda_tiwa <- exp(intercept_log + zone_effect_log)

# Calcular a diferença percentual
percent_diff <- (lambda_tiwa - lambda_inpa) / lambda_inpa * 100

# Mostrar o resultado
cat("Estimated abundance INPA:", round(lambda_inpa, 3), "\n")
cat("Estimated abundance TIWA:", round(lambda_tiwa, 3), "\n")
cat("Percentual difference between TIWA and INPA:", round(percent_diff, 1), "%\n")



############### Bayesian model

install.packages("ubms")
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()  # roda só uma vez

library(ubms)
library(cmdstanr)
library(ggplot2)

# Correct order of factors in covariantes site
site_covs$distance <- factor(site_covs$distance,
                             levels = c("10m", "100m", "250m", "500m", "1000m"))

# Create unmarked object
umf_cotia <- unmarkedFramePCount(
  y = matriz_cotia,
  siteCovs = site_covs
)

# Bayesian model (stan)
modelo_bayes_cotia <- stan_pcount(~1 ~ distance,
                                  data = umf_cotia,
                                  chains = 4,
                                  iter = 2000,
                                  seed = 123)

# Sumary 
summary(modelo_bayes_cotia, submodel = "state")

# Create new table with distances
new_data <- data.frame(distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                                         levels = c("10m", "100m", "250m", "500m", "1000m")))

# Predict abundance
pred_cotia_bayes <- predict(modelo_bayes_cotia, submodel = "state", newdata = new_data)

# Verify columns names
names(pred_cotia_bayes)

# Plot data
plot_data_cotia_bayes <- cbind(new_data, pred_cotia_bayes)

# Making the graph
library(ggplot2)
ggplot(plot_data_cotia_bayes, aes(x = distance, y = Predicted)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1) +
  labs(
    title = "Estimated Abundance of Agouti (Bayesian model)",
    x = "Distance from road (m)",
    y = "Estimated Abundance (λ)"
  ) +
  theme_minimal()


View(site_covs)

