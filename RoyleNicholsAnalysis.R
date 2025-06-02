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

result <- make_detection_matrix(pkg, species = "Tapirus terrestris")

matriz_tapir <- resultado$matrix[["Tapirus terrestris"]]  

# View table

View(matriz_tapir)  


#
#
#
# Result for Puma

result <- make_detection_matrix(pkg, species = "Puma concolor")

matriz_puma <- result$matrix[["Puma concolor"]]  

# View table

View(matriz_puma)  


#
#
#
# Result for Agouti Dasyprocta

result <- make_detection_matrix(pkg, species = "Dasyprocta punctata")

matriz_dasyprocta <- result$matrix[["Dasyprocta punctata"]]  

# View table

View(matriz_dasyprocta)  



unique(pkg$data$observations$scientificName)

matriz_tapir <- result$matrix[["Tapirus terrestris"]]


install.packages("unmarked")
library(unmarked)

# Verifique se a matriz está no formato correto (sites x occasions)
dim(matriz_tapir)  # Deve ser [n_sites x n_occasions]

# Criar o objeto unmarked
umf <- unmarkedFramePCount(y = matriz_tapir)
modelo_rn <- pcount(~1 ~1, data = umf)
summary(modelo_rn)

# Create covariates table

site_covs <- pkg$data$deployments %>%
  dplyr::distinct(locationName, .keep_all = TRUE) %>%  # um por site
  dplyr::select(locationName, distance, inside_reserve = zone) %>%  # escolha as covariáveis
  dplyr::mutate(
    distance = as.numeric(distance),                  # numérica
    inside_reserve = as.factor(inside_reserve)        # categórica (fator)
  ) %>%
  tibble::column_to_rownames("locationName")


# Confirm parametres

ls()

# View table

View(pkg$data$deployments)

# Confirm names

names(pkg$data$deployments)

# Create site covariate

library(dplyr)

site_covs <- pkg$data$deployments %>%
  distinct(locationName, .keep_all = TRUE) %>%  # Garante um por site
  select(
    deploymentID = locationName,   # coluna 1: ID do site
    distance,                # coluna 2: distância da estrada
    transect,                # coluna 3: transecto
    zone,                    # coluna 4: dentro/fora da TI
    phase                    # coluna 5: campanha
  )

View(site_covs)


## Até aqui está certo



# Create again Tapir Matriz

result <- make_detection_matrix(pkg, species = "Tapirus terrestris")
matriz_tapir <- result$matrix[["Tapirus terrestris"]]
View(matriz_tapir)


# Create object at unmarked

library(unmarked)

umf <- unmarkedFramePCount(
  y = matriz_tapir,
  siteCovs = site_covs
)


# Simple Model

modelo_simples <- pcount(~1 ~1, data = umf)
summary(modelo_simples)


# Include covariates

modelo_rn <- pcount(~1 ~ distance + transect + zone + phase, data = umf)
summary(modelo_rn)


# Graph distance from road
install.packages("ggplot2")
library(ggplot2)


# Correcting order of the factors

site_covs$distance <- factor(site_covs$distance, levels = c("10m", "100m", "250m", "500m", "1000m"))

# Model considering distances

modelo_dist <- pcount(~1 ~ distance, data = umf)
summary(modelo_dist)


# Criar tabela com os níveis de distância
new_data <- data.frame(distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                                         levels = c("10m", "100m", "250m", "500m", "1000m")))


# Prever abundância estimada pelo modelo
pred <- predict(modelo_dist, type = "state", newdata = new_data)

# Juntar dados para plotagem
plot_data <- cbind(new_data, pred)

library(ggplot2)

ggplot(plot_data, aes(x = distance, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(
    title = "Abundância estimada de anta por distância da estrada",
    x = "Distância da estrada (m)",
    y = "Abundância estimada (λ)"
  ) +
  theme_minimal()


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
    title = "Abundância estimada de cotia por distância da estrada",
    x = "Distância da estrada (m)",
    y = "Abundância estimada (λ)"
  ) +
  theme_minimal()

# Same thing for Tayassu pecari
# Create detection matriz

result_tayassu <- make_detection_matrix(pkg, species = "Tayassu pecari")
matriz_tayassu <- result_tayassu$matrix[["Tayassu pecari"]]
View(matriz_tayassu)

# 1. Criar a matriz de detecção
result_queixada <- make_detection_matrix(pkg, species = "Tayassu pecari")
matriz_queixada <- result_queixada$matrix[["Tayassu pecari"]]

# 2. Alinhar os sites da matriz com os da tabela de covariáveis
site_covs_queixada <- site_covs[rownames(matriz_queixada), ]

# 3. Criar o objeto unmarked
umf_queixada <- unmarkedFramePCount(
  y = matriz_queixada,
  siteCovs = site_covs_queixada
)

# 4. Garantir que distance esteja como fator na ordem certa
site_covs_queixada$distance <- factor(site_covs_queixada$distance,
                                      levels = c("10m", "100m", "250m", "500m", "1000m"))

# 5. Rodar o modelo considerando a distância
modelo_dist_queixada <- pcount(~1 ~ distance, data = umf_queixada, K = 20)
summary(modelo_dist_queixada)

# 6. Prever abundância estimada para cada distância
new_data <- data.frame(distance = factor(c("10m", "100m", "250m", "500m", "1000m"),
                                         levels = c("10m", "100m", "250m", "500m", "1000m")))

pred_queixada <- predict(modelo_dist_queixada, type = "state", newdata = new_data)
plot_data_queixada <- cbind(new_data, pred_queixada)

# 7. Plotar o gráfico
library(ggplot2)

ggplot(plot_data_queixada, aes(x = distance, y = Predicted)) +
  geom_point(size = 3, color = "firebrick") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  labs(
    title = "Abundância estimada de queixada (Tayassu pecari) por distância da estrada",
    x = "Distância da estrada (m)",
    y = "Abundância estimada (λ)"
  ) +
  theme_minimal()

# Align with covariates site table

site_covs_tayassu <- site_covs[rownames(matriz_tayassu), ]

# Create object at unmarked

umf_tayassu <- unmarkedFramePCount(
  y = matriz_tayassu,
  siteCovs = site_covs_tayassu
)

# Generate model with distance

modelo_dist_tayassu <- pcount(~1 ~ distance, data = umf_tayassu, K = 20)
summary(modelo_dist_tayassu)

sum(matriz_tayassu, na.rm = TRUE)
