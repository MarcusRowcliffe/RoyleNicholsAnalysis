# install packages
install.packages(c("devtools", "tidyverse", "unmarked"))
devtools::install_github("inbo/camtraptor")
devtools::install_github("MarcusRowcliffe/camtrapDensity")

# load packages
library(camtrapDensity)
library(dplyr)
library(purrr)

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
