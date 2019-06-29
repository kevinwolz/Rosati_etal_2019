### MAIN
### Programmer: Kevin Wolz (kevin@savannainstitute.org)
### Originally Created: 15 Apr 2017
### Last Updated: 28 Jun 2019

## DIRECTORY SETUP
modelPath      <- "./model/"
rawDataPath    <- "./raw_data/"
outputDataPath <- "./output/data/"
dir.create(outputDataPath, showWarnings = FALSE, recursive = TRUE)

## REQUIRED LIBRARIES
library(tidyverse)
library(lubridate)
library(foreach)
library(snow)
library(doSNOW)

## SOURCE SORTIE MODEL FUNCTIONS
source(paste0(modelPath, "solar_geometry.R"))
source(paste0(modelPath, "tree_geometry.R"))
source(paste0(modelPath, "light_model.R"))

## READ IN RAW INPUT DATA
scenarios     <- read_csv(paste0(rawDataPath, "scenarios.csv"), col_types = cols())
orchard.trees <- read_csv(paste0(rawDataPath, "trees.csv"),     col_types = cols())
sensors       <- read_csv(paste0(rawDataPath, "sensors.csv"),   col_types = cols())

## SPECIFY SCENARIOS TO RUN
SCENARIOS.TO.RUN <- scenarios$scenario.id

## RUN SCENARIOS IN PARALLEL
## Running the points in serial occurs 10% slower than for running the points in parallel,
## but the scenarios can be run in parallel, which quadruples the speed at that level.
scenario_call <- function(scenario) {

  ## PULL OUT DATA FOR A GIVEN SCENARIO
  sub.scenarios <- dplyr::filter(scenarios,     scenario.id == scenario)
  sub.trees     <- dplyr::filter(orchard.trees, orchard.id  == sub.scenarios$orchard.id)
  sub.sensors   <- dplyr::filter(sensors,       orchard.id  == sub.scenarios$orchard.id)
  shape <- sub.scenarios$canopy.shape

  ## CREATE TIME SERIES AT WHICH TO RUN MODEL
  ## CONVERT TIME FROM CT (Local Time) to Solar Time
  ## (the Central Time tz takes into account DST, so no need to subtract an hour to correct for DST)
  start.time       <- lubridate::mdy_hm(sub.scenarios$time.start, tz = "US/Central")
  end.time         <- lubridate::mdy_hm(sub.scenarios$time.end,   tz = "US/Central")
  solar.start.time <- solaR::local2Solar(start.time, lon = sub.scenarios$long)
  solar.end.time   <- solaR::local2Solar(end.time,   lon = sub.scenarios$long)
  times            <- as.POSIXct(round(seq(solar.start.time, solar.end.time, "1 min"), "mins"))

  ## CREATE TREE FILE FOR MODEL
  ## THIS INCLUDES TWO RINGS OF TREES AROUND THE FOUR CORE TREES
  trees <- data.frame(tree.id = sub.trees$tree.id,
                      X       = sub.trees$X,
                      Y       = sub.trees$Y,
                      TRANS   = ifelse(shape == "C", sub.trees$trans.cy, sub.trees$trans.el),
                      SHAPE   = shape,
                      # Cylinder-only parameters
                      TOP      = ifelse(shape == "C", sub.trees$height.cy, NA),
                      BOTTOM   = ifelse(shape == "C", sub.trees$bole.height.cy, NA),
                      RADIUS   = ifelse(shape == "C", ((sub.trees$diam.x + sub.trees$diam.y) / 2), NA),
                      # Ellipsoid-only parameters
                      CENTER   = ifelse(shape == "E", ((sub.trees$height.el + sub.trees$bole.height.el) / 2), NA),
                      XRAD     = ifelse(shape == "E", sub.trees$diam.x, NA),
                      YRAD     = ifelse(shape == "E", sub.trees$diam.y, NA),
                      VRAD     = ifelse(shape == "E", (sub.trees$height.el - sub.trees$bole.height.el), NA))

  ## CREATE POINTS FILE FOR MODEL
  points <- expand.grid(sensor.id = sub.sensors$sensor.id, times = times) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(sub.sensors, by = "sensor.id") %>%
    dplyr::mutate(DATE = substr(times, 1,  10)) %>%
    dplyr::mutate(TIME = substr(times, 12, 16)) %>%
    dplyr::select(-times, -orchard.id) %>%
    dplyr::arrange(DATE, TIME, sensor.id) %>%
    dplyr::mutate(GLI = NA) # Add new column for output

  ## OTHER SCENARIO PARAMETERS
  latitude     <- sub.scenarios$lat
  azi_of_north <- sub.scenarios$y.orient
  overcast     <- sub.scenarios$overcast

  result <- purrr::map_df(1:nrow(points),
                          light_model,
                          trees        = trees,
                          points       = points,
                          latitude     = latitude,
                          azi_of_north = azi_of_north,
                          overcast     = overcast)

  ## WRITE OUT MODELED DATA FOR SCENARIO
  points <- result %>%
    dplyr::arrange(DATE, TIME, sensor.id) %>%
    dplyr::mutate(scenario.id = scenario) %>%
    dplyr::left_join(scenarios, by = "scenario.id") %>%
    dplyr::select(-time.start, -time.end)

  readr::write_csv(points, paste0(outputDataPath, "Model_Output_", scenario, ".csv"))

  return(data.frame(done = scenario))
}

## START TIME & PRINT STATUS
TimerSplitStart <- TimerStart <-  proc.time()[3] # Start timer
print("Starting scenarios")

myCluster <- makeCluster(spec = 4, type = "SOCK")
clusterExport(myCluster, c("scenarios", "orchard.trees", "sensors", "outputDataPath"))
registerDoSNOW(myCluster)

result <- foreach(i         = SCENARIOS.TO.RUN,
                  .combine  = "rbind",
                  .inorder  = FALSE,
                  .packages = c("solaR", "dplyr")) %dopar% scenario_call(i)

stopCluster(myCluster)
registerDoSEQ()

## CHECK TIMER & PRINT STATUS
TimerElapsed <- proc.time()[3] - TimerStart
print(paste("Scenarios completed in", round(TimerElapsed/60, 1), "minutes"))
