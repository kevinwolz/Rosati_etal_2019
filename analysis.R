### ANALYSIS
### Programmer: Kevin Wolz (kevin@savannainstitute.org)
### Originally Created: 15 Apr 2017
### Last Updated: 28 Jun 2019

## DIRECTORY SETUP
measuredDataPath <- "./raw_data/measured_data/"
outputDataPath   <- "./output/data/modeled_data/"
outputPlotPath   <- "./output/plots/"
dir.create(outputPlotPath, showWarnings = FALSE, recursive = TRUE)

## USER PARAMETERS
REDO.SENSITIVITY <- FALSE
BASE.SIZE <- 20
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## FARQUHAR MODEL
source("Farquhar.R")

#Guo et al 2015
V <- 125.1
J <- 209.75
R <- -1.15

## REQUIRED LIBRARIES
library(tidyverse)
library(DeLuciatoR)
library(lubridate)
library(stringr)
library(solaR)

## NAMES
orchard.names <- c("Young A", "Young B", "Mature A", "Mature B")
orchard.shapes <- 1:6

##### DATA PREP #####
## READ IN MODELED DATA & CREATE TIMESTAMP
modeled <- list.files(path = outputDataPath, full.names = T) %>%
  purrr::map(read.csv, header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::bind_rows() %>%
  dplyr::as_tibble() %>%
  dplyr::filter(orchard.id %in% orchard.names) %>%
  dplyr::select(-trans) %>%
  dplyr::mutate(timestamp = as.POSIXct(strptime(paste(DATE, TIME), format = "%F %R", tz = "UTC")))

modeled$timestamp[which(is.na(modeled$timestamp))] <- as.POSIXct(strptime(paste(modeled$DATE[which(is.na(modeled$timestamp))],
                                                                                modeled$TIME[which(is.na(modeled$timestamp))]),
                                                                          format = "%m/%d/%y %R", tz = "UTC"))

modeled <- modeled %>%
  dplyr::mutate(DATE = substr(as.character(timestamp), 1,  10)) %>%
  dplyr::mutate(TIME = substr(as.character(timestamp), 12, 16)) %>%
  dplyr::rename(modeled.GLI = GLI) %>%
  dplyr::select(-scenario.id, -X, -Y, -Z, -lat, -long, -y.orient, -daily.PARi) %>%
  dplyr::arrange(orchard.id, canopy.shape, sensor.id, DATE, TIME)

## READ IN MEAURED DATA & CREATE TIMESTAMP (and convert from Local Time to Solar Time)
measured <- list.files(path = measuredDataPath, full.names = T) %>%
  purrr::map(readr::read_csv, col_types = readr::cols()) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(orchard.id %in% orchard.names) %>%
  dplyr::rename(incident.PAR = incident) %>%
  tidyr::gather(-orchard.id, -timestamp, -incident.PAR, key = "sensor.id", value = "measured.PAR") %>%
  dplyr::mutate(timestamp = lubridate::mdy_hm(timestamp, tz = "US/Central")) %>%
  dplyr::mutate(timestamp = solaR::local2Solar(timestamp, lon = -92.7656)) %>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp)) %>%
  dplyr::mutate(timestamp = lubridate::round_date(timestamp, "minutes")) %>%
  dplyr::mutate(sensor.id = stringr::str_remove(sensor.id, "sensor\\.")) %>%
  dplyr::mutate(sensor.id = as.numeric(sensor.id)) %>%
  dplyr::mutate(day       = lubridate::day(timestamp)) %>%
  dplyr::filter(!(sensor.id %in% 19:24 & orchard.id %in% c("Young A", "Young B"))) ## ERRONEOUS DATA


## INTERPOLATE MISSING MEASURED PAR DATA
measured$measured.PAR[measured$measured.PAR <  0] <- 0  # REMOVE NEGATIVE MEASURED PAR
measured$measured.PAR[measured$measured.PAR == 0] <- NA

for(o in unique(measured$orchard.id)){
  o.measured <- subset(measured, orchard.id == o)
  for(s in unique(o.measured$sensor.id)){
    s.measured <- subset(o.measured, sensor.id == s)
    for(d in unique(s.measured$day)){
      X <- as.numeric(measured$timestamp[which(measured$orchard.id == o & measured$sensor.id == s & measured$day == d)])
      Ym <- measured$measured.PAR[which(measured$orchard.id == o & measured$sensor.id == s & measured$day == d)]
      Yi <- measured$incident.PAR[which(measured$orchard.id == o & measured$sensor.id == s & measured$day == d)]
      if(any(is.na(Ym)) & !all(is.na(Ym))){
        yout <- approx(x = X, y = Ym, xout = X, rule = 2)$y
        measured$measured.PAR[which(measured$orchard.id == o & measured$sensor.id == s & measured$day == d)] <- yout
      }
      if(any(is.na(Yi)) & !all(is.na(Yi))){
        yout <- approx(x = X, y = Yi, xout = X, rule = 2)$y
        measured$incident.PAR[which(measured$orchard.id == o & measured$sensor.id == s & measured$day == d)] <- yout
      }
    }
  }
}

## COMBINE MEASURED & MODELED DATA
dat <- measured %>%
  dplyr::select(-day) %>%
  dplyr::full_join(modeled, by = c("orchard.id", "timestamp", "sensor.id")) %>%
  dplyr::filter(!is.na(modeled.GLI)) %>%
  dplyr::mutate(modeled.PAR = modeled.GLI * incident.PAR) %>%
  dplyr::mutate(hour        = as.numeric(substr(TIME, 1, 2))) %>%
  dplyr::mutate(minute      = as.numeric(substr(TIME, 4, 5))) %>%
  dplyr::mutate(halfhour    = as.numeric(minute < 30)) %>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp)) %>%
  dplyr::mutate(canopy.shape = factor(canopy.shape, levels = c("E", "C"), labels = c("Ellipsoid", "Cylinder"))) %>%
  dplyr::mutate(orchard.id = factor(orchard.id,
                                    levels = orchard.names,
                                    labels = orchard.names)) %>%
  dplyr::as_tibble()

##### DATA CALCULATIONS & CONVERSIONS #####
calc_A <- function(d, V, J) {
  d %>%
    dplyr::mutate(modeled.A  = Farquhar_VJ(modeled.PAR,  V, J, R)) %>%
    dplyr::mutate(measured.A = Farquhar_VJ(measured.PAR, V, J, R))
}

## CREATE DIFFERENTE TIME FRAMES
min.dat <- dat %>%
  dplyr::filter(hour %in% 6:16) %>%
  dplyr::select(orchard.id, sensor.id, DATE, hour, halfhour, TIME, canopy.shape, overcast,
                incident.PAR, modeled.PAR, measured.PAR) %>%
  calc_A(V = V, J = J) ## CALCULATE PHOTOSYNTHESIS

readr::write_csv(min.dat, "./output/data/HARC_Chestnuts_Processed_Data.csv")

hour.dat <- min.dat %>%
  dplyr::group_by(orchard.id, sensor.id, DATE, hour, canopy.shape, overcast) %>%
  dplyr::summarize(modeled.PAR  = mean(modeled.PAR,  na.rm = TRUE),
                   measured.PAR = mean(measured.PAR, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  calc_A(V = V, J = J) ## CALCULATE PHOTOSYNTHESIS

thirty.dat <- min.dat %>%
  dplyr::group_by(orchard.id, sensor.id, DATE, hour, halfhour, canopy.shape, overcast) %>%
  dplyr::summarize(modeled.PAR  = mean(modeled.PAR,  na.rm = TRUE),
                   measured.PAR = mean(measured.PAR, na.rm = TRUE)) %>%
  calc_A(V = V, J = J) ## CALCULATE PHOTOSYNTHESIS

day.dat <- min.dat %>%
  dplyr::group_by(orchard.id, sensor.id, DATE, canopy.shape, overcast) %>%
  dplyr::summarize(modeled.PAR  = mean(modeled.PAR,  na.rm = TRUE),
                   measured.PAR = mean(measured.PAR, na.rm = TRUE)) %>%
  calc_A(V = V, J = J) ## CALCULATE PHOTOSYNTHESIS

## CONVERT PAR & A TO CORRECT UNITS FOR EACH TIMEFRAME
cols.to.convert <- c("modeled.PAR", "measured.PAR", "modeled.A", "measured.A")

min.dat <- min.dat %>%
  dplyr::mutate_at(cols.to.convert, function(x) x * 60) ## CONVERT PAR & A FROM umol m-2 s-1 TO umol m-2 min-1

hour.dat <- hour.dat %>%
  dplyr::mutate_at(cols.to.convert, function(x) x * 60 * 60) ## CONVERT PAR & A FROM umol m-2 min-1 TO umol m-2 hr-1

thirty.dat <- thirty.dat %>%
  dplyr::mutate_at(cols.to.convert, function(x) x * 60 * 30) ## CONVERT PAR & A FROM umol m-2 min-1 TO umol m-2 half-hour-1

day.dat <- day.dat %>%
  dplyr::mutate_at(cols.to.convert, function(x) x * 60 * 60 *11) #* (11 + 49/60) ## CONVERT PAR & A FROM umol m-2 min-1 TO umol m-2 day-1

## AGGREGATE UP TO DAILY DATA FOR ALL AVERAGING TIMEFRAMES
aggregate <- function(d) {
  d %>%
    dplyr::group_by(orchard.id, sensor.id, DATE, canopy.shape, overcast) %>%
    dplyr::summarize(modeled.PAR  = sum(modeled.PAR,  na.rm = TRUE),
                     measured.PAR = sum(measured.PAR, na.rm = TRUE),
                     modeled.A    = sum(modeled.A,    na.rm = TRUE),
                     measured.A   = sum(measured.A,   na.rm = TRUE)) %>%
    ungroup()
}

convert_units <- function(d) {
  d %>%
    dplyr::mutate(modeled.PAR  = modeled.PAR  / 1e6) %>% ## CONVERT UNITS TO mol  m-2 day-1 FOR PAR
    dplyr::mutate(measured.PAR = measured.PAR / 1e6) %>% ## CONVERT UNITS TO mol  m-2 day-1 FOR PAR
    dplyr::mutate(modeled.A    = modeled.A    / 1e3) %>% ## CONVERT UNITS TO mmol m-2 day-1 FOR PAR
    dplyr::mutate(measured.A   = measured.A   / 1e3)     ## CONVERT UNITS TO mmol m-2 day-1 FOR PAR
}

timeframe.levels <- c("Minutely", "Half-hourly mean", "Hourly mean", "Daily mean")

day.min.dat    <- min.dat    %>% aggregate() %>% convert_units() %>% dplyr::mutate(timeframe = timeframe.levels[1])
day.thirty.dat <- thirty.dat %>% aggregate() %>% convert_units() %>% dplyr::mutate(timeframe = timeframe.levels[2])
day.hour.dat   <- hour.dat   %>% aggregate() %>% convert_units() %>% dplyr::mutate(timeframe = timeframe.levels[3])
day.day.dat    <- day.dat    %>%                 convert_units() %>% dplyr::mutate(timeframe = timeframe.levels[4])


##### CHECK FOR OPTIMAL POROSITY USING AVERAGE DAILY PAR COMPARISON #####
orchard.avg <- day.min.dat %>%
  dplyr::group_by(orchard.id, DATE, canopy.shape, overcast) %>%
  dplyr::summarize(modeled.PAR  = mean(modeled.PAR,  na.rm = TRUE),
                   measured.PAR = mean(measured.PAR, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(DOM = substr(DATE, 9, 10))

max.val <- max(pmax(orchard.avg$measured.PAR, orchard.avg$modeled.PAR))

orchard.avg.together.fig <- ggplot(orchard.avg, aes(y     = modeled.PAR,
                                                    x     = measured.PAR,
                                                    color = canopy.shape,
                                                    shape = orchard.id)) +
  labs(y     = bquote('Modeled daily PAR (mol ' ~m^-2~ ' '~day^-1*~')'),
       x     = bquote('Measured daily PAR (mol ' ~m^-2~ ' '~day^-1*~')'),
       color = NULL,
       shape = NULL) +
  scale_x_continuous(limits = c(0, max.val), sec.axis = sec_axis(~ ., labels=NULL)) +
  scale_y_continuous(limits = c(0, max.val), sec.axis = sec_axis(~ ., labels=NULL)) +
  geom_point(size = 4, stroke = 1) +
  scale_color_manual(values = c("black", "grey50")) +
  scale_shape_manual(values = orchard.shapes) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  guides(color = FALSE) +
  theme_ggEHD(base_size = BASE.SIZE) +
  theme(legend.position = c(0.2, 0.7),
        legend.key.width  = unit(0.25, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.background = element_rect(fill="transparent"))

ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_1_Daily_Orchard_Avg_PAR.jpg"),
              plot     = orchard.avg.together.fig,
              scale    = 1,
              dpi      = 300)


##### COMBINE & COMPARE DIFFERENT AVERAGING APPROACHES - PAR #####
combinded.dat <- dplyr::bind_rows(day.min.dat, day.hour.dat, day.thirty.dat, day.day.dat) %>%
  mutate(timeframe = factor(timeframe, levels = timeframe.levels))

combinded.dat.melted <- combinded.dat %>%
  tidyr::gather(modeled.A, measured.A, key = "type", value = "value") %>%
  dplyr::mutate(type = factor(type, levels = c("modeled.A", "measured.A"),
                                    labels = c("Modeled", "Measured"))) %>%
  dplyr::arrange(dplyr::desc(type)) %>%
  dplyr::mutate(respective.PAR = modeled.PAR)

combinded.dat.melted$respective.PAR[which(combinded.dat.melted$type == "Measured")] <- combinded.dat.melted$measured.PAR[which(combinded.dat.melted$type == "Measured")]

com_fig_respective <- function(dat, c) {

  dat <- dplyr::filter(dat, canopy.shape == c)

  plot.annotation <- data.frame(timeframe      = timeframe.levels,
                                respective.PAR = 0,
                                value          = max(combinded.dat.melted$value, na.rm = TRUE),
                                label          = paste(c("(a)", "(b)", "(c)", "(d)"), timeframe.levels))

  ggplot(dat, aes(x = respective.PAR, y = value)) + #shape = factor(sensor.id)
    labs(y     = bquote('Daily ' ~A[n]~ ' (mmol ' ~CO[2]~ ' '~m^-2~ ' '~day^-1*')'),
         x     = bquote('Daily incident PAR (mol ' ~m^-2~ ' '~day^-1*')'),
         color = "") +
    scale_x_continuous(sec.axis = sec_axis(~ ., labels=NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ ., labels=NULL)) +
    facet_wrap(~timeframe, ncol = 1) +
    geom_point(size = 2, na.rm = T, aes(color = type, shape = type)) +
    geom_text(data = plot.annotation, aes(label = label), hjust = 0, vjust = 1, size = 6) +
    scale_color_manual(values = c("black", "grey50")) +
    scale_shape_manual(values = c(21, 16)) +
    theme_ggEHD(base_size = BASE.SIZE) +
    guides(color = guide_legend(override.aes = list(size = 4, stroke = 1, shape = c(21, 16))), shape = FALSE) +
    theme(legend.position   = c(0.7, 0.1),
          legend.key.width  = unit(0.25, "lines"),
          legend.background = element_rect(fill = "transparent"),
          strip.background  = element_blank(),
          strip.text        = element_blank())
}

comp.fig.respective.el <- com_fig_respective(combinded.dat.melted, "Ellipsoid")
ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_3_Approach_Comparison-Respective_PAR-EL.jpg"),
              plot     = comp.fig.respective.el,
              scale    = 1.8,
              dpi      = 300)

comp.fig.respective.cy <- com_fig_respective(combinded.dat.melted, "Cylinder")
ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_S1_Approach_Comparison-Respective_PAR-CY.jpg"),
              plot     = comp.fig.respective.cy,
              scale    = 1.8,
              dpi      = 300)

##### RUE vs PAR #####
rue.par.data <- combinded.dat %>%
  dplyr::mutate(modeled.RUE = modeled.A / modeled.PAR) %>%
  dplyr::mutate(measured.RUE = measured.A / measured.PAR) %>%
  tidyr::gather(modeled.RUE, measured.RUE, key = "type", value = "value") %>%
  dplyr::mutate(type = factor(type, levels =  c("modeled.RUE", "measured.RUE"), labels = c("Modeled", "Measured"))) %>%
  dplyr::arrange(dplyr::desc(type)) %>%
  dplyr::mutate(respective.PAR = modeled.PAR)

rue.par.data$respective.PAR[which(rue.par.data$type == "Measured")] <- rue.par.data$measured.PAR[which(rue.par.data$type == "Measured")]

rue_com_fig_respective <- function(dat, c) {

  dat <- dplyr::filter(dat, canopy.shape == c)

  plot.annotation <- data.frame(timeframe      = timeframe.levels,
                                respective.PAR = 0,
                                value          = max(rue.par.data$value, na.rm = TRUE),
                                label          = paste(c("(a)", "(b)", "(c)", "(d)"), timeframe.levels))

  ggplot(dat, aes(x = respective.PAR, y = value)) +
    labs(y     = bquote('Daily RUE (mmol ' ~CO[2]~ ' mol ' ~PAR^-1~ ')'),
         x     = bquote('Daily incident PAR (mol ' ~m^-2~ ' '~day^-1*')'),
         color = "") +
    scale_x_continuous(sec.axis = sec_axis(~ ., labels = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ ., labels = NULL), limits = c(0, NA)) +
    facet_wrap(~timeframe, ncol = 1) +
    geom_point(size = 2, na.rm = T, aes(color = type, shape = type)) +
    geom_text(data = plot.annotation, aes(label = label, y = value + 6), hjust = 0, vjust = 1, size = 6) +
    scale_color_manual(values = c("black", "grey50")) +
    scale_shape_manual(values = c(21, 16)) +
    theme_ggEHD(base_size = BASE.SIZE) +
    guides(color = guide_legend(override.aes = list(size = 4, stroke = 1, shape = c(21, 16))), shape = FALSE) +
    theme(legend.position   = c(0.7, 0.1),
          legend.key.width  = unit(0.25, "lines"),
          legend.background = element_rect(fill = "transparent"),
          strip.background  = element_blank(),
          strip.text        = element_blank())
}

rue.comp.fig.respective.el <- rue_com_fig_respective(rue.par.data, "Ellipsoid")
ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_4_RUE_Approach_Comparison-Respective_PAR-EL.jpg"),
              plot     = rue.comp.fig.respective.el,
              scale    = 1.8,
              dpi      = 300)

rue.comp.fig.respective.cy <- rue_com_fig_respective(rue.par.data, "Cylinder")
ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_S2_RUE_Approach_Comparison-Respective_PAR-CY.jpg"),
              plot     = rue.comp.fig.respective.cy,
              scale    = 1.8,
              dpi      = 300)

##### COMBINE & COMPARE DIFFERENT AVERAGING APPROACHES - An #####
## MODELED vs. MEASURED MINUTELY A using daily orchard averages
combinded.dat.orchard.avg <- combinded.dat %>%
  dplyr::group_by(orchard.id, DATE, canopy.shape, overcast, timeframe) %>%
  dplyr::summarize(modeled.A  = mean(modeled.A,  na.rm = TRUE),
                   measured.A = mean(measured.A, na.rm = TRUE)) %>%
  dplyr::ungroup()

MvM.plot.data <- combinded.dat.orchard.avg

measured.minutely.A <- MvM.plot.data %>%
  dplyr::filter(timeframe == "Minutely") %>%
  dplyr::select(-modeled.A, -timeframe) %>%
  dplyr::rename(minutely.measured.A = measured.A)

MvM.minutely <- MvM.plot.data %>%
  dplyr::left_join(measured.minutely.A, by = c("orchard.id", "DATE", "canopy.shape", "overcast")) %>%
  tidyr::gather(measured.A, modeled.A, key = "type", value = "A") %>%
  dplyr::mutate(type = factor(type, levels = c("measured.A", "modeled.A"), labels = c("Measured PAR", "Modeled PAR")))

min.val <- min(ceiling(min(pmin(MvM.minutely$A, MvM.minutely$minutely.measured.A, na.rm = TRUE), na.rm = TRUE)), 0)
max.val <- ceiling(max(pmax(MvM.minutely$A, MvM.minutely$minutely.measured.A, na.rm = TRUE), na.rm = TRUE))

make_label <- function(x) paste0("y = ", round(x[1], 2), "x")

MvM_minutely_fig <- function(d, c) {

  d <- d %>% dplyr::filter(canopy.shape == c)

  A.eq <- d %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(type, timeframe) %>%
    tidyr::nest() %>%
    dplyr::mutate(model = purrr::map(data, ~lm(A ~ minutely.measured.A -1, data = .))) %>%
    dplyr::mutate(coefs = purrr::map(model, coef)) %>%
    dplyr::mutate(lab   = purrr::map_chr(coefs, make_label))

  ggplot(d, aes(x = minutely.measured.A, y = A)) +
    labs(x = bquote('Daily ' ~A[n]~ ' (mmol ' ~CO[2]~ ' '~m^-2~ ' '~day^-1*') from measured minutely PAR'),
         y = bquote('Daily ' ~A[n]~ ' (mmol ' ~CO[2]~ ' '~m^-2~ ' '~day^-1*') from respective PAR'),
         shape = NULL) +
    scale_x_continuous(limits = c(min.val, max.val)) +
    scale_y_continuous(limits = c(min.val, max.val)) +
    facet_grid(timeframe ~ type) +
    geom_point(size = 3, na.rm = T, aes(shape = orchard.id)) +
    scale_shape_manual(values = orchard.shapes) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    geom_smooth(color = "black", linetype = "solid", method = "lm", se = F, formula = y ~ x - 1) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    geom_text(data = A.eq, aes(x = max.val-10, y = min.val + 25, label = lab), hjust = 1, vjust = 0, size = 5, fontface = "italic") +
    theme_ggEHD(base_size = BASE.SIZE) +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.position   = c(0.13, 0.94),
          legend.key.width  = unit(0.25, "lines"),
          strip.text        = element_text(size = 18))
}

combinded.MvM.fig.el <- MvM_minutely_fig(MvM.minutely, "Ellipsoid")
ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_5_A_v_MINUTELY_Measured-EL.jpg"),
              plot     = combinded.MvM.fig.el,
              scale    = 1.7,
              dpi      = 300)

combinded.MvM.fig.cy <- MvM_minutely_fig(MvM.minutely, "Cylinder")
ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_S3_A_v_MINUTELY_Measured-CY.jpg"),
              plot     = combinded.MvM.fig.cy,
              scale    = 1.7,
              dpi      = 300)


##### RUE FIGURE MVM #####
MvM.RUE.data <- combinded.dat %>%
  dplyr::mutate(modeled.RUE  = modeled.A / modeled.PAR) %>%
  dplyr::mutate(measured.RUE = measured.A / measured.PAR) %>%
  dplyr::group_by(orchard.id, DATE, canopy.shape, overcast, timeframe) %>%
  dplyr::summarize(modeled.RUE  = mean(modeled.RUE,  na.rm = TRUE),
                   measured.RUE = mean(measured.RUE, na.rm = TRUE))

measured.minutely.RUE <- MvM.RUE.data %>%
  dplyr::filter(timeframe ==  "Minutely") %>%
  dplyr::select(-modeled.RUE, -timeframe) %>%
  dplyr::rename(minutely.measured.RUE = measured.RUE)

MvM.RUE.data <- MvM.RUE.data %>%
  left_join(measured.minutely.RUE, by = c("orchard.id", "DATE", "canopy.shape", "overcast")) %>%
  gather(measured.RUE, modeled.RUE, key = "type", value = "RUE") %>%
  mutate(type = factor(type, levels = c("measured.RUE", "modeled.RUE"), labels = c("Measured PAR", "Modeled PAR")))

min.val <- ceiling(min(pmin(MvM.RUE.data$RUE, MvM.RUE.data$minutely.measured.RUE, na.rm = T), na.rm = T))
max.val <- ceiling(max(pmax(MvM.RUE.data$RUE, MvM.RUE.data$minutely.measured.RUE, na.rm = T), na.rm = T))

make_RUE_label <- function(x) paste0("y = ", round(x[2], 2), "x + ", round(x[1], 2))

MvM_RUE_fig <- function(dat, c) {

  dat <- dat %>% dplyr::filter(canopy.shape == c)

  RUE.eq <- dat %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(type, timeframe) %>%
    tidyr::nest() %>%
    dplyr::mutate(model = map(data, ~lm(RUE ~ minutely.measured.RUE, data = .))) %>%
    dplyr::mutate(coefs = map(model, coef)) %>%
    dplyr::mutate(lab = map_chr(coefs, make_RUE_label))

  RUE.eq$lab <- stringr::str_remove(RUE.eq$lab, " \\+ 0")

  OE.text <- dat %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(type, timeframe) %>%
    dplyr::mutate(perc.oe = (RUE - minutely.measured.RUE) / minutely.measured.RUE * 100) %>%
    dplyr::summarize(perc.oe = mean(perc.oe)) %>%
    dplyr::mutate(lab = as.character(paste0(round(perc.oe), "% overestimation")))

  OE.text$lab[OE.text$perc.oe == 0] <- ""

  ggplot(dat, aes(x = minutely.measured.RUE, y = RUE)) +
    labs(x = bquote('Daily ' ~RUE~ ' (mmol ' ~CO[2]~ ' mol ' ~PAR^-1~ ') from measured minutely PAR'),
         y = bquote('Daily ' ~RUE~ ' (mmol ' ~CO[2]~ ' mol ' ~PAR^-1~ ') from respective PAR'),
         shape = NULL) +
    scale_x_continuous(limits = c(min.val, max.val)) +
    scale_y_continuous(limits = c(min.val, max.val)) +
    facet_grid(timeframe ~ type) +
    geom_point(size = 3, na.rm = T, aes(shape = orchard.id)) +
    scale_shape_manual(values = orchard.shapes) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    geom_smooth(color = "black", linetype = "solid", method = "lm", se = F) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    geom_text(data = RUE.eq,  aes(x = max.val - 0.5, y = min.val + 1, label = lab), hjust = 1, vjust = 0, size = 5, fontface = "italic") +
    geom_text(data = OE.text, aes(x = max.val - 0.5, y = min.val + 2.5, label = lab), hjust = 1, vjust = 0, size = 5, fontface = "italic") +
    theme_ggEHD(base_size = BASE.SIZE) +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.position   = c(0.13, 0.94),
          legend.key.width  = unit(0.25, "lines"),
          strip.text        = element_text(size = 18))
}

RUE.MvM.fig.el <- MvM_RUE_fig(MvM.RUE.data, "Ellipsoid")
ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_6_RUE_v_MINUTELY_Measured-EL.jpg"),
              plot     = RUE.MvM.fig.el,
              scale    = 1.7,
              dpi      = 300)

RUE.MvM.fig.cy <- MvM_RUE_fig(MvM.RUE.data, "Cylinder")
ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_S4_RUE_v_MINUTELY_Measured-CY.jpg"),
              plot     = RUE.MvM.fig.cy,
              scale    = 1.7,
              dpi      = 300)

##### EXAMPLE MINUTELY TIMESERIES FROM ONE ORCHARD ONE DAY ######
ex.dat <- dat %>%
  dplyr::filter(hour %in% 6:16) %>%
  dplyr::select(orchard.id, sensor.id, DATE, hour, TIME, canopy.shape, overcast, modeled.PAR, measured.PAR, incident.PAR) %>%
  tidyr::gather(incident.PAR, modeled.PAR, measured.PAR, key = "type", value = "value") %>%
  dplyr::mutate(type = as.character(type)) %>%
  dplyr::mutate(type = paste(type, canopy.shape, sep = "."))

temp <- as.data.frame(str_split_fixed(ex.dat$type, fixed("."), 2))
ex.dat$approach <- temp$V1
ex.dat$metric   <- temp$V2

ex.dat$timestamp <- as.POSIXct(strptime(paste(ex.dat$DATE, ex.dat$TIME), format = '%F %R'))

ex.plot.dat <- ex.dat %>%
  dplyr::filter(orchard.id == "Young A" & DATE == "2016-09-29" & sensor.id == 4)

plot.annotation <- data.frame(canopy.shape = c("Ellipsoid", "Cylinder"),
                              timestamp    = as.POSIXct(strptime(paste("2016-09-29", "06:00"), format = '%F %R')),
                              value        = 1500,
                              label        = c("(a) Ellipsoid", "(b) Cylinder"))

ex.plot.dat.above <- dplyr::filter(ex.plot.dat, approach == "incident")
ex.plot.dat.below <- ex.plot.dat %>%
  dplyr::filter(approach != "incident") %>%
  dplyr::mutate(approach = factor(approach, levels = c("modeled", "measured"), labels = c("Modeled", "Measured"))) %>%
  dplyr::arrange(dplyr::desc(approach))

example.timeseries <- ggplot(ex.plot.dat.below, aes(y = value, x = timestamp)) +
  labs(y     = bquote('PAR (' ~mu~ 'mol ' ~m^-2~ ' '~s^-1*')'),
       x     = "Time of day",
       color = NULL) +
  facet_wrap(~canopy.shape, ncol = 1) +
  scale_x_datetime(sec.axis = sec_axis(~ ., labels=NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., labels=NULL)) +
  geom_area(data = ex.plot.dat.above, size = 1, na.rm = T, fill = "grey80") +
  geom_point(size = 2, na.rm = T, aes(color = approach, shape = approach)) +
  scale_color_manual(values = c("black", "grey50")) +
  scale_shape_manual(values = c(21, 16)) +
  geom_text(data = plot.annotation, aes(label = label), hjust = 0, vjust = 1, size = 6) +
  theme_ggEHD(base_size = BASE.SIZE) +
  guides(color = guide_legend(override.aes = list(size = 4, stroke = 1, shape = c(21, 16))), shape = FALSE) +
  theme(legend.position   = c(0.18, 0.9),
        legend.key.width  = unit(0.25, "lines"),
        legend.box        = "vertical",
        legend.background = element_rect(fill = "transparent"),
        strip.background  = element_blank(),
        strip.text        = element_blank())

ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_2_Example_Minutely_PAR_Timeseries.jpg"),
              plot     = example.timeseries,
              scale    = 1.1,
              dpi      = 300)


##### INCIDENT PAR PLOT #####
incident <- dat %>%
  dplyr::filter(canopy.shape == "Ellipsoid") %>%
  dplyr::filter(sensor.id == 1) %>%
  dplyr::filter(hour %in% 6:16) %>%
  dplyr::mutate(orchard.date = paste0(orchard.id, "   ", DATE)) %>%
  dplyr::mutate(TIME = lubridate::ymd_hm(paste0("2000-1-1 ", TIME)))

for(i in 1:18) {
sensor.1 <- dat %>%
  dplyr::filter(canopy.shape == "Ellipsoid") %>%
  dplyr::filter(sensor.id == i) %>%
  dplyr::filter(hour %in% 6:16) %>%
  dplyr::mutate(orchard.date = paste0(orchard.id, "   ", DATE)) %>%
  dplyr::mutate(TIME = lubridate::ymd_hm(paste0("2000-1-1 ", TIME)))

incident.par.plot <- ggplot(incident, aes(x = TIME)) +
  labs(x = "Time of day",
       y = bquote('Incident PAR above orchard (' ~mu~ 'mol ' ~m^-2~ ' '~s^-1*')')) +
  facet_wrap(~orchard.date) +
  geom_area(fill = "grey80", aes(y = incident.PAR)) +
  geom_line(data = sensor.1, size = 1, aes(y = measured.PAR)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., labels=NULL)) +
  theme_ggEHD()

ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_S1_", i, "_Incident_PAR_Orchard_Comp.jpg"),
              plot     = incident.par.plot,
              scale    = 1.7,
              dpi      = 300)
}


##### SENSITIVITY ANALYSIS #####
## SENSITIVITY FUNCTION
sensitivity <- function(V, R, d) {

  J <- exp(log(V) * 0.890 + 1.010)

  calc_A <- function(d, V, J, R) {
    d %>%
      dplyr::mutate(modeled.A  = Farquhar_VJ(modeled.PAR,  V, J, R)) %>%
      dplyr::mutate(measured.A = Farquhar_VJ(measured.PAR, V, J, R))
  }

  cols.to.convert <- c("modeled.PAR", "measured.PAR", "modeled.A", "measured.A")

  ## CREATE DIFFERENTE TIME FRAMES
  d$timestamp <- as.POSIXct(d$timestamp)
  d <- dplyr::as_tibble(d) %>%
    dplyr::filter(canopy.shape == "Ellipsoid")

  min.dat <- d %>%
    dplyr::filter(hour %in% 6:16) %>%
    dplyr::select(orchard.id, sensor.id, DATE, hour, halfhour, TIME, canopy.shape, overcast, modeled.PAR, measured.PAR) %>%
    calc_A(V = V, J = J, R = R) ## CALCULATE PHOTOSYNTHESIS

  hour.dat <- min.dat %>%
    dplyr::group_by(orchard.id, sensor.id, DATE, hour, canopy.shape, overcast) %>%
    dplyr::summarize(modeled.PAR  = mean(modeled.PAR,  na.rm = TRUE),
                     measured.PAR = mean(measured.PAR, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    calc_A(V = V, J = J, R = R) ## CALCULATE PHOTOSYNTHESIS

  thirty.dat <- min.dat %>%
    dplyr::group_by(orchard.id, sensor.id, DATE, hour, halfhour, canopy.shape, overcast) %>%
    dplyr::summarize(modeled.PAR  = mean(modeled.PAR,  na.rm = TRUE),
                     measured.PAR = mean(measured.PAR, na.rm = TRUE)) %>%
    calc_A(V = V, J = J, R = R) ## CALCULATE PHOTOSYNTHESIS

  day.dat <- min.dat %>%
    dplyr::group_by(orchard.id, sensor.id, DATE, canopy.shape, overcast) %>%
    dplyr::summarize(modeled.PAR  = mean(modeled.PAR,  na.rm = TRUE),
                     measured.PAR = mean(measured.PAR, na.rm = TRUE)) %>%
    calc_A(V = V, J = J, R = R) ## CALCULATE PHOTOSYNTHESIS


  min.dat <- min.dat %>%
    dplyr::mutate_at(cols.to.convert, function(x) x * 60) ## CONVERT PAR & A FROM umol m-2 s-1 TO umol m-2 min-1

  hour.dat <- hour.dat %>%
    dplyr::mutate_at(cols.to.convert, function(x) x * 60 * 60) ## CONVERT PAR & A FROM umol m-2 min-1 TO umol m-2 hr-1

  thirty.dat <- thirty.dat %>%
    dplyr::mutate_at(cols.to.convert, function(x) x * 60 * 30) ## CONVERT PAR & A FROM umol m-2 min-1 TO umol m-2 half-hour-1

  day.dat <- day.dat %>%
    dplyr::mutate_at(cols.to.convert, function(x) x * 60 * 60 *11) ## CONVERT PAR & A FROM umol m-2 min-1 TO umol m-2 day-1

  ## AGGREGATE UP TO DAILY DATA FOR ALL AVERAGING TIMEFRAMES
  aggregate <- function(d) {
    d %>%
      dplyr::group_by(orchard.id, sensor.id, DATE, canopy.shape, overcast) %>%
      dplyr::summarize(modeled.PAR  = sum(modeled.PAR,  na.rm = TRUE),
                       measured.PAR = sum(measured.PAR, na.rm = TRUE),
                       modeled.A    = sum(modeled.A,    na.rm = TRUE),
                       measured.A   = sum(measured.A,   na.rm = TRUE)) %>%
      ungroup()
  }

  convert_units <- function(d) {
    d %>%
      dplyr::mutate(modeled.PAR  = modeled.PAR  / 1e6) %>% ## CONVERT UNITS TO mol  m-2 day-1 FOR PAR
      dplyr::mutate(measured.PAR = measured.PAR / 1e6) %>% ## CONVERT UNITS TO mol  m-2 day-1 FOR PAR
      dplyr::mutate(modeled.A    = modeled.A    / 1e3) %>% ## CONVERT UNITS TO mmol m-2 day-1 FOR PAR
      dplyr::mutate(measured.A   = measured.A   / 1e3)     ## CONVERT UNITS TO mmol m-2 day-1 FOR PAR
  }

  day.min.dat    <- min.dat    %>% aggregate() %>% convert_units()
  day.hour.dat   <- hour.dat   %>% aggregate() %>% convert_units()
  day.thirty.dat <- thirty.dat %>% aggregate() %>% convert_units()
  day.day.dat    <- day.dat    %>% convert_units()

  ### COMBINE
  timeframe.levels <- c("Minutely", "Half-hourly mean", "Hourly mean", "Daily mean")

  day.min.dat$timeframe    <- timeframe.levels[1]
  day.thirty.dat$timeframe <- timeframe.levels[2]
  day.hour.dat$timeframe   <- timeframe.levels[3]
  day.day.dat$timeframe    <- timeframe.levels[4]

  combinded.dat <- dplyr::bind_rows(day.min.dat, day.hour.dat, day.thirty.dat, day.day.dat) %>%
    mutate(timeframe    = factor(timeframe,    levels = timeframe.levels)) %>%
    mutate(canopy.shape = factor(canopy.shape, labels = "Ellipsoid"))

  ## A PLOT
  combinded.dat.orchard.avg <- combinded.dat %>%
    dplyr::group_by(orchard.id, DATE, canopy.shape, overcast, timeframe) %>%
    dplyr::summarize(modeled.A  = mean(modeled.A,  na.rm = TRUE),
                     measured.A = mean(measured.A, na.rm = TRUE))

  MvM.plot.data <- combinded.dat.orchard.avg

  measured.minutely.A <- MvM.plot.data %>%
    dplyr::filter(timeframe == "Minutely") %>%
    dplyr::select(-modeled.A, -timeframe) %>%
    dplyr::rename(minutely.measured.A = measured.A)

  MvM.minutely <- merge(MvM.plot.data, measured.minutely.A) %>%
    tidyr::gather(measured.A, modeled.A, key = "type", value = "A") %>%
    dplyr::mutate(type = factor(type, levels = c("measured.A", "modeled.A"), labels = c("Measured PAR", "Modeled PAR")))

  max.val <- ceiling(max(pmax(MvM.minutely$A, MvM.minutely$minutely.measured.A, na.rm = T), na.rm = T))

  out <- MvM.minutely %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(type, timeframe) %>%
    tidyr::nest() %>%
    dplyr::mutate(model = purrr::map(data, ~lm(A ~ minutely.measured.A -1, data = .))) %>%
    dplyr::mutate(coefs = purrr::map_dbl(model, coef)) %>%
    dplyr::select(type, timeframe, coefs) %>%
    dplyr::mutate(Vcmax = V) %>%
    dplyr::mutate(Jmax = J) %>%
    dplyr::mutate(Rd = R)

  return(out)
}

## RUN SENSITIVITY ANALYSIS ON Vcmax & Jmax
Vcmax.range <- 60:160
# 60:160 range is a "sensible range of winter wheat Vxmax per Sun et al 2015's data range
# 10:185 is the full range of possible Vcmax across all PFTs from Walker et al 2014.

if(REDO.SENSITIVITY) {
  OUT <- purrr::map_df(Vcmax.range, sensitivity, R = R, d = dat)
  readr::write_csv(OUT, paste0("./output/data/", "Chestnut_Light_Model_Sensitivity_Results.csv"))
} else {
  OUT <- readr::read_csv(paste0("./output/data/", "Chestnut_Light_Model_Sensitivity_Results.csv"), col_types = readr::cols())
}

## SENSITIVITY ANALYSIS PLOT - Vcmax Jmax Separately
sensitivity.data <- OUT %>%
  dplyr::mutate(timeframe = factor(timeframe, levels =  timeframe.levels, labels = timeframe.levels)) %>%
  dplyr::mutate(type = factor(type, levels = c("Modeled PAR", "Measured PAR"), labels = c("Modeled", "Measured"))) %>%
  dplyr::mutate(coefs = (coefs - 1) * 100) %>%
  dplyr::filter(!(type == "Measured" & timeframe == "Minutely"))

plot.annotation <- sensitivity.data %>%
  group_by(timeframe) %>%
  summarize(coefs = max(coefs)) %>%
  mutate(Vcmax = 60,
         Jmax  = 100,
         label = paste(c("(a)", "(b)", "(c)", "(d)"), timeframe.levels[1:4]))

sensitivity.plot <- ggplot(sensitivity.data, aes(x = Vcmax, y = coefs)) +
  labs(x = bquote(~V[cmax]~ ' (' ~mu~ 'mol '~m^-2~ ' '~s^-1*')'),
       y = bquote('Overestimation in daily ' ~A[n]~ ' (%)')) +
  facet_wrap(~timeframe, scales = "free_y", ncol = 1) +
  scale_x_continuous(limits = range(Vcmax.range)) +
  geom_line(size = 1.5, aes(color = type)) +
  scale_color_manual(values = c("black", "grey50")) +
  geom_vline(xintercept = V, linetype = "dashed") +
  geom_text(data = plot.annotation, aes(label = label, y = coefs * 1.15), hjust = 0, vjust = 1, size = 6) +
  theme_ggEHD(base_size   = BASE.SIZE) +
  theme(strip.background  = element_blank(),
        strip.text        = element_blank(),
        legend.title      = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.position   = c(0.27, 0.8))

ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_7_Sensitivity_Analysis-VCMAX.jpg"),
              plot     = sensitivity.plot,
              scale    = 1.6,
              dpi      = 300)


Vcalc <- function(V) exp(log(V) * 0.890 + 1.010)

sensitivity.plot.J <- ggplot(sensitivity.data, aes(x = Jmax, y = coefs)) +
  labs(x = bquote(~J[max]~ ' (' ~mu~ 'mol '~m^-2~ ' '~s^-1*')'),
       y = bquote('Overestimation in daily ' ~A[n]~ ' (%)')) +
  facet_wrap(~timeframe, scales = "free", ncol = 1) +
  scale_x_continuous(limits = Vcalc(range(Vcmax.range))) +
  geom_blank() +
  geom_text(data = plot.annotation, aes(label = label, y = coefs * 1.1), hjust = 0, vjust = 1, size = 6) +
  theme_ggEHD(base_size   = BASE.SIZE) +
  theme(strip.background  = element_blank(),
        axis.ticks.y      = element_blank(),
        strip.text        = element_blank())

ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_7_Sensitivity_Analysis-JMAX.jpg"),
              plot     = sensitivity.plot.J,
              scale    = 1.6,
              dpi      = 300)

## SENSITIVITY ANALYSIS SUMMARY FOR TEXT
for.text <- OUT %>%
  dplyr::mutate(coefs = (coefs - 1) * 100) %>%
  dplyr::filter(Jmax > 100) %>%
  dplyr::group_by(timeframe) %>%
  dplyr::summarize(min.coefs = min(coefs),
                   max.coefs = max(coefs))


## RUN SENSITIVITY ANALYSIS ON Rd
Rd.range <- seq(-4, 0, 0.1)
# -10:0 range is the full range of meausred data in Sun et al 2015 (Figure 1C)
# -4:0 is a more sensible range from Sun et al 2015 (Figure 2C) that accounts for effects of gm.

if(REDO.SENSITIVITY) {
  OUT.Rd <- purrr::map_df(Rd.range, sensitivity, V = V, d = dat)
  readr::write_csv(OUT.Rd, paste0("./output/data/", "Chestnut_Light_Model_Sensitivity_Results_Rd.csv"))
} else {
  OUT.Rd <- readr::read_csv(paste0("./output/data/", "Chestnut_Light_Model_Sensitivity_Results_Rd.csv"), col_types = readr::cols())
}

## SENSITIVITY ANALYSIS PLOT - Vcmax Jmax Separately
sensitivity.data.Rd <- OUT.Rd %>%
  dplyr::mutate(timeframe = factor(timeframe, levels =  timeframe.levels, labels = timeframe.levels)) %>%
  dplyr::mutate(type = factor(type, levels = c("Modeled PAR", "Measured PAR"), labels = c("Modeled", "Measured"))) %>%
  dplyr::mutate(coefs = (coefs - 1) * 100) %>%
  dplyr::filter(!(type == "Measured" & timeframe == "Minutely"))

plot.annotation.Rd <- sensitivity.data.Rd %>%
  group_by(timeframe) %>%
  summarize(coefs = max(coefs)) %>%
  mutate(Rd    = -4,
         label = paste(c("(a)", "(b)", "(c)", "(d)"), timeframe.levels[1:4]))

sensitivity.plot.Rd <- ggplot(sensitivity.data.Rd, aes(x = Rd, y = coefs)) +
  labs(x = bquote(~R[d]~ ' (' ~mu~ 'mol '~m^-2~ ' '~s^-1*')'),
       y = bquote('Overestimation in daily ' ~A[n]~ ' (%)')) +
  facet_wrap(~timeframe, scales = "free_y", ncol = 1) +
  scale_x_continuous(limits = range(Rd.range)) +
  geom_line(size = 1.5, aes(color = type)) +
  scale_color_manual(values = c("black", "grey50")) +
  geom_vline(xintercept = R, linetype = "dashed") +
  geom_text(data = plot.annotation.Rd, aes(label = label, y = coefs * 1.2), hjust = 0, vjust = 1, size = 6) +
  theme_ggEHD(base_size   = BASE.SIZE) +
  theme(strip.background  = element_blank(),
        strip.text        = element_blank(),
        legend.title      = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.position   = c(0.27, 0.8))

ggsave_fitmax(filename = paste0(outputPlotPath, "Fig_8_Sensitivity_Analysis-Rd.jpg"),
              plot     = sensitivity.plot.Rd,
              scale    = 1.6,
              dpi      = 300)
