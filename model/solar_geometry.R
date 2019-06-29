##############################################################
# Solar geometry
#
# Azimuth: 0 = north, increasing clockwise, always positive, 0-2pi
#
# Created Feb 2, 2017
# Lora Murphy (murphyl@caryinstitute.org)
#
# Modified May 26, 2017
# Kevin Wolz (kevin@savannainstitute.org)
##############################################################
library(solaR)

#---------------------------------------------------------------#
# Gets the sky radiation array. This is the amount of the total
# sky's radiation that appears in each 1-degree by 1-degree sky
# grid cell (the sum of the radiation array is 1). Radiation is
# split into two components: direct and diffuse. Direct (beam)
# radiation is what is coming from the sun, and all of it is
# placed in the grid cell in which the sun is found. Diffuse
# radiation is isotropic and is added to all cells equally. Since
# the cells are not equal in area, the diffuse radiation is scaled
# based on cell area.
#
# Arguments:
# latitude: Latitude.
# azi_of_north: The azimuth of north for the solar sky. The plot's
#               north is always azimuth 0. If this is set to
#               something other than 0, the sky is rotated
#               relative to the plot. Values are between 0 and
#               359 degrees.
# date: Date for which to get radiation array, as "YYYY-MM-DD".
# time: Time for which to get radiation array, as "HH:MM" (24 hour).
#
# Returns:
# Radiation array as a matrix. Rows = altitude, columns = azimuth.
# Will be all 0s if it's nighttime.
#---------------------------------------------------------------#
get_radiation_array <- function(latitude, azi_of_north, overcast, date, time) {

  # Create the sky array where we will put radiation values.
  radiation_array <- matrix(0, nrow = 90, ncol = 360)

  # Date of solar calculation, as a POSIXct object
  soldat <- as.POSIXct(date)

  # Calculate daily apparent movement of sun for chosen latitude and date
  solD <- fSolD(latitude, soldat)

  # Compute the angles for the intradaily movement for specific times. This
  # does the whole day at once, so we have to do it on a minute basis and
  # then extract out the time we want.
  solI <- fSolI(solD, sample = "min", keep.night = FALSE)

  # Find the index of the desired time
  time.index <- which(index(solI) == paste0(date, " ", time, ":00 UTC"))
  if (length(time.index) != 1) return(radiation_array)

  # This is an argument with total daily global irradiation. It's a zoo
  # object whose time index matches the day we're working with. The value of
  # 5000 doesn't matter because we are working with relative amounts of
  # radiation rather than total values.
  g0d <- zoo(5000, soldat)

  # Compute the diffuse and direct components from daily global irradiation
  # on a horizontal surface.
  compD <- fCompD(solD, g0d, corr = "CPR")

  # Use the daily proportions to calculate intradaily direct and diffuse
  # amounts for all the times in solI
  compI <- fCompI(solI, compD)

  convert_to_radians <- pi / 180
  convert_to_degrees <- 180 / pi

  #-------------------------------------------------------------#
  # Add diffuse radiation to solar array
  # Diffuse radiation is isotropic. Since our grid cells are on a hemisphere,
  # they're not equal in area. In order to equally divide diffuse radiation
  # among the cells, we have to do it on an equal-area basis.
  #-------------------------------------------------------------#

  # Grid cell area differs by altitude but not by azimuth. So take a single
  # row of cells from horizon to zenith, and get the area of each. On a unit
  # sphere, the area of a cell is (azi2 - azi1)(sin(alt2)-sin(alt1)). Since
  # azi is the same for each cell, we can drop that term. Calculate the
  # sin(alt2)-sin(alt1) bit for each cell.
  alt_areas <- rep(0, 90)
  for (i in 1:90) {
    alt_areas[i] <- sin(convert_to_radians * i) - sin(convert_to_radians * (i - 1))
  }

  # Relativize areas so each is a proportion of total sky area
  alt_areas <- alt_areas / (sum(alt_areas) * 360)

  # Add the area into the radiation array so we can do easy matrix multiplication
  for (i in 1:360) {
    radiation_array[,i] <- alt_areas
  }

  #-------------------------------------------------------------#
  # If the day is overcast, get the global (direct + diffuse) radiation, and apply this across the whole sky evenly
  # A completely overcast sky will scatter ALL radiation (direct + diffuse) evenly
  #-------------------------------------------------------------#
  if(overcast == "Y"){
    tot_radiation   <- compI[time.index, "G0"]
    radiation_array <- radiation_array * as.numeric(tot_radiation)
  } else{

    # Continue with calculating diffuse and direct beam radiation fractions if the day is NOT overcast
    # Get the amount of diffuse radiation
    tot_diffuse <- compI[time.index, "D0"]

    # Multiply it into the array. Have to use the as.numeric or the dimensions
    # collapse.
    radiation_array <- radiation_array * as.numeric(tot_diffuse)

    #-------------------------------------------------------------#
    # Add direct beam radiation to solar array. All of it goes into
    # whatever cell the sun is in.
    #-------------------------------------------------------------#
    # Azimuth of sun. Eastward is negative, south is 0, and westward is positive.
    adjust_azi <- function(azi) {
      ifelse(azi <= 0, pi - abs(azi), pi + azi)
    }
    # Convert to an azimuth index in degrees
    azi <- trunc(adjust_azi(solI[time.index, "AzS"]) * convert_to_degrees)

    # Altitude of sun. Convert to index in degrees
    alt <- trunc(solI[time.index, "AlS"] * convert_to_degrees)

    # Put all of the direct radiation into that cell
    radiation_array[alt, azi] <- compI[time.index, "B0"]
  }

  #-------------------------------------------------------------#
  # Relativize array
  #-------------------------------------------------------------#
  radiation_array <- radiation_array / sum(radiation_array)

  #-------------------------------------------------------------#
  # Rotate radiation array if azimuth > 0
  #-------------------------------------------------------------#
  if (azi_of_north > 0) {

    # Copy the radiation array
    rad_copy <- radiation_array

    # Just in case this isn't an integer - truncate to determine
    # columns of offset
    offset <- trunc(azi_of_north)

    # Offset by that many columns
    for (j in 1:360) {
      new_col <- j + offset
      if (new_col > 360) new_col <- new_col - 360
      radiation_array[, new_col] <- rad_copy[, j]
    }

    rm(rad_copy)
  }
  return(radiation_array)
}
