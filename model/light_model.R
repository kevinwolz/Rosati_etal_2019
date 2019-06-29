##############################################################
# Light model
#
# Created Feb 22, 2017
# Lora Murphy (murphyl@caryinstitute.org)
#
# Modified Apr 20, 2017
# Kevin Wolz (kevin@savannainstitute.org)
##############################################################

light_model <- function(i, trees, points, latitude, azi_of_north, overcast){
  # Get the radiation array for this date and time - the fraction of
  # total sky radiation coming from each grid cell
  rad_array <- get_radiation_array(latitude, azi_of_north, overcast, points$DATE[i], points$TIME[i])

  # Add all the trees to the photo array, that says where
  # light is blocked by neighboring canopies
  photo_array <- get_photo_array(trees, points$X[i], points$Y[i], points$Z[i])

  # GLI is just the radiation array multiplied by the photo
  # array and summed
  points$GLI[i] <- sum(rad_array * photo_array)

  return(points[i,])
}