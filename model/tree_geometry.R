##############################################################
# Tree geometry
#
# Created Feb 22, 2017
# Lora Murphy (murphyl@caryinstitute.org)
#
# Modified May 28, 2017
# Kevin Wolz (kevin@savannainstitute.org)
##############################################################

#-----------------------------------------------------------------------#
# Things to calculate on sourcing
#-----------------------------------------------------------------------#
##
## Get the slope to all azimuth cells for adding cylindrical trees
##
# Make an array of the azimuth angle to the center of 1-degree cells,
# in radians
azi_cell_centers <- seq(from = 0.5, to = 359.5, by = 1) * (pi / 180.0)
all_azi_slope = 1/tan(azi_cell_centers)
rm(azi_cell_centers)

##
## Make a set of unit vectors to different azimuth and altitude settings
## for working with ellipsoidal shapes
##
unit_vectors <- matrix(list(c(1, 1, 1)), nrow = 90,  ncol = 360)
for (azi_d in 1:360) {
  # Make azi in radians, to center of 1 degree cell
  azi_r <- (azi_d - 0.5) * (pi / 180)

  for (alt in 1:90) {
    # Zenith angle needed, in radians; again, get the center of the
    # 1 degree cell
    zen <- (90.5 - alt) * (pi / 180)

    x <- sin(azi_r) * sin(zen)
    y <- cos(azi_r) * sin(zen)
    z <- cos(zen)

    unit_vectors[alt, azi_d][[1]] <- c(x, y, z)
  }
}
rm(alt, azi_d, azi_r, x, y, z, zen)



#-----------------------------------------------------------------------#
# Adds a set of trees to a photo array. The photo array is a grid where
# every 1-degree by 1-degree cell has the amount of light blocked by
# trees surrounding a particular point-of-view. 0 = all light blocked,
# 1 = no light blocked.
#
# Arguments
# trees: dataframe, one row per tree.
# X, Y: Coordinates
# TRANS: Crown light transmittance for each tree, from 0 (no
#        light transmitted) to 1 (all light transmitted).
# SHAPE: Shape of crown. "C" = cylinder, "E" = ellipsoid. (Paraboloid
#        not yet supported.) Shapes can be mixed.
# TOP: Top of crown, for trees that are cylinders. This can be NA
#      for other shapes.
# BOTTOM: Bottom of crown, for trees that are cylinders. This can
#         be NA for other shapes.
# RADIUS: Horizontal radius of crown.
# CENTER: Height off the ground of center of crown, for trees that
#         are ellipsoids. Can be NA for other shapes.
# VRAD: Vertical radius of crown, for trees that are ellipsoids. Can
#       be NA for other shapes.
#
# X, Y, Z: Coordinates of point-of-view for GLI photo.
#
# Returns: Photo array as a matrix. Rows = altitude, columns = azimuth.
#-----------------------------------------------------------------------#
get_photo_array <- function(trees, X, Y, Z) {

  # Declare photo array
  photo_array <- matrix(1, nrow = 90, ncol = 360)

  # Go through trees
  for (i in 1:nrow(trees)) {

    # What is the tree's shape?
    if (trees$SHAPE[i] == "C") { # Cylinder

      # Get coordinates of cylinder bottom RELATIVE TO POV
      bottom <- c(trees$X[i]      - X,
                  trees$Y[i]      - Y,
                  trees$BOTTOM[i] - Z)

      # Get cylinder top RELATIVE TO POV
      top <- trees$TOP[i] - Z

      p1 <- cylinder_photo_array(bottom        = bottom,
                                 top           = top,
                                 radius        = trees$RADIUS[i],
                                 transmittance = trees$TRANS[i])
      photo_array <- photo_array * p1

    } else if (trees$SHAPE[i] == "E") { # Ellipsoid

      # Get coordinates of center of ellipsoid RELATIVE TO POV
      center <- c(trees$X[i]      - X,
                  trees$Y[i]      - Y,
                  trees$CENTER[i] - Z)

      p1 <- ellipsoid_photo_array(center        = center,
                                  xrad          = trees$XRAD[i],
                                  yrad          = trees$YRAD[i], # X &Y radii of ellipsoids now distinct
                                  vrad          = trees$VRAD[i],
                                  transmittance = trees$TRANS[i])
      photo_array <- photo_array * p1
    }
  }
  return(photo_array)
}


#-----------------------------------------------------------------------#
# Get the azimuth angle of a point relative to the origin. Azimuth is 0
# north, east positive.
#
# This works on vectors.
#-----------------------------------------------------------------------#
get_azimuth <- function(X, Y) {

  azi <- rep(NA, length(X))

  #Calculate the azimuth - correct for quadrant of "to" relative to "from"
  #(counting clockwise from upper right)
  ## first quadrant
  azi <- ifelse (Y > 0 & X >= 0, atan(X / Y), azi)

  # second quadrant
  azi <- ifelse (X > 0 & Y <= 0, (pi / 2.0) + atan((-1.0 * Y) / X), azi)

  # third quadrant
  azi <- ifelse (X <= 0 & Y < 0, pi + atan((-1.0 * X) / (-1.0 * Y)), azi)

  # fourth quadrant
  azi <- ifelse (X < 0 & Y >= 0, ((1.5 * pi) + atan(Y / (-1.0 * X))), azi)

  return(azi)
}


#-----------------------------------------------------------------------#
# Gets the photo array with the blocking from a single cylinder shape.
#
# Assumptions:
# - The point-of-view of the photo array is at position 0, 0, 0
#
# Arguments:
# bottom: the center of the bottom endcap of the cylinder, as
#         a length 3 vector (X, Y, height of bottom). Height of
#        bottom might be negative if it is below the point-of-view.
# top: the Z value of the top of the cylinder. So height of
#      the cylinder is cylinder_top - cylinder_bottom[3].
# radius: radius of cylinder
# transmittance: light transmittance of canopy, from 0 (no light transmitted) to
#                1 (lets through all light)
#
# Returns: Photo array as a matrix. Rows = altitude, columns = azimuth.
# The cells blocked by the cylinder will have a value of [transmittance].
# All other cells will be 1.
#-----------------------------------------------------------------------#
cylinder_photo_array <- function(bottom, top, radius, transmittance) {

  photo_array <- matrix(1, nrow = 90,  ncol = 360)

  # Is the top below the point, so that it can never shade?
  if (top <= 0) return(photo_array)

  # Basic method here: Find the center of the cylinder, walk left/right until
  # we fall off the edge, finding the top/bottom of each 1-degree azimuth slice

  # Find distance to center of cylinder
  dist <- sqrt(bottom[1]^2 + bottom[2]^2)

  # Is this point within the canopy?
  if (dist < radius && bottom[3] <= 0) {
    # Yes - reduce light for the whole photo array
    photo_array <- photo_array * transmittance
    return(photo_array)
  }

  # Find the azimuth to the center of the cylinder
  azi <- get_azimuth(bottom[1], bottom[2])

  # Turn this into an array index, by converting to degrees
  home_azi_index = trunc(azi * 180 / pi)

  # Calculate true slope to tree to make sure we get at least one hit in
  # the photo array, in case the canopy is very very skinny
  if(bottom[1] == 0) {
    slope <- 100000.0
  } else {
    slope <- bottom[2] / bottom[1]
  }

  # Start our calcs here
  azi_index <- home_azi_index

  # This controls searching. We'll use it to decrement the azimuth index
  # until we run out of tree, then we can flip the sign to increment
  # our way out the other side of the tree. By checking the sign we can
  # tell if we've searched both sides or not
  azi_inc <- -1

  # Count the cells to make sure we don't increment more than once
  # around for a tree that's very close
  cell_counter <- 0
  keep_going <- TRUE

  while (keep_going) {

    # We have to use algebra to solve the system of equations of the circle of
    # the neighbor's canopy and the line from the target to the center of the
    # azimuth cell in question. Remember the POV is at the origin

    # Plug in for a, b, and c in the quadratic formula
    a = 1.0 + slope * slope
    b = -2.0 * (bottom[1] + slope * bottom[2])
    c = bottom[2] * bottom[2] + bottom[1] * bottom[1] - radius * radius

    # If the bit inside the square root is positive, it means that
    # the system of equations has a real solution, the circle and line
    # intersect, and thus the azimuth chunk is blocked.
    if ((b * b - 4.0 * a * c) >= 0.0) {

      # Finish solving the quadratic formula so we can calculate altitude angles.
      # We'll have an X,Y point for the slope's intersection with the near and
      # far sides of the neighbor's crown
      c <- sqrt(b * b - 4.0 * a * c)
      nearX <- (-b - c) / (2.0 * a)
      farX  <- (-b + c) / (2.0 * a)
      nearY <- slope * nearX
      farY  <- slope * farX

      # Calculate the squares of the distances to the near and far points on the
      # neighbor's crown - we don't need to bother taking the square roots
      square_near_dist <- nearX * nearX + nearY * nearY
      square_far_dist  <- farX  * farX  + farY  * farY

      # Switch them in case we were wrong about which point was near
      if (square_near_dist > square_far_dist) {
        a = square_near_dist
        square_near_dist = square_far_dist
        square_far_dist = a
      }

      # Calculate the starting and ending altitudes for this azimuth
      theta1 <- asin(bottom[3] / sqrt(square_far_dist + bottom[3]^2))
      theta2 <- asin(top / sqrt(square_near_dist + top^2))

      alt_start <- max(trunc(theta1 * 180 / pi), 0)
      alt_end   <- min(trunc(theta2 * 180 / pi), 90)

      # Now multiply the light extinction coefficient into the photo array for
      # the altitude cells for this azimuth angle
      photo_array[alt_start:alt_end, azi_index] <- photo_array[alt_start:alt_end, azi_index] * transmittance

    } # end of if ((b * b - 4.0 * a * c) >= 0.0)

    # If that azimuth cell wasn't blocked, we're done searching in this
    # direction.  If we've also searched the other direction, quit
    else {
      if (azi_inc == 1) {
        keep_going = FALSE

        # We still need to search the other direction
      } else {

        # Go back to the home azimuth cell - we'll increment past it below
        azi_index <- home_azi_index
        # Flip the incrementer to the other direction
        azi_inc <- 1
        # Start the chunk counter over
        cell_counter <- 0
      }
    }

    # Increment the azimuth chunk in the appropriate direction
    azi_index <- azi_index + azi_inc
    cell_counter <- cell_counter + 1

    # If we've done half the azimuth angles, quit even if we're still
    # finding solutions
    if (cell_counter == 180) keep_going = FALSE

    # If we've gone too far in either direction, wrap around
    if (azi_index < 1)
      azi_index <- 360
    if(azi_index == 361)
      azi_index <- 1

    slope <- all_azi_slope[azi_index]
  } # end of while (keep_going)

  return(photo_array)
}


#-----------------------------------------------------------------------#
# Gets the photo array with the blocking from a single ellipsoid shape.
#
# Assumptions:
# - The point-of-view of the photo array is at position 0, 0, 0
#
# Arguments:
# center: the center of the ellipsoid, as a length 3 vector (X, Y, Z [height]).
#         Z might be negative if it is below the point-of-view.
# hrad: the horizontal radius of the ellipsoid
# vrad: the vertical radius of the ellipsoid
# transmittance: light transmittance of canopy, from 0 (no light transmitted) to
#                1 (lets through all light)
#
# Returns: Photo array as a matrix. Rows = altitude, columns = azimuth.
# The cells blocked by the ellipsoid will have a value of [transmittance].
# All other cells will be 1.
#-----------------------------------------------------------------------#
ellipsoid_photo_array <- function(center, xrad, yrad, vrad, transmittance) {

  photo_array <- matrix(1, nrow = 90,  ncol = 360)

  # Is the top of the ellipsoid below the point, so that it can never shade?
  if (center[3] + vrad <= 0) return(photo_array)

  # Is this point within the canopy? Use ellipsoid equation x^2/a^2 + y^2/b^2 +
  # z^2/c^2 = 1. Plug in and if solution < 1 point is inside. (It's true that
  # this equation wants the ellipsoid centered at (0, 0, 0). But recentering
  # would have us subtracting ellipsoid coordinates from our POV at (0, 0, 0),
  # and what with squaring them at all, it comes out exactly the same.)
  if (center[1]^2 / xrad^2 + center[2]^2 / yrad^2 + center[3]^2 / vrad^2 < 1) {
    # Yes - reduce light for the whole photo array
    photo_array <- photo_array * transmittance
    return(photo_array)
  }

  # Basic method here: Find the center of the ellipsoid, walk left/right until
  # we fall off the edge, walking up/down each 1-degree azimuth slice

  # Find the azimuth to the center of the ellipsoid
  azi <- get_azimuth(center[1], center[2])

  # Turn this into an array index, by converting to degrees
  home_azi_index = trunc(azi * 180/pi)
  if (home_azi_index == 0) {
    home_azi_index = 360
  }

  # Start our calcs here
  azi_index <- home_azi_index

  # This controls searching. We'll use it to decrement the azimuth index
  # until we run out of tree, then we can flip the sign to increment
  # our way out the other side of the tree. By checking the sign we can
  # tell if we've searched both sides or not
  azi_inc <- -1

  # Count the cells to make sure we don't increment more than once
  # around for a tree that's very close
  cell_counter <- 0
  keep_going <- TRUE

  while (keep_going) {

    # Count up the altitude angles from 0 to 90. Keep a flag to tell
    # us if we got a hit
    hit <- FALSE

    for (alt in 1:90) {

      if (ellipsoid_intersect(unit_vectors[alt, azi_index][[1]], # ray direction
                              center, xrad, yrad, vrad)) {

        hit <- TRUE

        # Now multiply the light extinction coefficient into the photo array for
        # this cell
        # Note special assignment operator - makes it happen in the global scope
        photo_array[alt, azi_index] <- photo_array[alt, azi_index] * transmittance

      }
    } # end of for (alt in 1:90)

    # If no altitude cells were blocked, we're done searching in this
    # direction.  If we've also searched the other direction, quit
    if (!hit) {
      if (azi_inc == 1) {
        keep_going = FALSE

        # We still need to search the other direction
      } else {

        # Go back to the home azimuth cell - we'll increment past it below
        azi_index <- home_azi_index
        # Flip the incrementer to the other direction
        azi_inc <- 1
        # Start the chunk counter over
        cell_counter <- 0
      }
    }

    # Increment the azimuth chunk in the appropriate direction
    azi_index <- azi_index + azi_inc
    cell_counter <- cell_counter + 1

    # If we've done half the azimuth angles, quit even if we're still
    # finding solutions
    if (cell_counter == 180) keep_going = FALSE

    # If we've gone too far in either direction, wrap around
    if (azi_index < 1)
      azi_index <- 360
    if(azi_index == 361)
      azi_index <- 1
  } # end of while (keep_going)

  return(photo_array)
}


#-----------------------------------------------------------------------#
# Function for determining the intersection of a line and an ellipsoid.
# Note that this is a line, not a ray, although it uses ray
# notation. If the ray is pointing away from the ellipsoid, but the
# direction vector is on a line that intersects the ellipsoid, this will
# return TRUE. Beware!
#
# Arguments:
# ray_origin: origin of ray, as length 3 vector of X,Y,Z coords
# ray_direction: direction vector of ray as length 3 vector
#                (normalization not required)
# center: coordinates of ellipsoid center, as length 3 vector
#         of X,Y,Z coords
# hrad: horizontal radius of ellipsoid (ellipsoid XY cross section assumed
#       to be a circle)
# vrad: vertical radius of ellipsoid
# Returns: TRUE if there is an intersection, FALSE otherwise
#-----------------------------------------------------------------------#
ellipsoid_intersect <- function(ray_direction, center, xrad, yrad, vrad) {
  rads <- c(xrad, yrad, vrad)^2
  a <- sum(ray_direction^2 / rads)
  b <- -2 * sum(ray_direction * center / rads)
  c <- sum(center^2 / rads) - 1

  # We don't care what the solutions actually are, just whether or
  # not they exist. So check the value inside the square root of
  # the quadratic formula. If it's 0 or positive, real solutions exist.
  if ((b^2 - 4 * a * c) >= 0) return(TRUE)

  return(FALSE)
}
