# Function to overlay maps on a fundus image
Overlay_Maps <- function(Fundus, MapR, Replace = FALSE, LimMap = c(200, 400), Pedestal = TRUE) {
  
  # Load the Parula color map
  load("R_Matching_Suite/Parula_Map")
  
  # Handle NaN and NA values in the MapR image
  MapR[is.nan(MapR)] <- 0
  MapR[is.na(MapR)] <- 0
  
  # Define a function to scale the color values
  cscale <- function(v) {
    CMap <- c(hcl.colors(64, "gray"), as.character(Parula_Map$Hex))
    v <- round(v*length(CMap))
    v[v < 1] <- 1
    as.character(CMap[v])
  }
  
  # Normalize the Fundus image to a range of 0 to 0.5
  Final_Fundus <- as.cimg(as.array(Fundus[,,1,1]))
  Final_Fundus <- (Final_Fundus - min(Final_Fundus, na.rm = TRUE))/diff(range(Final_Fundus, na.rm = TRUE))*.5
  
  # Determine the offset for the MapR image
  Offset <- .5 - min(Final_Fundus[MapR > 0])
  
  # Normalize the MapR image to a range of 0.05 to 0.75
  Final_MapR <- (MapR - LimMap[1])/diff(range(LimMap))*(.75 - Offset) + Offset
  Final_MapR[MapR == 0] <- 0
  
  # Replace the Fundus image with the MapR image if specified
  if (Replace) {
    if (Pedestal) {
      Final_Fundus[MapR != 0] <- .5
      Final_MapR <- (MapR - LimMap[1])/diff(range(LimMap))*(.5 - .05) + .05
      Final_MapR[MapR == 0] <- 0
    } else {
      Final_Fundus[MapR != 0] <- 0
      Final_MapR <- (MapR - min(MapR, na.rm = TRUE))/diff(range(MapR,  na.rm = TRUE))*.5
    }
  }
  
  # Combine the Fundus and MapR images
  Final_Im <- Final_Fundus + Final_MapR
  
  # Convert the combined image to a raster and handle NAs
  Final_Im <- as.raster(Final_Im, colourscale = cscale, rescale = FALSE)
  Final_Im[is.na(Final_Im)] <- cscale(1)
  
  # Return the final overlaid image
  return(Final_Im)
}

# Function to plot a color bar
color.barv <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  # Create a plot with the color bar
  plot(c(0,1), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,1,y+1/scale, col=lut[i], border=NA)
  }
}

# Function to plot a horizontal color bar
color.barh <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  # Create a plot with the horizontal color bar
  plot(c(min,max), c(0,1), type='n', bty='n', xaxt='n', yaxt='n', ylab='', xlab=title)
  axis(1, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(y, 0, y+1/scale, 1, col=lut[i], border=NA)
  }
}