# Function to match segmentation of the optic nerve head (ONH)
MatchSegmentationVolONH <- function(BScanHeader, Header) {
  
  # Check the scan pattern
  if (Header$ScanPattern == 2) {
    # Old peripapillary scans
    
    # Calculate the ONH coordinates
    ONH_X <- BScanHeader$EndX/Header$ScaleXSlo
    ONH_Y <- BScanHeader$EndY/Header$ScaleYSlo
    ONH_RX <- BScanHeader$StartX/Header$ScaleXSlo
    ONH_RY <- BScanHeader$StartY/Header$ScaleYSlo
    
    # Calculate the ONH radius
    ONH_Rad <- sqrt((ONH_RX - ONH_X)^2 + (ONH_RY - ONH_Y)^2)
    
    # Create a vector of angles around the ONH
    ONH_Angles <- rev(seq(0, 360 - 360/768, length.out = 768))
    
    # Determine the scan position (OS or OD)
    if (Header$ScanPosition == "OS") {
      Inv <- 1
    } else {
      Inv <- -1
    }
    
    # Calculate the nerve fiber layer (NFL) thickness
    NFL_Thick <- (c(BScanHeader$NFL) - c(BScanHeader$ILM))*Header$ScaleZ*1000
    
    # Convert the polar coordinates to Cartesian coordinates
    NFL_Coords <- useful::pol2cart(theta = pracma::deg2rad(ONH_Angles), r = rep(ONH_Rad, 768))
    
    # Create the NFL segmentation list
    NFL_Seg <- list(x = Inv*NFL_Coords$x + ONH_X, y = NFL_Coords$y + ONH_Y,
                    NFL = NFL_Thick, theta = NFL_Coords$theta, ONH_X = ONH_X, ONH_Y = ONH_Y)
    
  } else if (Header$ScanPattern == 6) {
    # New ONH scans with glaucoma module - Assumes SLO is 768 x 768 pixels
    
    # Select the 25th scan. Compatible with previous scans (~6 degrees)
    Sel <- 25
    
    # Calculate the ONH coordinates
    ONH_X <- (BScanHeader$EndX/Header$ScaleXSlo)[Sel]
    ONH_Y <- (BScanHeader$EndY/Header$ScaleYSlo)[Sel]
    ONH_RX <- (BScanHeader$StartX/Header$ScaleXSlo)[Sel]
    ONH_RY <- (BScanHeader$StartY/Header$ScaleYSlo)[Sel]
    
    # Calculate the ONH radius
    ONH_Rad <- sqrt((ONH_RX - ONH_X)^2 + (ONH_RY - ONH_Y)^2)
    
    # Create a vector of angles around the ONH
    ONH_Angles <- rev(seq(0, 360 - 360/768, length.out = 768))
    
    # Determine the scan position (OS or OD)
    if (Header$ScanPosition == "OS") {
      Inv <- 1
    } else {
      Inv <- -1
    }
    
    # Calculate the nerve fiber layer (NFL) thickness
    NFL_Thick <- (c(BScanHeader$NFL[Sel,]) - c(BScanHeader$ILM[Sel,]))*Header$ScaleZ*1000
    
    # Convert the polar coordinates to Cartesian coordinates
    NFL_Coords <- useful::pol2cart(theta = pracma::deg2rad(ONH_Angles), r = rep(ONH_Rad, 768))
    
    # Create the NFL segmentation list
    NFL_Seg <- list(x = Inv*NFL_Coords$x + ONH_X, y = NFL_Coords$y + ONH_Y,
                    NFL = NFL_Thick, theta = NFL_Coords$theta, ONH_X = ONH_X, ONH_Y = ONH_Y)
  }
  
  # Return the NFL segmentation list
  return(NFL_Seg)
}