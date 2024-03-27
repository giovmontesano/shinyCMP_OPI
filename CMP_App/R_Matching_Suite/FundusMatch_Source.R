# Function to match a source image to a target image
FundusMatch <- function(source, target, CropSize = 480, type = "affine", Center = dim(target)[1:2]/2, tol = .2) {
  
  # Recenter the target image
  target <- imshift(target, dim(target)[1]/2 - Center[1], dim(target)[2]/2 - Center[2])
  
  # Crop the central location of the target image
  if (CropSize > 0) {
    targetC <- as.array(crop.borders(target, nx = CropSize, ny = CropSize))[,,1,1]
  } else {
    targetC <- as.array(target)[,,1,1]
  }
  
  # Resize the source image to match the size of the cropped target image
  sourceR <- as.array(resize(source, dim(targetC)[1], dim(targetC)[2]))[,,1,1]
  
  # Calculate the scale factor between the source and target images
  ScF <- dim(targetC)[1]/dim(source)[1]
  
  # Adjust the tolerance value to have a round CropSize
  if(CropSize > 0) {
    tol <- 1 - round(CropSize*(1 - tol))/CropSize
  } else {
    tol <- 0
  }
  
  # Crop the target image with the adjusted tolerance
  targetC <- as.array(crop.borders(target, nx = round(CropSize*(1 - tol)), 
                                   ny = round(CropSize*(1 - tol))))[,,1,1]
  
  # Perform an affine transformation to match the source image to the target image
  Match <- niftyreg(sourceR, targetC, scope = type, internal = FALSE) 
  
  # Store additional information in the Match object
  Match$CropTX <- CropSize
  Match$CropTY <- CropSize
  Match$Offset <- dim(target)[1:2]/2 - Center
  Match$ScF <- ScF
  Match$sourceR <- sourceR
  Match$targetC <- targetC
  Match$tol <- tol
  
  # Return the Match object
  return(Match)
}