# Function to detect the fovea in an image
FoveaDetection <- function(WRT, CX = 385, CY = 385) {
  
  # Load the fovea template image
  load("Files/FoveaTemplate.RData")
  FoveaTemplate <- as.cimg(FoveaTemplate[91:271,91:271])
  
  # Crop the WRT image to a region of interest (ROI) around the expected fovea location
  VecCropX <- -160:160 + CX
  VecCropY <- -160:160 + CY
  WRT_Fovea <- as.cimg(WRT[VecCropX, VecCropY,,])
  
  # Replace any NA values in the cropped image with 0
  WRT_Fovea[is.na(WRT_Fovea)] <- 0
  
  # Perform cross-correlation between the cropped image and the fovea template
  # Normalize the cross-correlation values to be between 0 and 1
  CC <- correlate(WRT_Fovea, FoveaTemplate, normalise = TRUE)
  
  # Find the location of the maximum cross-correlation value, which corresponds to the fovea
  FoveaXY <- c(min(VecCropX), min(VecCropY)) + which(CC == max(CC), arr.ind = TRUE)[1:2]
  
  # Return the x and y coordinates of the detected fovea
  return(FoveaXY)
}