# Function to map a map image to a CMP fundus image
SLO_to_CMP_Maps <- function(Map, CMP_Fundus, Match, fullCMP = TRUE) {
  
  # Resize the map image to match the size of the source image in the Match object
  MapR <- as.array(resize(Map, dim(Match$sourceR)[1], dim(Match$sourceR)[2]))[,,1,1]
  
  # Apply the forward transformation from the Match object to the resized map image
  MapR <- as.cimg(applyTransform(forward(Match), MapR))
  
  # If fullCMP is TRUE, pad the transformed map image and shift it to the original CMP fundus position
  if(fullCMP) {
    MapR <- pad(MapR + 1, Match$CropTX*2*(1 - Match$tol), 0, axes = "xy", val = NaN) - 1
    MapR <- imshift(MapR, -Match$Offset[1], -Match$Offset[2])
    MapR[MapR < 0] <- NaN
  }
  
  # Return the transformed map image
  return(MapR)
}

# Function to transform OCT coordinates to CMP coordinates
OCT_to_CMP <- function(X, Y, Match) {
  # Create a matrix of the input X and Y coordinates
  XY <- as.matrix(data.frame(X = as.numeric(X), Y = as.numeric(Y)))
  
  # Apply the forward transformation from the Match object to the input coordinates
  Out <- (Match$CropTX*(1 - Match$tol)) + applyTransform(forward(Match), XY*Match$ScF)
  
  # Shift the transformed coordinates to the original CMP fundus position
  Out <- matrix(Out, ncol = 2)
  Out[,1] <- Out[,1] - Match$Offset[1]
  Out[,2] <- Out[,2] - Match$Offset[2]
  
  # Return the transformed coordinates as a data frame
  colnames(Out) <- c("X", "Y")
  return(as.data.frame(Out))
}

# Function to transform CMP coordinates to OCT coordinates
CMP_to_OCT <- function(X, Y, Match) {
  # Create a data frame of the input X and Y coordinates with the offset added
  XY <- data.frame(X = as.numeric(X + Match$Offset[1]), 
                   Y = as.numeric(Y + Match$Offset[2]))
  XY <- as.matrix(XY - Match$CropTX*(1 - Match$tol))
  
  # Apply the reverse transformation from the Match object to the input coordinates
  Out <- applyTransform(reverse(Match), XY)/Match$ScF
  
  # Return the transformed coordinates as a data frame
  Out <- matrix(Out, ncol = 2)
  colnames(Out) <- c("X", "Y")
  return(as.data.frame(Out))
}