# Function to match segmentation of OCT images
# (Assumes SLO is 768 x 768 pixels)
MatchSegmentation <- function(BScanHeader, Header, SegList = "WRT", CrSmooth = TRUE){
  
  # Calculate the x and y ranges of the OCT image
  XRange <- range(c(BScanHeader$StartX, BScanHeader$EndX)/Header$ScaleX)
  YRange <- range(c(BScanHeader$StartY, BScanHeader$EndY)/Header$ScaleY)
  
  # Calculate the length of the OCT image
  LengthX <- length(BScanHeader$StartX)
  
  # Calculate the bounding box of the OCT image
  BBoxX <- c(BScanHeader$StartX[c(1, LengthX)], rev(BScanHeader$EndX[c(1, LengthX)]))/Header$ScaleX
  BBoxY <- c(BScanHeader$StartY[c(1, LengthX)], rev(BScanHeader$EndY[c(1, LengthX)]))/Header$ScaleY 
  BBoxX <- c(BBoxX, BBoxX[1])
  BBoxY <- c(BBoxY, BBoxY[1])
  
  # Calculate the width and height of the bounding box
  XWidth <- sqrt(diff(BBoxX[c(2,3)])^2 + diff(BBoxY[c(2,3)])^2)
  YWidth <- sqrt(diff(BBoxX[c(1,2)])^2 + diff(BBoxY[c(1,2)])^2)
  
  # Calculate the rotation angle of the bounding bo
  BBoxRot <- pracma::rad2deg(atan(diff((BBoxY[c(2,3)]))/diff((BBoxX[c(2,3)]))))
  
  # Initialize an empty list to store the segmentation results
  SEG_List <- list()
  
  # Initialize variables for the loop
  TryStep <- 0
  ReCheck <- TRUE
  
  # Loop through the segmentation layers
  for (i in 1:length(SegList)){
    Limits <- switch(SegList[i], 
                     "WRT" = c("ILM", "BM"),
                     "NFL" = c("ILM","NFL"),
                     "GCL" = c("NFL","GCL"),
                     "IPL" = c("GCL","IPL"),
                     "OUTER" = c("IPL","BM"))
    
    # Reshape the segmentation data
    SEG <- flipdim(t(BScanHeader[[Limits[2]]] - BScanHeader[[Limits[1]]]),2)
    
    # Smooth the segmentation data across slices
    if (CrSmooth){
      # Store the locations of NAs
      NAIdex <- is.na(SEG)
      # Fill in the NAs
      SEG <- inpaint(as.cimg(SEG), sigma = 2)
      # Select the axis for smoothing
      AxisSm <- c("x", "y")[which.min(dim(SEG)[1:2])]
      # Apply one-dimensional blur
      SEG <- vanvliet(SEG, sigma=2, order=0, axis=AxisSm)
      # Replace the NAs back in
      SEG[NAIdex] <- NaN
    }
    
    # Resize the segmentation data to the bounding box size
    SEG <- resize(as.cimg(SEG + 1), round(XWidth), round(YWidth), interpolation_type = 3)
    
    ### Following code is necessary as imrotate does weird approximations ###
    ###### resorted to rounding the expected size after rotation and calculating back ##
    ###### the angle needed for a round result. 
    w <- width(SEG)
    h <- height(SEG)
    
    Check <- FALSE
    
    #Checks different rounding methods (floor, round, ceil), to see which one works
    while (!Check){
      if (ReCheck){TryStep <- TryStep + 1 }
      
      if (TryStep == 1){
        W <- floor(w*cos(pracma::deg2rad((BBoxRot))) + h*sin(pracma::deg2rad((BBoxRot))))
      }
      
      if (TryStep == 2){
        W <- round(w*cos(pracma::deg2rad((BBoxRot))) + h*sin(pracma::deg2rad((BBoxRot))))
      }
      
      if (TryStep == 3){
        W <- ceil(w*cos(pracma::deg2rad((BBoxRot))) + h*sin(pracma::deg2rad((BBoxRot))))
      }
      
      # Rotate the segmentation data
      COS <- (W*w + h*(- W^2 + h^2 + w^2)^(1/2))/(h^2 + w^2)
      SEGR <- imrotate(SEG, sign(BBoxRot)*(pracma::rad2deg(acos(COS))), interpolation = 0, boundary = 0) - 1 
      
      # Handle negative values in the segmentation data
      SEGR[SEGR < 0] <- NaN
      
      ######################################################################################
      
      # Determine the x and y limits of the segmentation data
      XLimSEG <- round(XRange)
      XLimSEG[XLimSEG < 0] <- abs(XLimSEG[XLimSEG < 0]) - 1
      if(sum(XLimSEG < 1)){XLimSEG <- XLimSEG + 1}
      XLimSEG[2] <- XLimSEG[1] + min(767, diff(XLimSEG))
      if(XRange[1] > 1){XLimSEG <- XLimSEG - XLimSEG[1] + c(1,0)}
      
      YLimSEG <- round(YRange)
      YLimSEG[YLimSEG < 0] <- abs(YLimSEG[YLimSEG < 0]) - 1
      if(sum(YLimSEG < 1)){YLimSEG <- YLimSEG + 1}
      YLimSEG[2] <- YLimSEG[1] + min(767, diff(YLimSEG))
      if(YRange[1] > 1){YLimSEG <- YLimSEG - YLimSEG[1] + c(1,0)}
      
      XLimSEG <- XLimSEG[1]:XLimSEG[2]
      YLimSEG <- YLimSEG[1]:YLimSEG[2]
      
      # Check if the segmentation data is within the image bounds
      Check <- sum(c(min(XLimSEG), min(YLimSEG)) > 0) & sum(c(max(XLimSEG) <= dim(SEG)[1], max(YLimSEG) <= dim(SEGR)[2]))
      if (Check){
        ReCheck <- FALSE
        SEG <- SEGR
      }
    }
    
    # Finalize the segmentation data
    SEG <- as.array(SEG, dim(SEG),1,1)
    SEG <- as.cimg(SEG[c(XLimSEG), c(YLimSEG),1,1])
    SEG <- pad(SEG, round(XRange[1]), pos = -1, axes = "x",  val = NaN)
    SEG <- pad(SEG, 768 - round(dim(SEG)[1]), 1, axes = "x", val = NaN)
    SEG <- pad(SEG, round(YRange[1]), pos = -1, axes = "y",  val = NaN)
    SEG <- pad(SEG, 768 - round(dim(SEG)[2]), 1, axes = "y", val = NaN)
    
    # Handle negative and NaN values in the segmentation data
    SEG[is.na(SEG)] <- NaN
    SEG[SEG < 0] <- 0
    
    # Scale the segmentation data to micrometers
    SEG <- SEG*Header$ScaleZ*1000
    
    # Add the segmentation data to the list
    SEG_List[[i]] <- SEG
  }
  
  # Name the segmentation layers in the list
  names(SEG_List) <- SegList
  
  # Return the segmentation list
  return(SEG_List)
}


