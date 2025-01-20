MatchSegmentationVolMac <- function(BScanHeader, Header, SegList = "WRT", CrSmooth = TRUE){
  
  ##### Assumes SLO is 768 x 768 pixels ########
  
  XRange <- range(c(BScanHeader$StartX, BScanHeader$EndX)/Header$ScaleX)
  YRange <- range(c(BScanHeader$StartY, BScanHeader$EndY)/Header$ScaleY)
  
  LengthX <- length(BScanHeader$StartX)
  BBoxX <- c(BScanHeader$StartX[c(1, LengthX)], rev(BScanHeader$EndX[c(1, LengthX)]))/Header$ScaleX
  BBoxY <- c(BScanHeader$StartY[c(1, LengthX)], rev(BScanHeader$EndY[c(1, LengthX)]))/Header$ScaleY 
  BBoxX <- c(BBoxX, BBoxX[1])
  BBoxY <- c(BBoxY, BBoxY[1])
  
  XWidth <- sqrt(diff(BBoxX[c(2,3)])^2 + diff(BBoxY[c(2,3)])^2)
  YWidth <- sqrt(diff(BBoxX[c(1,2)])^2 + diff(BBoxY[c(1,2)])^2)
  
  BBoxRot <- rad2deg(atan(diff((BBoxY[c(2,3)]))/diff((BBoxX[c(2,3)]))))
  if (BBoxRot > 45 & BBoxRot < 90){BBoxRot <- BBoxRot + 180}
  SEG_List <- list()
  TryStep <- 0
  ReCheck <- TRUE
  
  for (i in 1:length(SegList)){
    Limits <- switch(SegList[i], 
                     "WRT" = c("ILM", "BM"),
                     "NFL" = c("ILM","NFL"),
                     "GCL" = c("NFL","GCL"),
                     "IPL" = c("GCL","IPL"),
                     "OUTER" = c("IPL","BM"))
    
    ######################################################################################
    #Reshapes
    SEG <- flipdim(t(BScanHeader[[Limits[2]]] - BScanHeader[[Limits[1]]]),2)
    
    #Smooths map across slices
    if (CrSmooth){
      #Stores NA locations
      NAIdex <- is.na(SEG)
      #Fills NAs
      SEG <- inpaint(as.cimg(SEG), sigma = 2)
      #Selects cross direction
      AxisSm <- c("x", "y")[which.min(dim(SEG)[1:2])]
      #One dimensional blur
      SEG <- vanvliet(SEG, sigma=2, order=0, axis=AxisSm)
      #Replaces NaNs back in
      SEG[NAIdex] <- NaN
    }
    
    # Resizes to fit expected image size (768 x 768)
    SEG <- resize(as.cimg(SEG + 1), round(XWidth), round(YWidth), interpolation_type = 3)
    
    # Rotates to match scan angle
    SEG <- imrotate(SEG, BBoxRot)
    
    # Pads to match scan bounding box (might me larger than SLO image)
    if(XRange[1] > 0){SEG <- pad(SEG, XRange[1], "x", 1)}
    if(YRange[1] > 0){SEG <- pad(SEG, YRange[1], "y", 1)}
    
    # Shifts scan map so that the beginning of the scan aligns with the SLO
    # For example, if the min X and Y range is negative, it will shift the scan 
    # to the left, outside the boundaries of the SLO image
    SEG <- imshift(SEG, XRange[1], YRange[1])
    
    # Cuts out excess
    SEG <- imsub(SEG, x <= 768, y <= 768)
    
    # Pads again if necessary
    SEG <- pad(SEG, 768 - dim(SEG)[1], "x", pos = 1)
    SEG <- pad(SEG, 768 - dim(SEG)[2], "y", pos = 1)
    
    # Changes NA to NaN (works better with imager)
    SEG[is.na(SEG)] <- NaN
    # Sets < 0 to 0 (it can happen when taking differences between layers)
    SEG[SEG < 0] <- 0
    
    # Changes z-scale to microns
    SEG <- SEG*Header$ScaleZ*1000
    
    # Stores matched layers as lists
    SEG_List[[i]] <- SEG
  }
  
  # Final list 
  names(SEG_List) <- SegList
  # Sccan rotation
  SEG_List$Angle <- BBoxRot
  
  return(SEG_List)
}

