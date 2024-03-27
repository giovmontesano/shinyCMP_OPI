# Function to open an OCT volume file
openVol <- function(path, quickFovea = FALSE, quickSegm = FALSE) {
  # Open the file for reading in binary mode
  fid <- file(path, "rb")
  
  # Initialize the header list
  Header <- list()
  
  # Read the header information from the file
  Header[[1]] <- readBin(fid, "integer", n = 12, size = 1)
  Header[[2]] <- readBin(fid, "integer", n = 1, size = 4)
  Header[[3]] <- readBin(fid, "integer", n = 1, size = 4)
  Header[[4]] <- readBin(fid, "integer", 1, size = 4)
  Header[[5]] <- readBin(fid, "double", 1)
  Header[[6]] <- readBin(fid, "double", 1)
  Header[[7]] <- readBin(fid, "double", 1)
  Header[[8]] <- readBin(fid, "integer", 1, size = 4)
  Header[[9]] <- readBin(fid, "integer", 1, size = 4)
  Header[[10]] <- readBin(fid, "double", 1)
  Header[[11]] <- readBin(fid, "double", 1)
  Header[[12]] <- readBin(fid, integer(32), 1)
  Header[[13]] <- readBin(fid, "double", 1)
  Header[[14]] <- rawToChar(readBin(fid, raw(), n = 4, size = 1, signed = F))
  Header[[15]] <- readBin(fid, "integer", n = 1, size = 8)
  Header[[16]] <- readBin(fid, "integer", n = 1, size = 4)
  Header[[17]] <- readBin(fid, "integer", n = 1, size = 4)
  Header[[18]] <- rawToChar(readBin(fid, raw(), n = 16, size = 1, signed = F))
  Header[[19]] <- rawToChar(readBin(fid, raw(), n = 16, size = 1, signed = F))
  Header[[20]] <- readBin(fid, "integer", n = 1, size = 4)
  Header[[21]] <- rawToChar(readBin(fid, raw(), n = 21, size = 1, signed = F))
  Header[[22]] <- readBin(fid, "integer", n = 3, size = 1)
  Header[[23]] <- readBin(fid, "double", n = 1)
  Header[[24]] <- readBin(fid, "integer", n = 1, size = 4)
  Header[[25]] <- rawToChar(readBin(fid, raw(), n = 24, size = 1, signed = F))
  Header[[26]] <- readBin(fid, "double", n = 1)
  Header[[27]] <- readBin(fid, "integer", n = 1, size = 4)
  Header[[28]] <- readBin(fid, "integer", n = 1, size = 4)
  #empty, skip
  readBin(fid, "integer", n = 1832, size = 1)
  
  # Name the header list elements
  names(Header) <- c("Version", "SizeX", "NumBScans", "SizeZ", "ScaleX", "Distance", "ScaleZ",
                     "SizeXSlo", "SizeYSlo", "ScaleXSlo", "ScaleYSlo", "FieldSizeSlo", "ScanFocus", "ScanPosition",
                     "ExamTime", "ScanPattern", "BScanHdrSize", "ID", "ReferenceID", "PID", "PatientID",
                     "Padding", "DOB", "VID", "VisitID", "VisitDate", "GridType","GridOffset")
  
  # Read the SLO image
  SLO <- matrix(readBin(fid, "integer", n = Header$SizeXSlo*Header$SizeYSlo, size = 1, signed = F),
                Header$SizeXSlo, Header$SizeYSlo)
  
  #Inverted to display correctly in imager
  BScans <- array(0, c(Header$SizeX, Header$SizeZ ,Header$NumBScans))
  
  # Initialize the B-scan array and the B-scan header list
  BScanHeader <- list()
  BScanHeader[[1]] <- numeric(Header$NumBScans)
  BScanHeader[[2]] <- numeric(Header$NumBScans)
  BScanHeader[[3]] <- numeric(Header$NumBScans)
  BScanHeader[[4]] <- numeric(Header$NumBScans)
  BScanHeader[[5]] <- numeric(Header$NumBScans)
  BScanHeader[[6]] <- numeric(Header$NumBScans)
  BScanHeader[[7]] <- numeric(Header$NumBScans)
  
  BScanHeader[[8]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[9]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[10]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[11]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[12]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[13]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[14]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[15]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[16]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[17]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  BScanHeader[[18]] <- matrix(NA, Header$NumBScans, Header$SizeX)
  
  
  names(BScanHeader) <- c("StartX","StartY","EndX","EndY","NumSeg","Quality","Shift",
                          "ILM","BM","NFL","GCL","IPL","INL","OPL","ONL","IS","OS","RPE")
  
  # Read the B-scan header information
  for (zz in 1:Header$NumBScans){
    
    #Resets origin
    seek(fid, 16+2048+(Header$SizeXSlo*Header$SizeYSlo)+
           ((zz-1)*(Header$BScanHdrSize+Header$SizeX*Header$SizeZ*4)), origin = "start")
    
    BScanHeader$StartX[zz] <- readBin(fid, "double", 1)
    BScanHeader$StartY[zz] <- readBin(fid, "double", 1)
    BScanHeader$EndX[zz] <- readBin(fid, "double", 1)
    BScanHeader$EndY[zz] <- readBin(fid, "double", 1)
    BScanHeader$NumSeg[zz] <- readBin(fid, "integer", 1, size = 4)
    
    #skip int32 (OffSeg)
    readBin(fid, "integer", 1, size = 4)
    
    BScanHeader$Quality[zz] <- readBin(fid, "numeric", 1, size = 4)
    BScanHeader$Shif[zz] <- readBin(fid, "integer", 1, size = 4)
    
    # Read the B-scan image data if not in quickFovea or quickSegm mode
    if (!quickFovea & !quickSegm){
      seek(fid, Header$BScanHdrSize+2048+(Header$SizeXSlo*Header$SizeYSlo)+
             ((zz - 1)*(Header$BScanHdrSize+Header$SizeX*Header$SizeZ*4)), origin = "start")
      oct <- matrix(readBin(fid, "numeric", Header$SizeX*Header$SizeZ, size = 4), 
                    Header$SizeX, Header$SizeZ)
      oct[oct > 1000] <- 0
      BScans[,,zz] <- oct^(1/4)
    }
    
    #Read segmentation images
    seek(fid, 256+2048+(Header$SizeXSlo*Header$SizeYSlo) +
           ((zz - 1)*(Header$BScanHdrSize + Header$SizeX*Header$SizeZ*4)), origin = "start");
    seg <- readBin(fid, "numeric", BScanHeader$NumSeg[zz]*Header$SizeX, size = 4)
    seg[seg > 496] <- NA
    
    # Store the segmentation data in the BScanHeader list
    if (BScanHeader$NumSeg[zz] >= 2) {
      q <- 1
      BScanHeader$ILM[zz,] <- seg[1:Header$SizeX];
      BScanHeader$BM[zz,] <- seg[(Header$SizeX*q+1):(Header$SizeX*(q+1))]
      
      if (BScanHeader$NumSeg[zz] > 2 & !quickFovea){
        #Required segmentations are not sequential
        ListSeg <- c(3:7,9,15:17)
        for (ss in 3:11){
          q <- ListSeg[ss-2] - 1
          BScanHeader[[ss + 7]][zz,] <- seg[(Header$SizeX*q+1):(Header$SizeX*(q+1))];
        }
      }
    }
  }
  
  # Close the file
  close(fid)
  
  # Return the data in a list
  if (!quickFovea) {
    ListOut <- list(SLO, Header, BScanHeader, BScans)
    names(ListOut) <- c("SLO", "Header", "BScanHeader", "BScans")
  } else {
    ListOut <- list(SLO, Header, BScanHeader)
    names(ListOut) <- c("SLO", "Header", "BScanHeader")
  }
  return(ListOut)
}


### Function to fix HR scans (rescales to 768 x 768) ###
FixScanSize <- function(OCT_Data){
  if (dim(OCT_Data$SLO)[1] > 768){
    
    # Store the initial size of the SLO image
    InitialSize <- dim(OCT_Data$SLO)
    
    # Rescale the SLO image to 768x768
    OCT_Data$SLO <- as.matrix(imresize(as.cimg(OCT_Data$SLO), 768/InitialSize[1]))
    
    # Rescales BSCANS
    suppressWarnings({
      OCT_Data$BScans <- as.array(resize(as.cimg(OCT_Data$BScans), size_x = 768, 
                                         size_y = dim(OCT_Data$BScans)[2], 
                                         dim(OCT_Data$BScans)[3]))
    })
    
    # Rescales SEGMENTATIONS
    ListSegm <- c("ILM","BM","NFL","GCL","IPL","INL","OPL","ONL","IS","OS","RPE")
    NScans <- dim(OCT_Data$BScanHeader$ILM)
    for (jj in 1:length(ListSegm)){
      OCT_Data$BScanHeader[[ListSegm[jj]]] <- as.matrix(resize(as.cimg(OCT_Data$BScanHeader[[ListSegm[jj]]]), 
                                                               size_x = NScans[1], size_y = 768))
    }
    
    # Update the header information
    OCT_Data$Header$SizeX <- 768
    OCT_Data$Header$SizeXSlo <- 768
    OCT_Data$Header$SizeXSlo <- 768
    OCT_Data$Header$ScaleX <- OCT_Data$Header$ScaleX*InitialSize[1]/768
    OCT_Data$Header$ScaleXSlo <- OCT_Data$Header$ScaleXSlo*InitialSize[1]/768
    OCT_Data$Header$ScaleYSlo <- OCT_Data$Header$ScaleYSlo*InitialSize[1]/768
    
  }
  
  return(OCT_Data)
}


