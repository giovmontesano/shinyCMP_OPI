
######################### Global functions ############################
#Get number of presentations
getPres <- function(state) {
  NPres <- state$numPresentations
  return(NPres)
}

#Get number of seen
getYes <- function(state) {
  NPres <- sum(state$responses)
  return(NPres)
}

#Get number of not seen
getNo <- function(state) {
  NPres <- sum(!state$responses)
  return(NPres)
}

#Extract response info
RespInfoExtract <- function(RespInfo){
  RespInfoMat <- NULL
  ColNames <- names(RespInfo[[1]][[1]])
  ColNames <- ColNames[ColNames != "err"]
  for (i in 1:length(RespInfo)){
    MT <- matrix(unlist(RespInfo[[i]]), nrow = length(RespInfo[[i]]), byrow = TRUE, 
                 dimnames = list(NULL, ColNames))
    RespInfoMat <- rbind(RespInfoMat, cbind(Location = i, MT))
  }
  return(RespInfoMat)
}

#Calculate normative values (based on CMP database)
CalcNorm <- function(Age, Ecc){
  NormModel <- read.csv("Files/NormModel.csv")
  NormData <- NormModel$Int + NormModel$AgeSl*Age
  NormVal<- approx(x = sort(unique(NormModel$Ecc)), 
                   y = tapply(NormData, NormModel$Ecc, mean), 
                   xout = Ecc)$y
  NormVal[is.na(NormVal)] <- max(c(0, NormVal), na.rm = TRUE)
  
  return(NormVal)
}

# Calculate neighborhood (only works for regular grids)
# The default distance limit is for a 10-2 grid, you can set a different
# limit for other grids
NHood <- function(X, Y, DistLim = 2.9){
  ListOut <- list()
  for (i in 1:length(X)){
    # Calculates distance
    Dist <- sqrt((X[i] - X)^2 + (Y[i] - Y)^2)
    # Identifies hemifield
    Hemi <- sign(Y[i])*sign(Y) > 0
    # Disconnects locations across different hemifields and 
    # further apart than limit distance
    NHs <- which(Dist < DistLim & Hemi) 
    ListOut[[i]] <- NHs[NHs != i]
  }
  return(ListOut)
}
################################################################

# Sources Drasdo and matching/reading functions for OCT data
source("DrasdoFunctions.R")
source("R_Matching_Suite/openVol_Source.R")
source("R_Matching_Suite/MatchSegmentation_Source.R")
source("R_Matching_Suite/MatchSegmentationONH_Source.R")
source("R_Matching_Suite/FoveaDetection_Source.R")
source("R_Matching_Suite/FundusMatch_Source.R")
source("R_Matching_Suite/SLO_to_CMP_Maps_Source.R")
source("R_Matching_Suite/Create_Overlay_Source.R")
source("R_Matching_Suite/Overlay_Maps_Source.R")

# Color map
load("R_Matching_Suite/Parula_Map")
# Archive of standard perimetric grids
load("Files/Standard_Grids")

#HACKS
#eventReactive elements don't output a NULL object as initial state, dependent reactive objects need
#to be managed with a try-catch statement. If not, the error will stop the execution with NO WARNING!
#The app prints a dot in the console every 10 seconds to avoid timeouts 
#(inactive sessions are disconnected, risky when running the test)


#################################################
# Controlling parameters
#################################################
TRACKING <<- TRUE

startTime <<- NULL # allocates NULL start time (to be set upon test start)
pauseDuration <<- 0 # allocates no pause time until a pause is done

# name of image file for first CMP image capture
OUTPUT_IMAGE <<- "data_Temp_Exam/Fundus_Image.jpg"

# Downscaling ratio and assumed original picture size for COMPASS
# Warning! Alignment is much more likely to fail at higher resolutions
Ratio <- 1/4
FullSizeCMP <- 1920
CMP_FoV <- 60


##### List of all .vol files in the data_Temp_Exam folder ########
filelist <- list.files(paste(getwd(),"/data_Temp_Exam/",sep = ""), 
                       pattern = ".vol", ignore.case = TRUE)
