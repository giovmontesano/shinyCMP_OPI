# Set the global variables
RW <<- 1500  # Response window (in milliseconds)
Pause <<- 200  # Pause duration (in milliseconds)
Domain <<- -5:40  # Stimulus domain (in decibels)

#################################### Functions #########
# Function to create a stimulus
makeStimHelper <- NULL
if (!SIMULATION) {
  # Function to create a stimulus for the real OPI device
  makeStimHelper <- function(db, n, x, y) {  # returns a function of (db, n)
    ff <- function(db, n) db + n
    
    body(ff) <- substitute(
      {s <- list(x = x, y = y, level = dbTocd(db), responseWindow = RW)
      class(s) <- "opiStaticStimulus"
      return(s)
      },
      list(x = x, y = y))
    return(ff)
  }
  
  # Create the false positive (FP) stimulus
  stimFP <- list(x = onh[1], y = onh[2], level = dbTocd(50), responseWindow = Pause)
  class(stimFP) <- "opiStaticStimulus"
} else {
  # Function to create a stimulus for the simulation
  makeStimHelper <- function(db, n, x, y) {  # returns a function of (db, n)
    ff <- function(db, n) db + n
    
    body(ff) <- substitute(
      {s <- list(x = x, y = y, level = dbTocd(db), size = 0.43, color = "white",
                 duration = 200, responseWindow = RW)
      class(s) <- "opiStaticStimulus"
      return(s)
      },
      list(x = x, y = y))
    return(ff)
  }
  
  # Create the false positive (FP) stimulus for the simulation
  stimFP <- list(x = onh[1], y = onh[2], level = dbTocd(50), size = 0.43, color = "white",
                 duration = 0, responseWindow = Pause)
  class(stimFP) <- "opiStaticStimulus"
}

# Function to create a bimodal ZEST PDF
makeBiModalPDF <- function(normalModePoint, weight, pdf.floor, domain = -5:40, normalize = TRUE) {
  # Load the ZEST curves
  ZEST_Curves <- read.table('OPI_Functions/ZEST_Curves.csv', header = TRUE, sep = ";")
  
  # Create the glaucoma PDF
  glaucomaPDF <- c(rep(0.001, -domain[1]), ZEST_Curves$Abnorm, rep(0.001, domain[length(domain) - 40]))
  
  # Create the shiftable healthy PDF
  healthyPDF <- function(normalModePoint) {
    temp <- c(rep(0.001, 50), ZEST_Curves$Norm, rep(0.001, 50))
    mode <- which.max(temp)
    return(temp[(mode - (normalModePoint) + domain[1]):(mode - (normalModePoint) + domain[length(domain)])])
  }
  
  # Combine the healthy and glaucoma PDFs
  npdf <- healthyPDF(normalModePoint)
  glaucomaPDF[domain > normalModePoint] <- 0
  cpdf <- npdf * weight + glaucomaPDF
  cpdf[which(cpdf < pdf.floor)] = pdf.floor
  if (normalize) { cpdf <- cpdf / sum(cpdf) }
  return(cpdf)
}

# Gets the total number of presentations per location
getPres <- function(state) {
  NPres <- state$numPresentations
  return(NPres)
}

###################################################################

###################################################################
#Starting states for ZEST
if (!SIMULATION){
  # Setup starting states for each location
  states <<- apply(locations, MARGIN = 1, function(loc) {
    ZESTBuild <- ZEST.start(
      domain=Domain,
      minStimulus=0,
      maxStimulus=40,
      maxPresentations=10,
      minNotSeenLimit=1,
      makeStim=makeStimHelper(db,n,loc[["X"]],loc[["Y"]]),
      stopType="S", stopValue= 1.5)
    ZESTBuild$Eccentricity <- loc[["Ecc"]]
    ZESTBuild$NormVals <- loc[["NormVal"]]
    return(ZESTBuild)
  })
}else{
  # Setup starting states for each location
  states <<- apply(locations, MARGIN = 1, function(loc) {
    ZESTBuild <- ZEST.start(
      domain=Domain,
      minStimulus=0,
      maxStimulus=40,
      makeStim=makeStimHelper(db,n,loc[["X"]], loc[["Y"]]),
      maxPresentations= 10,
      minNotSeenLimit=1,
      stopType="S", stopValue= 1.5, tt=loc[["Threshold"]], fpr=FPR, fnr=FNR)
    ZESTBuild$Eccentricity <- loc[["Ecc"]]
    ZESTBuild$NormVals <- loc[["NormVal"]]
    return(ZESTBuild)
  })
}

####### sets the starting priors
weigth <- 4
for (ll in 1:length(states)) {
  FixedStart <- states[[ll]]$NormVals  # starting dB for normal peak
  states[[ll]]$pdf <<- makeBiModalPDF(FixedStart, weigth, 0.001, Domain)
}
########################

# Calculates all starting means and SDs for future update
SMeans <<- numeric(dim(locations)[1])
SDs <<- numeric(dim(locations)[1])
for (i in 1:dim(locations)[1]) {
  PP <<- states[[i]]$pdf
  SMeans[i] <<- ZEST.final(states[[i]])
  SDs[i] <<- sqrt(sum(((Domain - SMeans[i])^2*PP)))
}

# Store false positive (FP) responses
FPsC <<- NULL
TotLisTime <<- 0
ProgPercSD <<- 0
EstFPR <<- 0

# Initialize the correlation list
CorrList <<- list()
for (ll in 1:length(locations$Threshold)) {
  CorrList[[ll]] <<- rep(Corr, length(NHs[[ll]]))
}