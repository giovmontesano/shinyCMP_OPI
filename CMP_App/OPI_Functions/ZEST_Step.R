#GLOBAL VARIABLES ARE NEEDED (<<-) OTHERWISE THE observe({}) will occlude them from view
#for other processes

# Check if all states are stable
if (!all(st)){ #additional check for stability
  # Select an unstable state
  i <- which((!st))  #Selects unstopped locations from the current stage                     
  if(length(i) > 1){i <- sample(i, 1)}  # unstopped state
  
  # Step the selected state
  r <- ZEST.step(states[[i]])           
  
  # Initialize the false positive (FP) trial counter
  FPTrial <- 0
  
  #### RESPONSE VALIDATION ######
  # Assign maximum response window (RW) if no response is seen
  if (!r$resp$seen){r$resp$time <- RW}
  
  # If the response time is greater than 180, update the state
  if (r$resp$time  > 180){
    # Update the state
    states[[i]] <<- r$state                
    # Extract the probability density function (PDF) and update the mean and standard deviation
    PP <- states[[i]]$pdf                   
    SMeans[i] <<- ZEST.final(states[[i]]) 
    SDs[i] <<- sqrt(sum(((Domain - SMeans[i])^2*PP))) 
    
    
    ###### Updates neighbors if spatial enhancement is enabled ########
    if (SPATIAL) {
      ListNHs <- NHs[[i]]
      
      # Compute the likelihood function
      LF <- states[[i]]$likelihood[states[[i]]$domain == tail(states[[i]]$stimuli, n = 1),]
      if (!tail(states[[i]]$responses, n = 1)){
        LF <- 1 - LF
      }
      LF <- (LF - .5)
      
      # Update the neighboring states
      for (nh in ListNHs){
        CorrNH <- CorrList[[i]][ListNHs == nh]
        LFNH <- LF*CorrNH + 0.5
        states[[nh]]$pdf <<- states[[nh]]$pdf*LFNH
        states[[nh]]$pdf <<- states[[nh]]$pdf/sum(states[[nh]]$pdf)
        
        PP <<- states[[nh]]$pdf                    #extracts pdf
        SMeans[nh] <<- ZEST.final(states[[nh]]) #updates mean
        SDs[nh] <<- sqrt(sum(((Domain - SMeans[nh])^2*PP))) #updates SD
      }
    }
    ######     ######     ######     ######     ######     #
    
    # Add a small delay if in simulation mode
    if(SIMULATION){
      Sys.sleep(.01)
    }
    
  }else{
    # If the response time is less than or equal to 180, it's a false positive trial
    FPTrial <- 1
  }
  
  # Pause/FP trial
  if (SIMULATION){
    FPPause <- opiPresent(stimFP, tt=0, fpr=FPR, fnr=FNR)
    SeenFP <- FPPause$seen & FPPause$time <= Pause
    
    ############################################################
    #Comment to make FP projections on the ONH (default setup)
    #Not recommended, waste of time for motor repositioning.
    
    #This will project the FP at the same location (50 dB projection)
    stimFP <- list(x=locations$X[i], y=locations$Y[i], level=dbTocd(50), size=0.43, color="white",
                   duration=0, responseWindow = Pause)
    class(stimFP) <- "opiStaticStimulus"
    ############################################################
    
  }else{
    
    ############################################################
    #Comment to make FP projections on the ONH (default setup)
    #Not recommended, waste of time for motor repositioning.
    
    #This will project the FP at the same location (50 dB projection)
    stimFP <- list(x=locations$X[i], y=locations$Y[i], level=dbTocd(50), size=0.43, color="white",
                   duration=0, responseWindow = Pause)
    class(stimFP) <- "opiStaticStimulus"
    ############################################################
    
    FPPause <- opiPresent(stimFP)
    SeenFP <- FPPause$seen
  }
  if (SeenFP){FPTrial <- 1}
  
  # Update the false positive rate (FPR) estimate
  FPsC <<-  c(FPsC, FPTrial)
  TotLisTime <<- TotLisTime + (180 + Pause)
  EstFPR <<- sum(FPsC)/(TotLisTime)
  EstFPR <<- 1 - exp(-EstFPR*RW)
  
  # Check if all states are stable
  st <<- (unlist(lapply(states, ZEST.stop)))
  
  # Update the progress indicators
  SDs_Progress <- SDs
  SDs_Progress[st] <- 1.5
  
  ProgPercSD <<- unlist(lapply(states, getPres))/states[[1]]$maxPresentations
  ProgPercSD[st] <<- 1
  ProgPercSD[ProgPercSD > 1] <<- 1
  
}