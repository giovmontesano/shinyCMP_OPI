# Close all open connections
closeAllConnections()

# Wrap the entire code block in a try-catch block to handle errors
Error_Report <<- try({
  if (SIMULATION) {
    # Set the output file path for the simulation
    OUTPUT_TEXT  <<- "data_Temp_Exam/PRL_ONH_Sim.txt"   # trial by trial output
    
    # Choose the "SimHensonRT" OPI device
    chooseOpi("SimHensonRT")
    
    # Initialize the OPI device, stop if it fails
    if (!is.null(opiInitialize(type="C", cap=6, rtData = RtDbUnits)))
      stop("opiInitialize failed")
    
    # Read the simulation data from the output file
    PRL_ON <<- read.table(OUTPUT_TEXT, sep = "\t")
    PRL <<- PRL_ON[PRL_ON$V1 %like% "PRL", ]
    ONH <<- PRL_ON[PRL_ON$V1 %like% "ONH", ]
    
    # Extract the PRL and ONH data from the simulation
    prl <<- as.numeric(strsplit(as.character(PRL), " ")[[1]][3:length(strsplit(as.character(PRL), " ")[[1]])])
    prl <<- prl[!is.na(prl)]
    onh <<- as.numeric(strsplit(as.character(ONH), " ")[[1]][3:length(strsplit(as.character(ONH), " ")[[1]])])
    onh <<- onh[!is.na(onh)]
    
  } else {
    # Choose the "Compass" OPI device
    chooseOpi("Compass")
  }
  
  if (!SIMULATION) {
    # Set the output file path for the real device
    OUTPUT_TEXT  <<- "data_Temp_Exam/PRL_ONH.txt"   # trial by trial output
    
    # Redirect the output to the file and the console
    sink(OUTPUT_TEXT, split=TRUE)
    
    # Initialize the OPI device with the specified IP address and port
    res <<- opiInitialise(ip = CMP_Address, port=50001)
    
    # Set the error handler to close the OPI device on error
    options(error=opiClose)
    
    # Print the error message, if any
    print(res$err)
    
    # Print the PRL and ONH data
    cat("PRL "); print(res$prl)
    cat("ONH "); print(res$onh)
    
    # Store the PRL and ONH data
    prl <<- res$prl
    onh <<- res$onh
    
    # Write the image data to the OUTPUT_IMAGE file
    f <<- file(OUTPUT_IMAGE, "wb")
    writeBin(res$im, f)
    close(f)
    
    # Set the fixation and tracking options
    opires <<- opiSetBackground(fixation=c(0,0,0))
    res <<- opiSetBackground(tracking_on=TRACKING)
  }
}, silent = TRUE)