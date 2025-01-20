
server <- function(input, output, session) {
  
  # Continuously prints dot to avoid timeout 
  # (just a hack to keep the app running and avoid timeouts)
  autoInvalidate <- reactiveTimer(10000)
  observe({
    autoInvalidate()
    cat(".")
  })
  
  # Renders different tabs in "Fundus alignment" tab,
  # This menu is rendered only if at least one OCT and one CMP image has been loaded.
  # Different aliment tabs are rendered according to what OCT files have been loaded.
  output$recAlign <- renderMenu({
    if(((OCT_Check$Check | OCT_CheckONH$Check) & 
        input$OUTPUT_IMAGE_S == OUTPUT_IMAGE)){ 
      
      #Shows Macular tab if loaded
      if (OCT_Check$Check){
        showTab(inputId = "alignTabs", target = "Macular_scan", select = TRUE)
      }else{
        hideTab(inputId = "alignTabs", target = "Macular_scan")
      }
      
      #Shows ONH tab if loaded
      if (OCT_CheckONH$Check){
        showTab(inputId = "alignTabs", target = "ONH_scan", select = TRUE)
      }else{
        hideTab(inputId = "alignTabs", target = "ONH_scan")
      }
      
      menuItem("Fundus aligment", tabName = "fundusAlign", icon = icon("eye"))
    }else{
      menuItem(NULL, tabName = NULL, icon = NULL)
    }
  })
  
  # Options to render the "Test" tab. This is active if the connection to the CMP has been
  # established and the fundus image has been loaded.
  output$recTest <- renderMenu({
    if((input$OUTPUT_IMAGE_S == OUTPUT_IMAGE)){ #
      
      # Shows Macular tab if loaded - structural calculations can only be performed if the
      # 10-2 is selected as the main grid
      if (OCT_Check$Check & input$standGrid %in% c("10-2","Upload your perimetric grid")){
        showTab(inputId = "settingTabs", target = "struc_Setting")
      }else{
        hideTab(inputId = "settingTabs", target = "struc_Setting")
      }
      
      menuItem("Test", tabName = "testTab", icon = icon("braille"))
      
    }else{
      Test_Results$states = NULL
      menuItem(NULL, tabName = NULL, icon = NULL)
    }
  })
  
  #Decides if test is a baseline or a follow-up
  IsBaseline <- reactiveValues(IsBaseline = TRUE)
  FU_Match <- reactiveValues(Match = NULL)
  
  # Swicth to follow-up mode if not a baseline
  # This lines up the new fundus image with baseline fundus image
  # It carries over the fovea/PRL location, old structural mesurements if provided
  # and uses the same test locations (mapped to the new fundus image).
  # New structural data can be imported to replace old ones.
  observeEvent(input$baseline_Test, {
    IsBaseline$IsBaseline = FALSE
    
    disable("gridCenter")
    disable("selectEye")
    disable("standGrid")
    disable("alignON")
    disable("user_grid")
    disable("Pat_ID")
    disable("age")
    disable("study")
    disable("userName")
    
    #Loads demogr data
    TempEnv <- new.env()
    load(input$baseline_Test$datapath, envir = TempEnv)
    
    updateTextInput(session, "Pat_ID", value = TempEnv$PatName)
    updateTextInput(session, "study", value = TempEnv$Study)
    updateTextInput(session, "userName", value = TempEnv$User)
    updateNumericInput(session, "age", value = round(as.numeric(TempEnv$Age + (Sys.time() - TempEnv$DateTime)/360)))
    updateRadioButtons(session, "standGrid", selected = TempEnv$Grid)
    
    enable("baselineAlign")
    enable("baselineAlignONH")
    enable("resets_Baseline")
    updateCheckboxInput(session, "plotDrasdo", value = FALSE)
    enable("plotDrasdo")
  })
  
  # Settings for a new baseline
  observeEvent(input$resets_Baseline, {
    IsBaseline$IsBaseline = TRUE
    
    enable("gridCenter")
    enable("selectEye")
    enable("standGrid")
    enable("alignON")
    enable("user_grid")
    enable("Pat_ID")
    enable("age")
    enable("study")
    enable("userName")
    
    disable("baselineAlign")
    disable("baselineAlignONH")
    
    shinyjs::reset("baseline_Test")
    disable("resets_Baseline")
  })
  
  
  #####################################################################################################
  
  ###############################################################
  ########### COMPASS connect tab ##############################
  #Reads in CMP fundus image only when connected
  CMP_FundusOr <- reactive({
    #######################################################################################
    #Reads fundus image and rescales between 0-1, excluding the 0 areas
    if (reactiveString$String == OUTPUT_IMAGE){
      ReadF <- load.image(OUTPUT_IMAGE)[,,1,1]
      ReadF <- (ReadF - min(ReadF[ReadF > 0]))/diff(range(ReadF[ReadF > 0]))
      ReadF[ReadF < 0] <- 0
      
      #Original fundus image
      CMP_FundusOr <- list()
      CMP_FundusOr$Image <- as.cimg(ReadF)
      CMP_FundusOr$FoV <- 60
      
      CMP_FundusOr
      #######################################################################################
      
    }
  })
  
  # This is the CMP fundus image used as a reference
  # Is aligned with baseline if follow-up
  CMP_Fundus <- reactive({
    if (reactiveString$String == OUTPUT_IMAGE){
      CMP_FundusOr <- CMP_FundusOr()
      #Fudus image for alignment and calculations, can be smaller for faster processing
      #Points can be scaled back but if image scaling is required it needs to be recalculated 
      #on the original image
      CMP_Fundus <- CMP_FundusOr
      CMP_Fundus$Image <- resize(CMP_Fundus$Image, 
                                 dim(CMP_Fundus$Image)[1]*Ratio, 
                                 dim(CMP_Fundus$Image)[2]*Ratio)
      
      #NA mask
      Dims <- dim(CMP_Fundus$Image)
      MaskReadF <- expand.grid(X = (1:Dims[1]) - (Dims[1]/2 - 10*Ratio), 
                               Y = (1:Dims[2]) - (Dims[2]/2 + 23*Ratio))
      MaskReadF <- matrix(sqrt(MaskReadF$X^2 + MaskReadF$Y^2), Dims[1], Dims[2]) > (910*Ratio)
      
      #NaNs the outer dark margins so they are not considered for alignment
      Source <- CMP_Fundus$Image
      Source[MaskReadF] <- NA
      CMP_Fundus$ImageNA <- as.cimg(Source)
      ##############################################################################
      
      
      ##############################################################################
      # Matches with baseline if follow-up - imports data from baseline test
      if (!IsBaseline$IsBaseline){
        ############ Loads baseline data #######################
        TempEnv <- new.env()
        load(input$baseline_Test$datapath, envir = TempEnv)
        CMP_Fundus$BaselineF <- TempEnv$CMP_Fundus$Image
        
        ###############################################
        #Aligns with baseline
        Source <- CMP_Fundus$Image
        Target <- CMP_Fundus$BaselineF
        
        #NA mask
        Dims <- dim(CMP_Fundus$Image)
        MaskReadF <- expand.grid(X = (1:Dims[1]) - (Dims[1]/2 - 10*Ratio), 
                                 Y = (1:Dims[2]) - (Dims[2]/2 + 23*Ratio))
        MaskReadF <- matrix(sqrt(MaskReadF$X^2 + MaskReadF$Y^2), Dims[1], Dims[2]) > (910*Ratio)
        
        #NaNs the outer dark margins so they are not considered for alignment
        Source[MaskReadF] <- NA
        Source <- as.cimg(Source)
        Target[MaskReadF] <- NA
        Target <- as.cimg(Target)
        
        CMP_Fundus$MatchFU <- niftyreg(Source, Target, scope = "affine", internal = FALSE) 
        
        #################################################################
        
        updateRadioButtons(session, "gridCenter", #choices = c("Fovea","PRL"),
                           selected = TempEnv$GridCenter)
        updateRadioButtons(session, "selectEye",
                           selected = TempEnv$Eye)
        
        ##################################################################
        
        prlPix <- applyTransform(reverse(CMP_Fundus$MatchFU), (TempEnv$prl*32 + 960)*Ratio)
        onhPix <- applyTransform(reverse(CMP_Fundus$MatchFU), (TempEnv$onh*32 + 960)*Ratio)
        
        prl <<- (prlPix/Ratio - 960)/32
        onh <<- (onhPix/Ratio - 960)/32
      }
      CMP_Fundus
    }else{
      CMP_Fundus <- list()
      CMP_Fundus$Image <- as.cimg(matrix(1, 1, 1))
      CMP_Fundus$FoV <- 60
      CMP_Fundus
    }
  })
  
  # Renders CMP fundus image
  output$previewCMP <- renderPlot({
    CMP_Fundus <- CMP_Fundus()
    par(mar = c(0, 0, 0, 0))
    plot(CMP_Fundus$Image, axes = FALSE)
  })
  
  # Disables the IP address text-box if in simulation mode
  observeEvent(input$isSimulation, {
    if(input$isSimulation){
      disable("CMP_Address")
    }else{
      enable("CMP_Address")
    }
    
  })
  
  # Makes it reactive so that it updates every time 
  # connect is pressed (check based on string)
  reactiveString <- reactiveValues(
    String = "NULL"
  )
  
  # Connects to CMP and loads baseline image
  observeEvent(input$buttonConnect, {
    Pat_ID <<- input$Pat_ID
    Test_Results$states = NULL
    
    if (gsub("\\s", "", Pat_ID)   != ""){
      CMP_Address <<- input$CMP_Address
      SIMULATION <<- sum(input$isSimulation == 1)
      Aborted <<- FALSE
      reactiveString$String = "NULL"
      
      source("OPI_Functions/devConnect.R")
      
      if (is.null(attr(Error_Report, "condition"))){
        reactiveString$String = OUTPUT_IMAGE
        updateTextInput(session, inputId = "OUTPUT_IMAGE_S", value = OUTPUT_IMAGE)
        updateTextInput(session, inputId = "isAligned", value = "Not aligned!")
        
        #Guesses eye based on the position of the ONH
        if (onh[1] < 0){
          updateSelectInput(session, inputId = "selectEye", selected = "Left")
        }else{
          updateSelectInput(session, inputId = "selectEye", selected = "Right")
        }
        
      }else{
        updateTextInput(session, inputId = "OUTPUT_IMAGE_S", value = "")
        showModal(modalDialog(
          title = "Impossible to connect to the Compass",
          "Make sure that the connection is running, that the Compass address is correct, that the first 
          image has been acquired, that the ONH has been located and that PRL has been intialized ",
          easyClose = TRUE))
      }
    }else{
      updateTextInput(session, inputId = "OUTPUT_IMAGE_S", value = "")
      showModal(modalDialog(
        title = "Patient ID missing",
        "Please input a Subject ID and check subject's data",
        easyClose = TRUE))
    }
  })
  
  #####################################################################################################
  #####################################################################################################
  
  
  
  
  
  
  
  
  ################################## OCT Data #####################################
  # Macular scans ##################################
  #####################################################################################################
  
  ##################################################
  #Loads Spectralis .vol macular data when pressed
  OCT_Data <- eventReactive(input$buttonStruct, {
    FoveaXY_SLO_Store$XY = NULL
    disable("buttonSegm")
    disable("buttonFovea")
    if(IsBaseline$IsBaseline){
      updateRadioButtons(session, "standGrid", selected = "None")
    }
    updateCheckboxInput(session, "plotDrasdo", value = FALSE)
    disable("plotDrasdo")
    OCT_Check$Check = FALSE
    CMP_MapCheck$Check = FALSE
    updateSelectInput(session, "segmList", choices = "None")
    updateSelectInput(session, "segmListCMPTest", choices = "None")
    OCT_Data <- FixScanSize(openVol(input$octLIst, quickFovea = FALSE, quickSegm = TRUE))
    if (OCT_Data$Header$ScanPattern != 3){
      showModal(modalDialog(
        title = "Wrong scan pattern!",
        "Not a volume scan",
        easyClose = TRUE))
      NULL
    }else{
      enable("buttonSegm")
      enable("buttonFovea")
      OCT_Data
    }
  })
  
  # Spectralis SLO image
  SLO <- reactive({
    OCT_Data <- OCT_Data()
    if (!is.null(OCT_Data)){
      #Reads SLO image and rescales between 0-1, excluding the 0 areas
      SLO <- list()
      SLO$Image <- OCT_Data$SLO
      SLO$Image <- (SLO$Image - min(SLO$Image[SLO$Image > 0]))/diff(range(SLO$Image[SLO$Image > 0]))
      SLO$Image[SLO$Image < 0] <- 0
      SLO$Image <- as.cimg(SLO$Image)
      SLO$FoV <- OCT_Data$Header$FieldSizeSlo
      
      ScSizeSLO <- (FullSizeCMP/CMP_FoV)*SLO$FoV*Ratio
      SLO$ImageSc <- resize(SLO$Image, ScSizeSLO, ScSizeSLO)
      SLO
    }else{
      NULL
    }
  })
  
  
  # Matches  structural maps to SLO and detects fovea
  Structural_Maps <- eventReactive(input$buttonSegm, {
    OCT_Data <- OCT_Data()
    SLO <- SLO()
    
    Layer_List <- c("WRT","GCL","NFL")
    updateSelectInput(session, "segmList", choices = c(Layer_List, "None"))
    updateSelectInput(session, "segmListStrat", choices = c(Layer_List, "None"), selected = "GCL")
    
    SEG_List <- MatchSegmentationVolMac(OCT_Data$BScanHeader, OCT_Data$Header, Layer_List)
    updateTextInput(session, inputId = "isAligned", value = "Not aligned!")
    
    #Calculates overlays
    Overlays <- list()
    for (i in 1:length(Layer_List)){
      #Scaling limits for colormap (microns)
      
      LimitsTh <- switch(Layer_List[i], 
                         "WRT" = c(200, 350),
                         "NFL" = c(0, 120),
                         "GCL" = c(0, 65),
                         "IPL" = c(0, 65),
                         "OUTER" = c(120, 200 - 120))
      
      Overlays[[i]] <- as.cimg(Overlay_Maps(SLO$Image, SEG_List[[Layer_List[i]]], LimMap = LimitsTh, Replace = TRUE))
    }
    
    names(Overlays) <- Layer_List
    disable("buttonSegm")
    list(SEG_List = SEG_List, Overlays = Overlays)
  })
  
  #Stores manual centering fovea
  CenterFov_Search <- reactiveValues(CX = 385, CY = 385)
  observeEvent(input$plot_click_fovea, {
    CX <- input$plot_click_fovea$x 
    CY <- input$plot_click_fovea$y
    
    if (CX < 161){CX <- 161}
    if (CY < 161){CY <- 161}
    if (CX > 608){CX <- 608}
    if (CY > 608){CY <- 608}
    
    CenterFov_Search$CX = CX
    CenterFov_Search$CY = CY
    enable("buttonFovea")
  })
  
  #Finds fovea
  FoveaXY_SLO_Store <- reactiveValues(XY = NULL)
  FoveaXY_SLO <- reactive({FoveaXY_SLO_Store$XY})
  observeEvent(input$buttonFovea, {
    Structural_Maps <- Structural_Maps()
    OCT_Check$Check = TRUE
    FoveaXY_SLO <- FoveaDetection(WRT = Structural_Maps$SEG_List$WRT, CX = CenterFov_Search$CX, CY = CenterFov_Search$CY)
    disable("buttonFovea")
    FoveaXY_SLO_Store$XY = FoveaXY_SLO
  })
  
  
  #Preview macular structural data
  output$previewOCT <- renderPlot({
    SLO <- SLO()
    if (!is.null(SLO)){
      Structural_Maps <- try({Structural_Maps <- Structural_Maps()}, silent = TRUE)
      if(input$segmList != "None" & is.null(attr(Structural_Maps, "condition"))){
        par(mar = c(0, 0, 0, 0))
        Image <- Structural_Maps$Overlays[[input$segmList]]
        plot(Image, axes = FALSE)
        
        RectX <- CenterFov_Search$CX + c(-160, 160, 160, -160, -160)
        RectY <- CenterFov_Search$CY + c(160, 160, -160, -160, 160)
        lines(RectX, RectY, col = "red", lwd = 3)
        
      }else{
        par(mar = c(0, 0, 0, 0))
        plot(SLO$Image, axes = FALSE)
        
        RectX <- CenterFov_Search$CX + c(-160, 160, 160, -160, -160)
        RectY <- CenterFov_Search$CY + c(160, 160, -160, -160, 160)
        lines(RectX, RectY, col = "red", lwd = 3)
      }
      
      #Checks if the fovea has been detected
      FoveaXY_SLO <- try({FoveaXY_SLO <- FoveaXY_SLO()}, silent = TRUE)
      if(is.null(attr(FoveaXY_SLO, "condition"))){
        points(FoveaXY_SLO[1], FoveaXY_SLO[2], col = 2, pch = 10, cex = 3)
      }
    }
  })
  
  # Colormap for thickness
  output$previewOCTcolorbar <- renderPlot({
    Segm <- input$segmList
    if(Segm != "None"){
      LimitsTh <- switch(Segm, 
                         "WRT" = c(200, 350),
                         "NFL" = c(0, 120),
                         "GCL" = c(0, 65),
                         "IPL" = c(0, 65),
                         "OUTER" = c(120, 200))
      par(mar = c(1, 3, 3, 2))
      color.barv(colorRampPalette(Parula_Map$Hex)(100), LimitsTh[1], LimitsTh[2], 
                 title = expression(atop("Thickness ", paste("(",mu,"m)", sep =""))))
    }else{
      NULL
    }
  })
  
  
  
  
  #####################################################################################################
  # ONH scans ##################################
  #####################################################################################################
  
  #Loads Spectralis peripapillary scan when pressed
  OCT_DataONH <- eventReactive(input$buttonStructONH, {
    disable("buttonSegm")
    
    OCT_CheckONH$Check = FALSE 
    
    updateSelectInput(session, "segmListONH", choices = "None")
    OCT_Data <- FixScanSize(openVol(input$octLIstONH, quickFovea = FALSE, quickSegm = TRUE))
    if (!(OCT_Data$Header$ScanPattern  %in% c(2,6))){
      showModal(modalDialog(
        title = "Wrong scan pattern!",
        "Not a peripapillary  scan",
        easyClose = TRUE))
      NULL
    }else{
      enable("buttonSegmONH")
      OCT_Data
    }
  })
  
  # Loads SLO Spectralis image for the peripillary scans
  SLOONH <- reactive({
    OCT_Data <- OCT_DataONH()
    if (!is.null(OCT_Data)){
      #Reads SLO image and rescales between 0-1, excluding the 0 areas
      SLO <- list()
      SLO$Image <- OCT_Data$SLO
      SLO$Image <- (SLO$Image - min(SLO$Image[SLO$Image > 0]))/diff(range(SLO$Image[SLO$Image > 0]))
      SLO$Image[SLO$Image < 0] <- 0
      SLO$Image <- as.cimg(SLO$Image)
      SLO$FoV <- OCT_Data$Header$FieldSizeSlo
      
      ScSizeSLO <- (FullSizeCMP/CMP_FoV)*SLO$FoV*Ratio
      SLO$ImageSc <- resize(SLO$Image, ScSizeSLO, ScSizeSLO)
      SLO
    }else{
      NULL
    }
  })
  
  
  # Matches  structural peripapillary thickness to SLO 
  Structural_MapsONH <- eventReactive(input$buttonSegmONH, {
    OCT_Data <- OCT_DataONH()
    OCT_CheckONH$Check = TRUE #!!!!!
    disable("buttonSegmONH")
    Layer_List <- c("NFL")
    updateSelectInput(session, "segmListONH", choices = c(Layer_List, "None"))
    updateTextInput(session, inputId = "isAlignedONH", value = "Not aligned!")
    
    NFL_Seg <- MatchSegmentationVolONH(OCT_Data$BScanHeader, OCT_Data$Header)
    
    ColMap <- data.frame(c = colorRampPalette(Parula_Map$Hex)(100),
                         Th = seq(0, 200, length.out = 100))
    CInd <- unlist(lapply(X = NFL_Seg$NFL, 
                          FUN = function(x, V) which.min(abs(x-V)), 
                          V = ColMap$Th))
    
    NFL_Seg$Col <- colorRampPalette(Parula_Map$Hex)(100)[CInd]
    
    NFL_Seg
  })
  
  # Renders matched peripapillary thickness map
  output$previewOCTONH <- renderPlot({
    SLO <- SLOONH()
    if (!is.null(SLO)){
      NFL_Seg <- try({NFL_Seg <- Structural_MapsONH()}, silent = TRUE)
      if(input$segmListONH != "None" & is.null(attr(Structural_MapsONH, "condition"))){
        Image <- SLO$Image
        par(mar = c(0, 0, 0, 0))
        plot(Image, axes = FALSE)
        
        points(NFL_Seg$x, NFL_Seg$y, pch = 20, cex = NFL_Seg$NFL/100,
               col = NFL_Seg$Col)
        
      }else{
        par(mar = c(0, 0, 0, 0))
        plot(SLO$Image, axes = FALSE)
      }
    }
  })
  
  # Colormap for thickness
  output$previewOCTcolorbarONH <- renderPlot({
    Segm <- input$segmListONH
    if(Segm != "None"){
      LimitsTh <- c(0,200)
      par(mar = c(1, 3, 3, 2))
      color.barv(colorRampPalette(Parula_Map$Hex)(100), LimitsTh[1], LimitsTh[2], 
                 title = expression(atop("Thickness ", paste("(",mu,"m)", sep =""))))
    }else{
      NULL
    }
  })
  
  
  
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################
  
  # Reactive values for checks
  
  OCT_Check <- reactiveValues(Check = FALSE)
  CMP_MapCheck <- reactiveValues(Check = FALSE)
  CMP_Check <- reactiveValues(Check = FALSE)
  OCT_CheckONH <- reactiveValues(Check = FALSE)
  
  
  ###############################################################
  ###############################################################
  ############### Fundus alignment tab ##########################
  
  #####################################################################
  #################### Align macular scan #############################
  # Shows alignment process for the macular scans
  output$alingnedPlot <- renderPlot({
    if (input$move){
      #Manual shifts
      plot(MovingPlot(), axes = FALSE, rescale = FALSE)
    }else{
      #Final alignment/static manual
      FixPlot <- FixPlot()
      plot(FixPlot[[2]], axes = FALSE, rescale = FALSE)
    }
  })
  
  # This will trigger live move (just dragging the mouse will move the picture)
  observeEvent(input$move,{
    if(input$move){
      updateCheckboxInput(session, inputId = "showAlign", value = FALSE)
      updateTextInput(session, inputId = "isAligned", value = "Not aligned!")
      disable("showAlign")
      disable("align")
    }else{
      enable("align")
    }
  })
  
  
  # Stops moving the fundus image when clicked if in live move
  observeEvent(input$plot_click, {
    if (input$move){
      updateCheckboxInput(session, "move", value = FALSE)
    }
    Center$XY = c(input$plot_click$x, input$plot_click$y)
    updateTextInput(session, "isAligned", value = "Not aligned!")
    disable("showAlign")
    enable("align")
  })
  
  #Reads manual centering
  Center <- reactiveValues(XY = NULL)
  
  # Alignes SLO when "Align" is pressed
  SLO_Aligned <- eventReactive(input$align,{
    CMP_MapCheck$Check = FALSE
    updateSelectInput(session, "segmListCMP", choices = c("None"))
    
    if (input$baselineAlign & !IsBaseline$IsBaseline){
      
      CMP_Fundus <- CMP_Fundus()
      
      TempEnv <- new.env()
      load(input$baseline_Test$datapath, envir = TempEnv)
      Match <- TempEnv$SLO_Aligned
      
      TranslationM <- diag(1, 4, 4)
      TranslationM[1:2,4] <- -c((Match$CropTX*(1 - Match$tol)) - Match$Offset)
      Match$forwardTransforms[[1]] <- asAffine((Match$forwardTransforms[[1]] %*% TranslationM)
                                               %*% CMP_Fundus$MatchFU$reverseTransforms[[1]] %*% inv(TranslationM))
      Match$reverseTransforms[[1]] <- invertAffine(Match$forwardTransforms[[1]])
      Match$image <- as.cimg(applyTransform(forward(Match), Match$sourceR))
      
    }else{
      SLO <- SLO()
      Center <- Center$XY
      CMP_Fundus <- CMP_Fundus()
      if (!is.null(Center)){
        Center <- round(Center)
      }else{
        Center <- dim(CMP_Fundus$Image)[1:2]/2 #Default half size of Fundus Image
      }
      CropSize <- (dim(CMP_Fundus$Image)[1] - SLO$FoV*dim(CMP_Fundus$Image)[1]/CMP_Fundus$FoV)/2
      Match <- FundusMatch(SLO$Image, CMP_Fundus$ImageNA, CropSize, "affine", Center)  
    }
    Match
  })
  
  # Shows the final alignment 
  observeEvent(input$align, {
    enable("showAlign")
    updateCheckboxInput(session, "showAlign", value = TRUE)
    updateTextInput(session, "isAligned", value = "Aligned!")
  })
  
  # Frozen overlay for alignment, after clicking on on the CMP image (final position before
  # alignment)
  FixPlot <- reactive({
    SLO <- SLO()
    CMP_Fundus <- CMP_Fundus()
    
    #Reads the center
    Center <- Center$XY
    
    CenterX <- Center[1]
    CenterY <- Center[2]
    
    if (is.null(CenterX) | input$baselineAlign){
      CenterX <- dim(CMP_Fundus$Image)[1]/2
      CenterY <- dim(CMP_Fundus$Image)[1]/2
    }
    par(mar = c(0, 0, 0, 0))
    
    CropSize <- (dim(CMP_Fundus$Image)[1] - SLO$FoV*dim(CMP_Fundus$Image)[1]/CMP_Fundus$FoV)/2
    #Recenters
    target <- imshift(CMP_Fundus$Image, dim(CMP_Fundus$Image)[1]/2 - CenterX, dim(CMP_Fundus$Image)[2]/2 - CenterY)
    #Crops central location (guess similarity)
    targetC <- crop.borders(target, nx = CropSize, ny = CropSize)
    
    
    if (input$showAlign & input$isAligned == "Aligned!"){
      
      Match <- SLO_Aligned()
      
      #Composite with Aligned SLO
      SLO_Al <- resize(as.cimg(Match$image), dim(Match$targetC)[1], dim(Match$targetC)[2])
      SLO_Al <- crop.borders(SLO_Al, nx = round(CropSize*(Match$tol)), 
                             ny = round(CropSize*(Match$tol)))
      SLO_Al <- (SLO_Al - min(SLO_Al[SLO_Al > 0], na.rm = TRUE))/
        diff(range(SLO_Al[SLO_Al > 0], na.rm = TRUE))
      SLO_Al[SLO_Al < 0] <- 0
      
      #Multiplies by .5 for scaled sum
      OutputIm <- SLO_Al*.5 + targetC*.5
      
    }else{
      #Composite with Raw SLO
      #Multiplies by .5 for scaled sum
      OutputIm <- SLO$ImageSc*.5 + targetC*.5
    }
    
    FixPlot <- list(Center = c(input$plot_click$x, input$plot_click$y),
                    Plot = OutputIm)
  })
  
  # Moving overlay for alignment
  MovingPlot <- reactive({
    SLO <- SLO()
    CMP_Fundus <- CMP_Fundus()
    
    CenterX <- input$plot_hover$x
    CenterY <- input$plot_hover$y
    
    if (is.null(CenterX)){
      CenterX <- dim(CMP_Fundus$Image)[1]/2
      CenterY <- dim(CMP_Fundus$Image)[1]/2
    }
    par(mar = c(0, 0, 0, 0))
    
    #Recenters
    target <- imshift(CMP_Fundus$Image, dim(CMP_Fundus$Image)[1]/2 - CenterX, dim(CMP_Fundus$Image)[2]/2 - CenterY)
    #Crops central location (guess similarity)
    targetC <- crop.borders(target, nx = dim(CMP_Fundus$Image)[1]/4, ny = dim(CMP_Fundus$Image)[1]/4)
    
    #Multiplies by .5 for scaled sum
    OutputIm <- SLO$ImageSc*.5 + targetC*.5
  })
  
  # Calculates matched structural maps for CMP - uses the affine transformation
  # calculated during alignment to map OCT structural maps onto CMP fundus image
  Structural_MapsCMP <- eventReactive(input$matchMapsCMP, {
    #Loads OCT and Fundus images
    SLO <- SLO()
    CMP_Fundus <- CMP_Fundus()
    Match <- SLO_Aligned()
    Structural_Maps <- Structural_Maps()
    
    SEG_List <- Structural_Maps$SEG_List
    
    Layer_List <- names(Structural_Maps$Overlays)
    updateSelectInput(session, "segmListCMP", choices = c(Layer_List, "SLO", "None"))
    updateSelectInput(session, "segmListCMPTest", choices = c(Layer_List, "SLO", "None"))
    
    #Calculates overlays
    Overlays <- list()
    for (i in 1:length(Layer_List)){
      #Scaling limits for colormap (microns)
      
      LimitsTh <- switch(Layer_List[i], 
                         "WRT" = c(200, 350),
                         "NFL" = c(0, 120),
                         "GCL" = c(0, 65),
                         "IPL" = c(0, 65),
                         "OUTER" = c(120, 200))
      
      Overlays[[i]] <- as.cimg(Create_Overlay(SEG_List[[Layer_List[i]]], CMP_Fundus$Image, Match, LimMap = LimitsTh, Replace = TRUE))
    }
    names(Overlays) <- Layer_List
    Overlays$SLO <- as.cimg(Create_Overlay(SLO$Image, CMP_Fundus$Image, Match, Replace = TRUE, Pedestal = FALSE))
    
    CMP_MapCheck$Check = TRUE
    enable("plotDrasdo")
    
    list(SEG_List = SEG_List, Overlays = Overlays)
  })
  
  # Allows mapping of structural data onto CMP fundus image if aligned
  observe({
    if (input$isAligned == "Aligned!"){
      enable("matchMapsCMP")
    }else{
      CMP_MapCheck$Check = FALSE
      updateProgressBar(
        session = session,
        id = "pb1",
        value = 0, total = 100,
        title = "   "
      )
      disable("matchMapsCMP")
      updateSelectInput(session, "segmListCMP", choices = c("None"))
    }
  })
  
  # Fovea translation on CMP fundus image (scaled, not for original image)
  FoveaXY_CMP <- reactive({
    if (IsBaseline$IsBaseline){
      if (CMP_MapCheck$Check){
        updateRadioButtons(session, "gridCenter", 
                           choices = c("Fovea","PRL"),
                           selected = NULL)
        
        FoveaXY_SLO <- FoveaXY_SLO()
        Match <- SLO_Aligned()
        FoveaXY_CMP <- (Match$CropTX*(1 - Match$tol)) +
          applyTransform(forward(Match), FoveaXY_SLO*Match$ScF) - Match$Offset
      }else{
        updateRadioButtons(session, "gridCenter", 
                           choices = c("PRL"),
                           selected = NULL)
        c(NA, NA)
      }
    }else{
      CMP_Fundus <- CMP_Fundus()
      TempEnv <- new.env()
      load(input$baseline_Test$datapath, envir = TempEnv)
      
      FoveaXY_CMP <- applyTransform(reverse(CMP_Fundus$MatchFU), TempEnv$FoveaXY_CMP*TempEnv$Ratio)
      
      if (is.na(FoveaXY_CMP[1])){
        if (CMP_MapCheck$Check){
          updateRadioButtons(session, "gridCenter", 
                             choices = c("Fovea","PRL"),
                             selected = NULL)
          
          FoveaXY_SLO <- FoveaXY_SLO()
          Match <- SLO_Aligned()
          FoveaXY_CMP <- (Match$CropTX*(1 - Match$tol)) +
            applyTransform(forward(Match), FoveaXY_SLO*Match$ScF) - Match$Offset
        }else{
          updateRadioButtons(session, "gridCenter", 
                             choices = c("PRL"),
                             selected = NULL)
          c(NA, NA)
        }
      }else{
        FoveaXY_CMP 
      }
    }
  })
  
  
  # Coordinates for the ONH on CMP - based on OCT
  ONH_Coords <- eventReactive(input$useOCTONH,{
    if (IsBaseline$IsBaseline){
      if (input$useOCTONH){
        NFLSeg <- Structural_MapsCMPONH()
        ONH_Coords <- c(NFLSeg$NFL_Seg$ONH_X, NFLSeg$NFL_Seg$ONH_Y)
      }else{
        ONH_Coords <- (32*onh + 960)*Ratio
      }
    }else{
      CMP_Fundus <- CMP_Fundus()
      TempEnv <- new.env()
      load(input$baseline_Test$datapath, envir = TempEnv)
      disable("useOCTONH")
      ONH_Coords <- applyTransform(reverse(CMP_Fundus$MatchFU), TempEnv$ONH_Coords)
    }
  })
  
  #Shows fundus image + structural overlays if available
  output$fundusImage <- renderPlot({
    CMP_Fundus <- CMP_Fundus()
    Structural_MapsCMP <- try({Structural_MapsCMP()}, silent = TRUE)
    
    if(input$segmListCMP != "None" & is.null(attr(Structural_MapsCMP, "condition"))){
      FoveaXY_CMP <- FoveaXY_CMP()
      par(mar = c(0, 0, 0, 0))
      plot(Structural_MapsCMP$Overlays[[input$segmListCMP]], axes = FALSE)
      points(FoveaXY_CMP[1], FoveaXY_CMP[2], col = 2, pch = 10, cex = 3)
    }else{
      par(mar = c(0, 0, 0, 0))
      plot(CMP_Fundus$Image, axes = FALSE)
    }
  })
  
  # Colorbar for structural maps in the CMP fundus image
  output$fundusImagecolorbar <- renderPlot({
    Segm <- input$segmListCMP
    if(Segm != "None" & Segm != "SLO"){
      LimitsTh <- switch(Segm, 
                         "WRT" = c(200, 350),
                         "NFL" = c(0, 120),
                         "GCL" = c(0, 65),
                         "IPL" = c(0, 65),
                         "OUTER" = c(120, 200))
      par(mar=c(4,2,0,2))
      color.barh(colorRampPalette(Parula_Map$Hex)(100), LimitsTh[1], LimitsTh[2], 
                 title = expression(paste("Thickness ","(",mu,"m)")))
    }else{
      NULL
    }
  })
  
  ###############################################################################
  #### ONH alignment #######################
  # Shows the alignment process for ONH SLO
  output$alingnedPlotONH <- renderPlot({
    if (input$moveONH){
      #Manual shifts
      plot(MovingPlotONH(), axes = FALSE, rescale = FALSE)
    }else{
      #Final alignment/static manual
      FixPlot <- FixPlotONH()
      plot(FixPlot[[2]], axes = FALSE, rescale = FALSE)
    }
  })
  
  # Shows the alignment process for ONH SLO and CMP fundus image
  observeEvent(input$moveONH,{
    if(input$moveONH){
      updateCheckboxInput(session, inputId = "showAlignONH", value = FALSE)
      updateTextInput(session, inputId = "isAlignedONH", value = "Not aligned!")
      disable("showAlignONH")
      disable("alignONH")
    }else{
      enable("alignONH")
    }
  })
  
  # Stops moving the fundus image when clicked if in live move
  observeEvent(input$plot_clickONH, {
    if (input$moveONH){
      updateCheckboxInput(session, "moveONH", value = FALSE)
    }
    CenterONH$XY = c(input$plot_clickONH$x, input$plot_clickONH$y)
    CenterONHO <<- c(CenterONH$XY)
    updateTextInput(session, "isAlignedONH", value = "Not aligned!")
    disable("showAlignONH")
    enable("alignONH")
  })
  
  # Reads manual centering from the initial setting up of the CMP
  CenterONH <- reactiveValues(XY = NULL)
  
  # Alignes SLO when "Align" is pressed
  SLO_AlignedONH <- eventReactive(input$alignONH,{
    updateSelectInput(session, "segmListCMPONH", choices = c("None"))
    
    if (input$baselineAlignONH & !IsBaseline$IsBaseline){
      
      CMP_Fundus <- CMP_Fundus()
      
      TempEnv <- new.env()
      load(input$baseline_Test$datapath, envir = TempEnv)
      Match <- TempEnv$SLO_AlignedONH
      
      TranslationM <- diag(1, 4, 4)
      TranslationM[1:2,4] <- -c((Match$CropTX*(1 - Match$tol)) - Match$Offset)
      Match$forwardTransforms[[1]] <- asAffine((Match$forwardTransforms[[1]] %*% TranslationM)
                                               %*% CMP_Fundus$MatchFU$reverseTransforms[[1]] %*% inv(TranslationM))
      Match$image <- as.cimg(applyTransform(forward(Match), Match$sourceR))
      
    }else{
      
      NFL_Seg <- Structural_MapsONH()
      SLO <- SLOONH()
      Center <- CenterONH$XY
      CMP_Fundus <- CMP_Fundus()
      if (!is.null(Center)){
        Center <- round(Center)
      }else{
        Center[1] <- (32*(onh[1]) + 960)*Ratio
        Center[2] <- (32*(onh[2]) + 960)*Ratio
      }
      ONH_OffSpec <- c(NFL_Seg$ONH_X - 384, NFL_Seg$ONH_Y - 384)*30/768*32*Ratio
      Center[1] <- Center[1] - ONH_OffSpec[1]
      Center[2] <- Center[2] - ONH_OffSpec[2]
      
      CropSize <- (dim(CMP_Fundus$Image)[1] - SLO$FoV*dim(CMP_Fundus$Image)[1]/CMP_Fundus$FoV)/2
      Match <- FundusMatch(SLO$Image, CMP_Fundus$ImageNA, CropSize, "affine", Center)
    }
    Match
  })
  
  # Shows aligned ONH SLO and CMP fundus image
  observeEvent(input$alignONH, {
    enable("showAlignONH")
    updateCheckboxInput(session, "showAlignONH", value = TRUE)
    updateTextInput(session, "isAlignedONH", value = "Aligned!")
  })
  
  # Frozen overlay for alignment, after clicking on on the CMP image (final position before
  # alignment)
  FixPlotONH <- reactive({
    SLO <- SLOONH()
    CMP_Fundus <- CMP_Fundus()
    
    #Reads the center
    Center <- CenterONH$XY
    
    CenterX <- Center[1]
    CenterY <- Center[2]
    
    NFL_Seg <- Structural_MapsONH()
    
    if (is.null(CenterX)){
      Center[1] <- (32*(onh[1]) + 960)*Ratio
      Center[2] <- (32*(onh[2]) + 960)*Ratio
    }
    ONH_OffSpec <- c(NFL_Seg$ONH_X - 384, NFL_Seg$ONH_Y - 384)*30/768*32*Ratio
    CenterX <- Center[1] - ONH_OffSpec[1]
    CenterY <- Center[2] - ONH_OffSpec[2]
    
    par(mar = c(0, 0, 0, 0))
    
    CropSize <- (dim(CMP_Fundus$Image)[1] - SLO$FoV*dim(CMP_Fundus$Image)[1]/CMP_Fundus$FoV)/2
    #Recenters
    target <- imshift(CMP_Fundus$Image, dim(CMP_Fundus$Image)[1]/2 - CenterX, dim(CMP_Fundus$Image)[2]/2 - CenterY)
    #Crops central location (guess similarity)
    targetC <- crop.borders(target, nx = CropSize, ny = CropSize)
    
    
    if (input$showAlignONH & input$isAlignedONH == "Aligned!"){
      
      Match <- SLO_AlignedONH()
      
      #Composite with Aligned SLO
      SLO_Al <- resize(as.cimg(Match$image), dim(Match$targetC)[1], dim(Match$targetC)[2])
      SLO_Al <- crop.borders(SLO_Al, nx = round(CropSize*(Match$tol)), 
                             ny = round(CropSize*(Match$tol)))
      SLO_Al <- (SLO_Al - min(SLO_Al[SLO_Al > 0], na.rm = TRUE))/
        diff(range(SLO_Al[SLO_Al > 0], na.rm = TRUE))
      SLO_Al[SLO_Al < 0] <- 0
      
      #Multiplies by .5 for scaled sum
      OutputIm <- SLO_Al*.5 + targetC*.5
      
    }else{
      #Composite with Raw SLO
      #Multiplies by .5 for scaled sum
      OutputIm <- SLO$ImageSc*.5 + targetC*.5
    }
    
    FixPlot <- list(Center = c(input$plot_clickONH$x, input$plot_clickONH$y),
                    Plot = OutputIm)
  })
  
  # Move SLO overlay for alignment
  MovingPlotONH <- reactive({
    SLO <- SLOONH()
    CMP_Fundus <- CMP_Fundus()
    
    CenterX <- input$plot_hoverONH$x
    CenterY <- input$plot_hoverONH$y
    
    NFL_Seg <- Structural_MapsONH()
    if (is.null(CenterX)){
      ONH_X <- (32*(onh[1]) + 960)*Ratio
      ONH_Y <- (32*(onh[2]) + 960)*Ratio
      
      CenterX <- ONH_X
      CenterY <- ONH_Y
    }
    
    ONH_OffSpec <- c(NFL_Seg$ONH_X - 384, NFL_Seg$ONH_Y - 384)*30/768*32*Ratio
    CenterX <- CenterX - ONH_OffSpec[1]
    CenterY <- CenterY - ONH_OffSpec[2]
    
    par(mar = c(0, 0, 0, 0))
    
    #Recenters
    target <- imshift(CMP_Fundus$Image, dim(CMP_Fundus$Image)[1]/2 - CenterX, dim(CMP_Fundus$Image)[2]/2 - CenterY)
    #Crops central location (guess similarity)
    targetC <- crop.borders(target, nx = dim(CMP_Fundus$Image)[1]/4, ny = dim(CMP_Fundus$Image)[1]/4)
    
    #Multiplies by .5 for scaled sum
    OutputIm <- SLO$ImageSc*.5 + targetC*.5
  })
  
  # Calculates matched structural maps for CMP - uses the affine transformation
  # calculated during alignment to map OCT structural maps onto CMP fundus image
  Structural_MapsCMPONH <- eventReactive(input$matchMapsCMPONH, {
    #Loads OCT and Fundus images
    SLO <- SLOONH()
    CMP_Fundus <- CMP_Fundus()
    Match <- SLO_AlignedONH()
    NFL_Seg <- Structural_MapsONH()
    
    Layer_List <- c("NFL")
    updateSelectInput(session, "segmListCMPONH", choices = c(Layer_List, "SLO", "None"))
    updateSelectInput(session, "segmListCMPTestONH", choices = c(Layer_List, "SLO", "None"))
    
    #Calculates overlays
    Overlays <- list()
    Overlays$SLO <- as.cimg(Create_Overlay(SLO$Image, CMP_Fundus$Image, Match, Replace = TRUE, Pedestal = FALSE))
    
    
    NFL_SegXY_CMP <- cbind(NFL_Seg$x, NFL_Seg$y)
    NFL_SegXY_CMP <- (Match$CropTX*(1 - Match$tol)) + 
      applyTransform(forward(Match), NFL_SegXY_CMP*Match$ScF)
    NFL_SegXY_CMP[,1] <- NFL_SegXY_CMP[,1] - Match$Offset[1]
    NFL_SegXY_CMP[,2] <- NFL_SegXY_CMP[,2] - Match$Offset[2]
    
    NFL_Seg$x <- NFL_SegXY_CMP[,1]
    NFL_Seg$y <- NFL_SegXY_CMP[,2]
    
    ONH_SegXY_CMP <- cbind(NFL_Seg$ONH_X, NFL_Seg$ONH_Y)
    ONH_SegXY_CMP <- (Match$CropTX*(1 - Match$tol)) + 
      applyTransform(forward(Match), ONH_SegXY_CMP*Match$ScF)
    ONH_SegXY_CMP[1] <- ONH_SegXY_CMP[1] - Match$Offset[1]
    ONH_SegXY_CMP[2] <- ONH_SegXY_CMP[2] - Match$Offset[2]
    
    NFL_Seg$ONH_X <- ONH_SegXY_CMP[1]
    NFL_Seg$ONH_Y <- ONH_SegXY_CMP[2]
    
    if (IsBaseline$IsBaseline){
      enable("useOCTONH") 
    }
    
    list(NFL_Seg = NFL_Seg, Overlays = Overlays)
  })
  
  # Allows mapping of structural data onto CMP fundus image if aligned
  observe({
    if (input$isAlignedONH == "Aligned!"){
      enable("matchMapsCMPONH")
    }else{
      #CMP_MapCheckONH$Check = FALSE
      updateProgressBar(
        session = session,
        id = "pb1",
        value = 0, total = 100,
        title = "   "
      )
      disable("matchMapsCMPONH")
      disable("useOCTONH")
      updateCheckboxInput(session, "useOCTONH", value = FALSE)
      updateSelectInput(session, "segmListCMPONH", choices = c("None"))
    }
  })
  
  # Shows fundus image + structural overlays if available
  output$fundusImageONH <- renderPlot({
    CMP_Fundus <- CMP_Fundus()
    Structural_MapsCMP <- try({Structural_MapsCMPONH()}, silent = TRUE)
    
    if(input$segmListCMPONH != "None" & is.null(attr(Structural_MapsCMP, "condition"))){
      FoveaXY_CMP <- FoveaXY_CMP()
      par(mar = c(0, 0, 0, 0))
      if (input$segmListCMPONH == "SLO"){
        plot(Structural_MapsCMP$Overlays[[input$segmListCMPONH]], axes = FALSE)
      }else{
        plot(CMP_Fundus$Image, axes = FALSE)
        NFL_Seg <- Structural_MapsCMP$NFL_Seg
        points(NFL_Seg$x, NFL_Seg$y, pch = 20, cex = NFL_Seg$NFL/100,
               col = NFL_Seg$Col)
      }
    }else{
      par(mar = c(0, 0, 0, 0))
      plot(CMP_Fundus$Image, axes = FALSE)
    }
  })
  
  # Renders color map for peripapillary thickness
  output$fundusImagecolorbar <- renderPlot({
    Segm <- input$segmListCMPONH
    if(Segm != "None" & Segm != "SLO"){
      LimitsTh <- switch(Segm, 
                         "WRT" = c(200, 350),
                         "NFL" = c(0, 120),
                         "GCL" = c(0, 65),
                         "IPL" = c(0, 65),
                         "OUTER" = c(120, 200))
      par(mar=c(4,2,0,2))
      color.barh(colorRampPalette(Parula_Map$Hex)(100), LimitsTh[1], LimitsTh[2], 
                 title = expression(paste("Thickness ","(",mu,"m)")))
    }else{
      NULL
    }
  })
  
  ###############################################################
  ############### Test tab ##########################
  #Shows fundus image + Overlays 
  output$fundusImageTest <- renderPlot({
    CMP_Fundus <- CMP_Fundus()
    Structural_MapsCMP <- try({Structural_MapsCMP()}, silent = TRUE)
    FoveaXY_CMP <- FoveaXY_CMP()
    
    if(input$segmListCMPTest != "None" & 
       is.null(attr(Structural_MapsCMP, "condition"))){
      par(mar = c(0, 0, 0, 0))
      Image <- Structural_MapsCMP$Overlays[[input$segmListCMPTest]]
      plot(Image, axes = FALSE)
      points(FoveaXY_CMP[1], FoveaXY_CMP[2], col = 2, pch = 10, cex = 3)
      
      LimitsTh <- switch(input$segmListCMPTest, 
                         "WRT" = c(200, 350),
                         "NFL" = c(0, 120),
                         "GCL" = c(0, 65),
                         "IPL" = c(0, 65),
                         "OUTER" = c(120, 200))
      
    }else{
      par(mar = c(0, 0, 0, 0))
      plot(CMP_Fundus$Image, axes = FALSE)
      points(FoveaXY_CMP[1], FoveaXY_CMP[2], col = 2, pch = 10, cex = 3)
    }
    
    # Projects grid on image
    Grid_Coords_Pix <- Grid_Coords_Pix()
    if (dim(Grid_Coords_Pix)[1] > 0){
      if (input$plotDrasdo){
        LocX <- Grid_Coords_Pix$DrasdoX
        LocY <- Grid_Coords_Pix$DrasdoY
      }else{
        LocX <- Grid_Coords_Pix$X
        LocY <- Grid_Coords_Pix$Y
      }
      
      IsRGC <- Grid_Coords_Pix$DrasdoD
      
      if (!is.null(Test_Results$states)){
        finals <<- unlist(lapply(Test_Results$states, ZEST.final)) # get final estimates of threshold
        points(LocX, LocY, pch = 16, cex = (round(finals))/10,
               col = alpha(rgb(0,0,0), 0.4))
      }
      
      if (input$isSimulation){
        Sim_Th <- Sim_Th()
      }else{
        points(LocX, LocY, col = 1, pch = 16, cex = 1.5)
      }
      
      # Manually added locations defined under RGC displacement a in blue
      points(LocX[!IsRGC], LocY[!IsRGC], col = "red", pch = 16)
      points(LocX[!IsRGC], LocY[!IsRGC], col = "black", pch = 1)
      points(LocX[IsRGC], LocY[IsRGC], col = "blue", pch = 16)
      points(LocX[IsRGC], LocY[IsRGC], col = "black", pch = 1)
    }
    ########################################
    
    #Projects ONH
    ONH_Coords <- ONH_Coords()
    ONH_X <- ONH_Coords[1] #(32*onh[1] + 960)*Ratio
    ONH_Y <- ONH_Coords[2] #(32*onh[2] + 960)*Ratio
    points(ONH_X, ONH_Y, col = 1, pch = 1, cex = 2.5)
    points(ONH_X, ONH_Y, col = 2, pch = 1, cex = 2)
    
  })
  
  
  ######################   ######################   ######################   ######################   ###################### 
  ###################### Setting for the structural strategy (only for 10-2)
  # Preview macular structural data
  output$previewOCTStrat <- renderPlot({
    SLO <- SLO()
    if (!is.null(SLO)){
      Structural_Maps <- try({Structural_Maps <- Structural_Maps()}, silent = TRUE)
      if(input$segmListStrat != "None" & is.null(attr(Structural_Maps, "condition"))){
        par(mar = c(0, 0, 0, 0))
        Image <- Structural_Maps$Overlays[[input$segmListStrat]]
        plot(Image, axes = FALSE)
        
      }else{
        par(mar = c(0, 0, 0, 0))
        plot(SLO$Image, axes = FALSE)
      }
      
      SLO_TestCoords <- SLO_TestCoords()
      
      points(SLO_TestCoords$SLO_Test$X, SLO_TestCoords$SLO_Test$Y, col = 1, pch = 16)
      points(SLO_TestCoords$SLO_Test$DrasdoX, SLO_TestCoords$SLO_Test$DrasdoY, col = 3, pch = 16)
      points(SLO_TestCoords$SLO_ONH$X, SLO_TestCoords$SLO_ONH$Y, col = 2, pch = 16)
      
      #Checks if the fovea has been detected
      FoveaXY_SLO <- try({FoveaXY_SLO <- FoveaXY_SLO()}, silent = TRUE)
      if(is.null(attr(FoveaXY_SLO, "condition"))){
        points(FoveaXY_SLO[1], FoveaXY_SLO[2], col = 2, pch = 10, cex = 3)
      }
      
      #Checks if drasdo masks have been calculated
      DispMaskEdges <- SLO_StructData$DispMasksEdges
      if(!is.null(DispMaskEdges)){
        points(SLO_StructData$DrasdoX, SLO_StructData$DrasdoY, col = 2, pch = 16)
        for (i in unique(DispMaskEdges$Ind)){lines(DispMaskEdges$X[DispMaskEdges$Ind == i], 
                                                   DispMaskEdges$Y[DispMaskEdges$Ind == i])}
      }
    }
  })
  
  # Colorbar for macular map
  output$previewOCTcolorbarStrat <- renderPlot({
    Segm <- input$segmListStrat
    if(Segm != "None"){
      LimitsTh <- switch(Segm, 
                         "WRT" = c(200, 350),
                         "NFL" = c(0, 120),
                         "GCL" = c(0, 65),
                         "IPL" = c(0, 65),
                         "OUTER" = c(120, 200))
      par(mar = c(1, 3, 3, 2))
      color.barv(colorRampPalette(Parula_Map$Hex)(100), LimitsTh[1], LimitsTh[2], 
                 title = expression(atop("Thickness ", paste("(",mu,"m)", sep =""))))
    }else{
      NULL
    }
  })
  
  # Test coordinates translated onto the SLO image - used for structural calculations
  SLO_TestCoords <- reactive({
    #Pixel coordinates from CMP
    Grid_Coords <- Grid_Coords_Pix()
    
    #Matching
    Match <- SLO_Aligned()
    
    
    SLO_Test <- cbind(CMP_to_OCT(Grid_Coords$X, Grid_Coords$Y, Match), 
                      CMP_to_OCT(Grid_Coords$DrasdoX, Grid_Coords$DrasdoY, Match))
    
    colnames(SLO_Test) <- c("X", "Y", "DrasdoX", "DrasdoY")
    
    ONH_Coords <- ONH_Coords()
    SLO_ONH <- CMP_to_OCT(ONH_Coords[1], ONH_Coords[2], Match)
    
    list(SLO_Test = SLO_Test, SLO_ONH = SLO_ONH)
    
  })
  
  #Structural calculations
  SLO_StructData <- reactiveValues(DispMasksEdges = NULL,
                                   Av_Th = NULL,
                                   X = NULL,
                                   Y = NULL,
                                   DrasdoX = NULL,
                                   DrasdoY = NULL)
  
  # Calculates structural data - average thickness of Drasdo displaced locations
  observeEvent(input$calcStruct,{
    
    disable('calcStruct')
    enable('STRUCTURE')
    
    SLOIm <- SLO()
    
    SLO_TestCoords <- SLO_TestCoords()
    SLO_Test <- SLO_TestCoords$SLO_Test
    SLO_ONH <- SLO_TestCoords$SLO_ONH
    
    FoveaXY_SLO <- FoveaXY_SLO()
    
    Structural_Maps <- Structural_Maps()
    ThMap <- Structural_Maps$SEG_List[[input$segmListStrat]]
    
    Eye <- input$selectEye
    
    
    ############################# SF Calculations #############################
    SLO_Sz <- dim(SLOIm$Image)[1]
    DegToRad <- SLOIm$FoV/SLO_Sz
    Fovea_Cent <- data.frame(X = (SLO_Test$X - FoveaXY_SLO[1])*DegToRad, 
                             Y = (SLO_Test$Y - FoveaXY_SLO[2])*DegToRad)
    Fovea_Cent$Ecc <- sqrt(Fovea_Cent$X^2 + Fovea_Cent$Y^2)
    ONH_FC <- (as.numeric(SLO_ONH) - FoveaXY_SLO)*DegToRad
    
    MaskEdges <- create_StimMaskEdge(Fovea_Cent$X, Fovea_Cent$Y)
    
    ###########################################################################
    
    #Drasdo displacement
    DrasdoDisp <- drasdoFast(Fovea_Cent$X, Fovea_Cent$Y, ONH = ONH_FC[1:2], 
                             SupPositive = FALSE, Eye = Eye)
    
    XDisp <- DrasdoDisp[[1]]$xdispDeg/DegToRad + FoveaXY_SLO[1]
    YDisp <- DrasdoDisp[[1]]$ydispDeg/DegToRad + FoveaXY_SLO[2]
    
    #Drasdo masks
    DispMasksEdges <- drasdoFast(MaskEdges$xF, MaskEdges$yF, ONH = ONH_FC[1:2],
                                 SupPositive = FALSE, Eye = Eye)
    DispMasksEdges <- DistMaskAsEdges(DispMasksEdges, MaskEdges)
    
    DispMasks <- create_StimMasks(DispMasksEdges[[1]], FoveaXY = FoveaXY_SLO, 
                                  ImSize = SLO_Sz, FoVX = SLOIm$FoV, FoVY = SLOIm$FoV)
    
    
    #Fast average calculation Combined masks
    Av_Th <- dim(DispMasks)[3]
    #CombMask <- matrix(0, dim(DispMasks)[1], dim(DispMasks)[2])
    for (j in 1:dim(DispMasks)[3]){
      Mask <- as.double(t(DispMasks[,,j]))
      Av_Th[j] <- mean(ThMap[Mask == 1], na.rm = TRUE)
      #CombMask <- CombMask + Mask*j
    }
    
    DispMasksEdges <- data.frame(X = DispMasksEdges[[1]]$xF/DegToRad + FoveaXY_SLO[1], 
                                 Y = DispMasksEdges[[1]]$yF/DegToRad + FoveaXY_SLO[2],
                                 Ind = DispMasksEdges[[1]]$Ind) 
    
    GridCoords <- GridCoords()
    
    SLO_StructData$DispMasksEdges = DispMasksEdges
    SLO_StructData$Av_Th = Av_Th
    SLO_StructData$X = GridCoords$X
    SLO_StructData$Y = GridCoords$Y
    SLO_StructData$DrasdoX = XDisp
    SLO_StructData$DrasdoY = YDisp
    
    SLO_StructData
  })
  
  ######################   ######################   ######################   ######################   ###################### 
  ######################   ######################   ######################   ######################   ###################### 
  
  # Colorbar for the thickness maps reported onto the fundus image
  output$fundusImagecolorbarTest <- renderPlot({
    Segm <- input$segmListCMPTest
    if(Segm != "None" & Segm != "SLO"){
      LimitsTh <- switch(Segm, 
                         "WRT" = c(200, 350),
                         "NFL" = c(0, 120),
                         "GCL" = c(0, 65),
                         "IPL" = c(0, 65),
                         "OUTER" = c(120, 200))
      par(mar=c(4,2,0,2))
      color.barh(colorRampPalette(Parula_Map$Hex)(100), LimitsTh[1], LimitsTh[2], 
                 title = expression(paste("Thickness ","(",mu,"m)")))
    }else{
      NULL
    }
  })
  
  
  # Original testing grid coordinates - in degrees and pixels
  # Changes according to whether the fovea or the PRL is taken as the grid centre
  GridCoords <-  reactive({
    
    Checks <- FALSE
    
    #Selects grids if new baseline
    if (IsBaseline$IsBaseline){
      
      if (input$standGrid == "None"){
        X_Perim <- double()
        Y_Perim <- double()
      }
      
      if (input$standGrid == "10-2"){
        X_Perim <- Coord_102$X
        Y_Perim <- Coord_102$Y
      }
      
      if (input$standGrid == "24-2"){
        X_Perim <- Coord_242$X
        Y_Perim <- Coord_242$Y
      }
      
      if (input$standGrid == "30-2"){
        X_Perim <- Coord_302$X
        Y_Perim <- Coord_302$Y
      }
      
      if (input$standGrid == "Upload your perimetric grid"){
        try({
          UserGrid <- input$user_grid
          UserCoords <- read.table(text = gsub(";",",", readLines(UserGrid$datapath)),
                                   header = TRUE, sep = ",")
          X_Perim <- UserCoords[,1]
          Y_Perim <- UserCoords[,2]
        }, silent = TRUE)
      }
      
      DF <- try({data.frame(X = X_Perim, Y = Y_Perim, DrasdoD = numeric(length(X_Perim)) > 0)})
      
      #Left eye flip horizontal
      if (input$selectEye == "Left"){
        DF$X <- -DF$X
      }
      
      #Uses ONH in degrees (retina convetion, so no invertion top-bottom)
      #Rotation happens in degrees, so no shift for fovea, which is done later
      if (input$alignON){
        CenterX <- 0
        CenterY <- 0
        
        if (input$gridCenter == "PRL"){
          CenterX <- prl[1]
          CenterY <- prl[2]
        }
        
        #Projects ONH
        ONH_Coords <- ONH_Coords()
        ONH_X <- (ONH_Coords[1]/Ratio - 960)/32 #(32*onh[1] + 960)*Ratio
        ONH_Y <- (ONH_Coords[2]/Ratio - 960)/32 #(32*onh[2] + 960)*Ratio
        
        FovDisc_Angle <- -foveaDisc_Angle(CenterX, CenterY, ONH_X, ONH_Y)
        R <- rbind(c(cos(FovDisc_Angle), -sin(FovDisc_Angle)), c(sin(FovDisc_Angle), cos(FovDisc_Angle)))
        XY <- R%*%rbind(DF$X, DF$Y)
        DF$X <- XY[1,]
        DF$Y <- XY[2,]
        
      }
      
      if (nrow(DF) > 0){DF$StructData <- NA}
      
      ###############
    }else{
      
      FoveaXY_CMP <- FoveaXY_CMP()
      
      #Loads previous if follow-up
      TempEnv <- new.env()
      load(input$baseline_Test$datapath, envir = TempEnv)
      
      CMP_Fundus <- CMP_Fundus()
      
      #Needs to convert from pixels back to degrees for consistency (after transforming from baseline image)
      CoordsPix_Base <- TempEnv$Grid_Coords_Pix[,1:2]
      CoordsPix <- applyTransform(reverse(CMP_Fundus$MatchFU), as.matrix(CoordsPix_Base))
      
      DF <- data.frame(X = CoordsPix[,1], Y = CoordsPix[,2])
      
      if(input$gridCenter == "Fovea"){
        DF$X <- (DF$X - FoveaXY_CMP[1])/32/Ratio
        DF$Y <- -(DF$Y - FoveaXY_CMP[2])/32/Ratio 
      }else{
        PRLOffset_X <- (prl[1]*32 + 960)*Ratio
        PRLOffset_Y <- (prl[2]*32 + 960)*Ratio
        DF$X <- (DF$X - PRLOffset_X)/32/Ratio
        DF$Y <- -(DF$Y - PRLOffset_Y)/32/Ratio 
      }
      
      DF$StructData <- TempEnv$StructData
      DF$DrasdoD <- TempEnv$Grid_Coords_Pix$DrasdoD
      
      if(is.null(DF$StructData)){DF$StructData <- NA}
      
      if (all(!is.na(DF$StructData))){
        enable('STRUCTURE') 
        Checks <- TRUE
      }
      
    }
    
    ####################################################################################
    ##############    #############      #############      #############      #############
    # Adds manual coordinates if present (no need to flip X axis, no rotation)
    # Fovea and ONH on CMP in scaled pixels
    ONH_Coords <- ONH_Coords()
    FoveaXY_CMP <- FoveaXY_CMP()
    AddLoc <- ManualCoords$DF
    
    # Recenters manual coordinates so they do not change based on slected centre
    if(input$gridCenter == "Fovea"){
      AddLoc$X <- (AddLoc$X - FoveaXY_CMP[1])/32/Ratio
      AddLoc$Y <- -(AddLoc$Y - FoveaXY_CMP[2])/32/Ratio 
    }else{
      PRLOffset_X <- (prl[1]*32 + 960)*Ratio
      PRLOffset_Y <- (prl[2]*32 + 960)*Ratio
      AddLoc$X <- (AddLoc$X - PRLOffset_X)/32/Ratio
      AddLoc$Y <- -(AddLoc$Y - PRLOffset_Y)/32/Ratio 
    }
    
    
    # Transforms Drasdo back into retinal coordinates
    ONH_CoordsC <- (ONH_Coords - FoveaXY_CMP)/Ratio/32
    DrasdoDisp <- drasdoFast(AddLoc$X[AddLoc$DrasdoD], AddLoc$Y[AddLoc$DrasdoD], ONH = ONH_CoordsC,
                             SupPositive = FALSE, Eye = input$selectEye, InvertD = TRUE)[[1]]
    
    # This adds the additional locations (in original retinal coordinates)
    AddLocOr <- AddLoc
    AddLocOr$DrasdoD <- FALSE
    AddLoc$X[AddLoc$DrasdoD] <- DrasdoDisp$xdispDeg
    AddLoc$Y[AddLoc$DrasdoD] <- DrasdoDisp$ydispDeg
    
    AddLoc$StructData <- NA
    AddLocOr$StructData <- NA
    
    DF <- rbind(DF, AddLoc)
    ##################################################################################################
    
    
    
    if(!is.null(SLO_StructData$Av_Th)){
      if (length(SLO_StructData$Av_Th) == length(DF$X)){
        if (all(SLO_StructData$X == DF$X) & all(SLO_StructData$Y == DF$Y)){
          DF$StructData <- SLO_StructData$Av_Th
          disable('calcStruct')
          enable('STRUCTURE')
          Checks <- TRUE
        }
      }
    }
    
    if (!Checks){
      SLO_StructData$Av_Th <- NULL
      SLO_StructData$DispMasksEdges <- NULL
      SLO_StructData$X <- NULL
      SLO_StructData$Y <- NULL
      enable('calcStruct')
      disable('STRUCTURE')
    }
    
    #Calculates interpolated normative value (based on eccentricity, in degrees)
    DF$Ecc <- sqrt(DF$X^2 + DF$Y^2)
    DF$NormVal <- CalcNorm(input$age, DF$Ecc) 
    
    #Enables startTest only if at least one location is tested
    if (dim(DF)[1] > 0){
      enable("startTest")
    }else{
      disable("startTest")
    }
    
    DF
  })
  
  # Compass coordinates in pixels (Centre reference 960,960)
  Grid_Coords_Pix <- reactive({
    #Grid coordinates in degrees
    GridCoords <- GridCoords()
    enable('calcStruct')
    
    if (dim(GridCoords)[2] > 0){
      
      #Fovea on CMP in scaled pixels
      FoveaXY_CMP <- FoveaXY_CMP()
      ONH_Coords <- ONH_Coords()
      
      #Uses PRL if fovea not provided
      if (is.na(FoveaXY_CMP[1])){
        FoveaXY_CMP <- (prl*32 + 960)*Ratio
      }
      
      #Scales by the chosen ratio (only for display, projections are in degrees)
      #Centers on the fovea
      Fovea_OffsetX <- (FoveaXY_CMP[1] - 960*Ratio)
      Fovea_OffsetY <- (FoveaXY_CMP[2] - 960*Ratio)
      
      if(input$gridCenter == "Fovea"){
        GridCoords$X <- (32*GridCoords$X + 960)*Ratio + Fovea_OffsetX
        GridCoords$Y <- (-32*GridCoords$Y + 960)*Ratio + Fovea_OffsetY 
      }else{
        GridCoords$X <- GridCoords$X + prl[1]
        GridCoords$Y <- -GridCoords$Y + prl[2]
        XPrlDeg <- GridCoords$X
        YPrlDeg <- GridCoords$Y
        
        GridCoords$X <- (32*GridCoords$X + 960)*Ratio
        GridCoords$Y <- (32*GridCoords$Y + 960)*Ratio
      }
      
      
      #Drasdo displacement (from anatomical fovea). 
      #They are assumed in retinal coordinates, right eye, in visual degrees
      
      GridCoords$XFoveaDeg <- (GridCoords$X - FoveaXY_CMP[1])/Ratio/32
      GridCoords$YFoveaDeg <- (GridCoords$Y - FoveaXY_CMP[2])/Ratio/32
      ONH_CoordsC <- (ONH_Coords - FoveaXY_CMP)/Ratio/32
      
      DrasdoDisp <- drasdoFast(GridCoords$XFoveaDeg, GridCoords$YFoveaDeg, ONH = ONH_CoordsC,
                               SupPositive = FALSE, Eye = input$selectEye)[[1]]
      
      GridCoords$DrasdoX <- DrasdoDisp$xdispDeg*32*Ratio + FoveaXY_CMP[1]
      GridCoords$DrasdoY <- DrasdoDisp$ydispDeg*32*Ratio + FoveaXY_CMP[2]
      GridCoords
    }
  })
  
  
  # This converts back to VF coordinates for OPI 
  Test_Locations <- reactive({
    #Pixel to VF conversion. The coordinates refer to the centre of the image (960, 960) -> (0,0)
    #No need to invert the Y axis, positive values are superior, negative values inferior.
    GridCoords <- Grid_Coords_Pix()
    
    GridCoords$X <- (GridCoords$X/Ratio - 960)/32
    GridCoords$Y <- (GridCoords$Y/Ratio - 960)/32
    
    if(input$isSimulation){
      Sim_Th <- Sim_Th()
      Sim_Th[is.na(Sim_Th)] <- 0
      GridCoords$Threshold <- Sim_Th
    }else{
      GridCoords$Threshold <- 0
    }
    
    #Reorders for simulations/projections
    GridCoords[,c("X","Y","Threshold","XFoveaDeg","YFoveaDeg","DrasdoX","DrasdoY","Ecc","NormVal","StructData")]
  })
  
  
  # Shows perimetric coordinates and corresponding structural data if available 
  output$grid_perim <- renderDataTable({
    GridCoords <- GridCoords() 
    if (nrow(GridCoords) > 0){
      GridCoords <- GridCoords[,c("X","Y","StructData")]
      GridCoords$X <- round(GridCoords$X, 2)
      GridCoords$Y <- round(GridCoords$Y, 2)
      GridCoords$StructData <- round(GridCoords$StructData, 2)
      GridCoords
    }else{
      GridCoords[,c("X","Y")]
    }
  }, options = list(pageLength = 200, dom = 't', scrollY = "100px"))
  
  
  #Manual location input
  ManualCoords <- reactiveValues(
    DF = NULL
  )
  
  # Clears all manually added locations
  observeEvent(input$clearManual, {
    ManualCoords$DF = NULL
  })
  
  # Adds manual locations via click
  observeEvent(input$plot_click_add, {
    #In scaled pixel coords
    AddLoc <-  data.frame(X = input$plot_click_add$x, Y = input$plot_click_add$y, DrasdoD = input$plotDrasdo)
    ManualCoords$DF = rbind(ManualCoords$DF, AddLoc)
    
  })
  
  # Removes individual locations via double click
  observeEvent(input$plot_dblclick_add, {
    
    #In scaled pixel coords
    RemoveLoc <-  data.frame(X = input$plot_dblclick_add$x, 
                             Y = input$plot_dblclick_add$y)
    #Tolearance in pixels
    Tolerance <- 30*Ratio
    
    if (!is.null(ManualCoords$DF)){
      MaualLocs <- ManualCoords$DF
      Sel_D <- sqrt((MaualLocs$X - RemoveLoc$X)^2 + (MaualLocs$Y - RemoveLoc$Y)^2)
      MaualLocs <- MaualLocs[Sel_D > Tolerance, ]
      ManualCoords$DF = MaualLocs
    }
    
  })
  
  # Shows manual coordinates
  output$manual_add <- renderDataTable({
    if (!is.null(ManualCoords$DF)){
      DF <- ManualCoords$DF
      DF$DrasdoD <- NULL
      round(DF, 2)
    }else{
      NULL
    }
  }, options = list(pageLength = 200, dom = 't', scrollY = "100px"))
  
  
  
  
  ########### CODE FOR PERIMETRIC TEST ###########
  Sim_Th <- reactive({
    
    Grid_Coords_Pix <- Grid_Coords_Pix()
    
    if (input$isSimulation & !sum(is.na(Grid_Coords_Pix$X))>0){
      
      if (IsBaseline$IsBaseline){
        
        #Left eye flip horizontal, data stored in RE convention
        simData <- read.csv("data_Temp_Exam/Test_Vals.csv", sep = ",")
        if (input$selectEye == "Left"){
          simData$X <- -simData$X
        }
        
        
        simData$X <- (32*(simData$X + prl[1]) + 960)*Ratio
        simData$Y <- (32*(simData$Y + prl[2]) + 960)*Ratio 
        
        suppressWarnings(Sim_Th <- akima::interp(simData$X, simData$Y, simData$Th))
        
        #Ensures at least one non NA
        HV <- Sim_Th$z[round(dim(Sim_Th$z)[1]/2),]
        VV <- Sim_Th$z[,round(dim(Sim_Th$z)[2]/2)]
        HV[is.na(HV)] <- HV[!is.na(HV)][c(1, length(HV[!is.na(HV)]))]
        VV[is.na(VV)] <- VV[!is.na(VV)][c(1, length(VV[!is.na(VV)]))]
        Sim_Th$z[round(dim(Sim_Th$z)[1]/2),] <- HV
        Sim_Th$z[,round(dim(Sim_Th$z)[2]/2)] <- VV
        
        Sim_Th$z <- na.spline(Sim_Th$z)
        
        #Z needs to be transposed for interpolation
        Sim_Th <- interp2(Sim_Th$x, Sim_Th$y, t(Sim_Th$z),  Grid_Coords_Pix$X,  Grid_Coords_Pix$Y)
      }else{
        TempEnv <- new.env()
        load(input$baseline_Test$datapath, envir = TempEnv)
        Sim_Th <- TempEnv$finals
      }
      
    }else{
      Sim_Th <- NULL
    }
    Sim_Th
  })
  
  ##### Testing bit #########################
  # Test results
  Test_Results <- reactiveValues(
    finals = NULL,
    states = NULL,
    presentations = NULL
  )
  
  # Control variables for user input
  Test_Controls <- reactiveValues(
    Run = FALSE,
    Pause = FALSE,
    Stop = FALSE,
    OngoingTest = FALSE,
    Aborted = FALSE,
    FinalString = "Testing"
  )
  
  # Test start
  observeEvent(input$startTest, {
    
    #Disables element in test settings
    disable("standGrid")
    disable("user_grid")
    disable("selectEye")
    disable("gridCenter")
    disable("alignON")
    disable("clearManual")
    
    
    disable("startTest")
    enable("pauseTest")
    enable("stopTest")
    
    #Date and time of test start (only if NULL)
    if (is.null(startTime)){startTime <<- Sys.time()}
    
    # Pause/restart if Paused, else starts a new test
    if(Test_Controls$Pause){
      Test_Controls$Run = TRUE
      Test_Controls$Pause = FALSE
      pauseDuration <<- pauseDuration + (Sys.time() - pauseStart)
    }else{
      #This bit sets up a new test
      
      Test_Results$finals = NULL
      Test_Results$presentations = NULL
      Test_Results$states = NULL
      
      Test_Controls$Run = TRUE
      Test_Controls$Pause = FALSE
      
      SIMULATION <- input$isSimulation
      
      #### Test parameters #######
      locations <<- Test_Locations()
      
      ######
      # These are set up only for simulations
      if (SIMULATION){locations$Threshold[locations$Threshold < 0] <<- 0}
      FPR <- .05
      FNR <- .05
      ######
      
      STRUCTURAL <<- input$STRUCTURE
      SAge <<- input$age
      
      #Calculates neighborhood - function defined in global.R
      NHs <<- NHood(locations$XFoveaDeg, locations$YFoveaDeg)
      
      if (!STRUCTURAL){
        #Spatial enhancement
        SPATIAL <<- FALSE
        #Correlation strength
        Corr <<- 0.2
        
        source("OPI_Functions/setupZEST.R", local = TRUE)
      }else{
        #Spatial enhancement
        SPATIAL <<- TRUE
        #Correlation strength
        Corr <<- 0.4
        
        source("OPI_Functions/setupZEST_Str.R", local = TRUE)
      }
      
      #Code for strategy. Iterates until all locations have reached stopping criterion
      #Declares useful variables as global to be accessed by the observer (next bunch of code)
      st <<- (unlist(lapply(states, ZEST.stop)))
      
      updateProgressBar(
        session = session,
        id = "pb1",
        value = 0, total = 100,
        title = paste("Testing", " (Closed locations: ", sum(st), "/",
                      length(st),", FP: ",round(EstFPR*100),"%)", sep = "")
      )
      
      #This is just to export the sequence of testing
      iout <<- NULL
      
      #This is a replacement for a while loop. Shiny cannot interact to stop while loops.
      #This hack is necessary to allow pauses and to stop the test interactively
      #The code ZEST_Step.R interacts with global variables, declared above.
      observe({
        
        #This is for one step of ZEST
        source("OPI_Functions/ZEST_Step.R", local = TRUE)
        
        #Updates progress bar
        updateProgressBar(
          session = session,
          id = "pb1", 
          value = round(mean(ProgPercSD)*100), total = 100,
          title = paste("Testing", " (Closed locations: ", sum(st), "/",
                        length(st),", FP: ",round(EstFPR*100),"%)", sep = "")
        )
        
        # Checks if all locations have reached the termination criterion
        if(!all(st) & Test_Controls$Run){
          invalidateLater(0, session)
        }else{
          
          #Updates progress bar
          updateProgressBar(
            session = session,
            id = "pb1", 
            value = round(mean(ProgPercSD)*100), total = 100,
            title = paste(Test_Controls$FinalString, " (Closed locations: ", sum(st), "/",
                          length(st),", FP: ",round(EstFPR*100),"%)", sep = "")
          )
          
          if(!Test_Controls$Pause){
            
            Test_Controls$Pause = FALSE
            Test_Controls$OngoingTest = FALSE
            Test_Controls$Run = FALSE
            
            if(!Test_Controls$Aborted){Test_Controls$FinalString = "Finished!"}
            
            #Calculates duration of the test
            TestDuration <<- Sys.time() - startTime
            
            updateProgressBar(
              session = session,
              id = "pb1",  
              value = round(mean(ProgPercSD)*100), total = 100,
              title = paste(Test_Controls$FinalString, " (Closed locations: ", sum(st), "/",
                            length(st),", FP: ",round(EstFPR*100),"%)", sep = "")
            )
            
            disable("pauseTest")
            disable("stopTest")
            enable("startTest")
            
            Test_Results$states = states
          }
        }
      })
    }
  })
  
  # Pauses the test
  observeEvent(input$pauseTest,{
    disable("pauseTest")
    enable("startTest")
    Test_Controls$Pause = TRUE
    Test_Controls$OngoingTest = TRUE
    Test_Controls$Run = FALSE
    Test_Controls$Aborted = FALSE
    Test_Controls$FinalString = "Paused"
    pauseStart <<- Sys.time()
  })
  
  # Stops the test (cannot be resumed after, saves current progress)
  observeEvent(input$stopTest,{
    confirmSweetAlert(
      session = session,
      inputId = "myconfirmation",
      title = "Do you want to abort the test?"
    )
  })
  
  # Confirms that test is to be stopped
  observeEvent(input$myconfirmation, {
    if (isTRUE(input$myconfirmation)) {
      disable("pauseTest")
      disable("stopTest")
      enable("startTest")
      Test_Controls$Pause = FALSE
      Test_Controls$OngoingTest = FALSE
      Test_Controls$Run = FALSE
      Test_Controls$Aborted = TRUE
      Test_Controls$FinalString = "Aborted!"
    }
  })
  
  # Shows final results as a table and saves final results
  output$showResults <- renderDataTable({
    GridCoords <- GridCoords()
    Grid_Coords_Pix <- Grid_Coords_Pix()
    
    if (!is.null(Test_Results$states)){
      finals <<- unlist(lapply(Test_Results$states, ZEST.final)) # get final estimates of threshold
      presentations <<- unlist(lapply(Test_Results$states, getPres)) # get final number of presentations
      
      ########################################################################################
      #This terminates the session. A new testing session needs to be started for more testing.
      #Uncomment to start multiple tests (only for development)
      Fixation <<- opiClose()
      ########################################################################################
      
      #Prepares data for saving
      CMP_FundusOr <- CMP_FundusOr()
      CMP_Fundus <- CMP_Fundus()
      Ratio <- Ratio
      FoveaXY_CMP <- FoveaXY_CMP()
      FoveaXY_SLO <- FoveaXY_SLO()
      if (input$isAligned == "Aligned!"){
        SLO_Aligned <- SLO_Aligned()
      }else{
        SLO_Aligned <- NULL
      }
      
      if (input$isAlignedONH == "Aligned!"){
        SLO_AlignedONH <- SLO_AlignedONH()
      }else{
        SLO_AlignedONH <- NULL
      }
      
      FoveaXY_CMP <- FoveaXY_CMP/Ratio
      
      PatName <- input$Pat_ID
      Eye <- input$selectEye
      
      DateTime <- Sys.time()
      
      User <- input$userName
      Age <- input$age
      Study <- input$study
      GridCenter <- input$gridCenter
      Grid <- input$standGrid
      ONH_Coords <- ONH_Coords()
      
      BaselineTest <- IsBaseline$IsBaseline
      
      StructData <- SLO_StructData$Av_Th
      
      
      #Saves the test
      FileName <- paste("data_Temp_Exam/","Test", "_",PatName,"_",Eye,"_",as.integer(DateTime), sep = "")
      save(file = paste(FileName,".RData", sep = ""), 
           PatName, Eye, DateTime, User, Study, Age, startTime, TestDuration, pauseDuration,
           finals, presentations, states, FPsC, TotLisTime, EstFPR, BaselineTest, Grid, STRUCTURAL,
           Fixation, prl, onh, ONH_Coords, FoveaXY_CMP, FoveaXY_SLO, GridCoords, Grid_Coords_Pix, 
           Ratio, SLO_Aligned, SLO_AlignedONH, GridCenter, CMP_FundusOr, CMP_Fundus, StructData)
      
      #Saves matching info if baseline for follow-up tests
      if (IsBaseline$IsBaseline){
        save(file = paste(FileName,".BSInfo", sep = ""), 
             PatName, Eye, DateTime, User, Study, Age, prl, onh, FoveaXY_CMP, FoveaXY_SLO, ONH_Coords,
             GridCoords, Grid_Coords_Pix, Ratio, CMP_Fundus, GridCenter, SLO_Aligned, SLO_AlignedONH,
             StructData, Grid = Grid, prevTh = finals)
      }
      
      GridCoords <- GridCoords[,c("X","Y","StructData")]
      GridCoords$X <- round(GridCoords$X, 2)
      GridCoords$Y <- round(GridCoords$Y, 2)
      GridCoords$StructData <- round(GridCoords$StructData, 2)
      
      cbind(round(GridCoords,2), data.frame(Thresholds = round(finals), 
                                            Presentations = presentations))
    }
  }, options = list(pageLength = 200, dom = 't', scrollY = "100px"))
  
}