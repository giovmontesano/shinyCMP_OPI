################## DASHBOARD ######################
ui <- dashboardPage(
  dashboardHeader(title = 'Compass OPI'),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Connect to Compass", tabName = "connectCMP", icon = icon("laptop-code")),
      menuItemOutput("recAlign"),
      menuItemOutput("recTest")
    )
  ),
  
  dashboardBody(
    
    useShinyjs(),
    
    tabItems(
      
      ########### Connect to CMP and load structure tab ####################
      tabItem(tabName = "connectCMP",
              fluidRow(
                box(title = NULL, width = 12,
                    box(title = "Test data", width = 12,
                        column(4,
                               textInput("Pat_ID", label = "Patient ID", value = NULL),
                               numericInput("age", label = "Age (years)", value = NULL)
                        ),
                        column(4,
                               textInput("study", label = "Study", value = NULL),
                               textInput("userName", label = "User", value = NULL)
                        ),
                        column(3,
                               fileInput("baseline_Test", label = "Load baseline test", accept = ".BSInfo"),
                               disabled(
                                 actionButton("resets_Baseline", label = "New baseline")
                               )
                        )
                    ),
                    
                    box(title = "Load Structural Data", width = 6, 
                        #Creates tables for upload of multiple structural data
                        tabsetPanel(type = "tabs",
                                    #Table for macular scan
                                    tabPanel("Macular scan",
                                             selectInput("octLIst", label = "Load Spectralis data", choices = paste("data_Temp_Exam/",filelist, sep = "")),
                                             fluidRow(
                                               column(3,
                                                      actionButton("buttonStruct", "Load Spectralis")
                                               ),
                                               column(3,
                                                      disabled(
                                                        actionButton("buttonSegm", "Segmentations")
                                                      )
                                               ),
                                               column(3,
                                                      selectInput("segmList", label = NULL, choices = "None")
                                               ),
                                               column(3,
                                                      disabled(
                                                        actionButton("buttonFovea", "Find fovea")
                                                      )
                                               )
                                             ),
                                             br(),
                                             fluidRow(
                                               column(3,
                                                      plotOutput("previewOCTcolorbar")
                                               ),
                                               column(9,
                                                      plotOutput("previewOCT", click = "plot_click_fovea") %>% withSpinner(color="#0dc5c1")
                                               )
                                             )
                                    ),
                                    
                                    #Table for the optic nerve head
                                    tabPanel("Optic Nerve Head",
                                             selectInput("octLIstONH", label = "Load Spectralis data", choices = paste("data_Temp_Exam/",filelist, sep = "")),
                                             fluidRow(
                                               column(3,
                                                      actionButton("buttonStructONH", "Load Spectralis")
                                               ),
                                               column(3,
                                                      disabled(
                                                        actionButton("buttonSegmONH", "Segmentations")
                                                      )
                                               ),
                                               column(3,
                                                      selectInput("segmListONH", label = NULL, choices = "None")
                                               )
                                             ),
                                             br(),
                                             fluidRow(
                                               column(3,
                                                      plotOutput("previewOCTcolorbarONH")
                                               ),
                                               column(9,
                                                      plotOutput("previewOCTONH") %>% withSpinner(color="#0dc5c1")
                                               )
                                             )
                                    )
                        )
                    ),
                    box(title = "Connect to Compass", width = 6,
                        fluidRow(
                          column(6,
                                 textInput("CMP_Address", label = "Compass Address", value = "192.168.50.1")
                          ),
                          column(2,
                                 br(),
                                 actionButton("buttonConnect", "Connect")
                          ),
                          column(3,
                                 br(),
                                 checkboxInput("isSimulation", "Simulation", value = TRUE)
                          )
                        ),
                        disabled(
                          textInput("OUTPUT_IMAGE_S", label = "Fundus image path", value = NULL)
                        ),
                        plotOutput("previewCMP") %>% withSpinner(color="#0dc5c1")
                    )
                )
              )
      ),
      
      ########### Fundus alignment tab ####################
      tabItem(tabName = "fundusAlign",
              #Three tabs, one for macular scans, one for ONH and one for FAF
              tabsetPanel(type = "tabs", id = "alignTabs",
                          #Table for macular scan
                          tabPanel("Macular scan", value = "Macular_scan",
                                   fluidRow(
                                     box(title = paste("Spectralis SLO + Compass fundus image (30°)", sep = ""), width = 6, style = "font-size: 18px",
                                         fluidRow(
                                           column(6,
                                                  actionButton("align", "Align!"),
                                                  disabled(
                                                    textInput("isAligned", label = NULL, value = "Not aligned!")
                                                  )
                                           ),
                                           column(6,
                                                  disabled(
                                                    checkboxInput("showAlign", "Show aligned", value = FALSE, width = NULL)
                                                  ),
                                                  disabled(
                                                    checkboxInput("baselineAlign", "Use baseline alignment", value = FALSE, width = NULL)
                                                  )
                                           )
                                         ),
                                         plotOutput("alingnedPlot", width = "100%", height = "600px")
                                     ),
                                     
                                     box(title = paste("Compass fundus image (60°)", sep = ""), width = 6, style = "font-size: 18px",
                                         bsPopover(id = "fundusImage", title = NULL,
                                                   content = "Clicking will move the center of Compass image",
                                                   placement = "auto", options = list(container = "body")),
                                         fluidRow(
                                           column(4,
                                                  checkboxInput("move", "Live move", value = FALSE, width = NULL)
                                           ),
                                           column(4,
                                                  disabled(
                                                    actionButton("matchMapsCMP", "Report OCT on Compass")
                                                  )
                                           ),
                                           column(3,
                                                  selectInput("segmListCMP", label = NULL, choices = "None")
                                           )
                                         ),
                                         
                                         
                                         plotOutput("fundusImage", width = "100%", height = "600px",
                                                    hover = "plot_hover",
                                                    click = "plot_click",
                                                    hoverOpts(id = "plot_hover", delay = 0, delayType = c("debounce"))),
                                         
                                         plotOutput("fundusImagecolorbar", width = "100%", height = "70px")
                                         
                                         
                                     )
                                   )
                          ),
                          
                          tabPanel("Optic nerve scan", value = "ONH_scan",
                                   fluidRow(
                                     box(title = paste("Spectralis SLO + Compass fundus image (30°)", sep = ""), width = 6, style = "font-size: 18px",
                                         fluidRow(
                                           column(6,
                                                  actionButton("alignONH", "Align!"),
                                                  disabled(
                                                    textInput("isAlignedONH", label = NULL, value = "Not aligned!")
                                                  )
                                           ),
                                           column(6,
                                                  disabled(
                                                    checkboxInput("showAlignONH", "Show aligned", value = FALSE, width = NULL)
                                                  ),
                                                  disabled(
                                                    checkboxInput("baselineAlignONH", "Use baseline alignment", value = FALSE, width = NULL)
                                                  )
                                           )
                                         ),
                                         plotOutput("alingnedPlotONH", width = "100%", height = "600px")
                                     ),
                                     box(title = paste("Compass fundus image (60°)", sep = ""), width = 6, style = "font-size: 18px",
                                         bsPopover(id = "fundusImageONH", title = NULL,
                                                   content = "Clicking will move the center of Compass image",
                                                   placement = "auto", options = list(container = "body")),
                                         fluidRow(
                                           column(4,
                                                  checkboxInput("moveONH", "Live move", value = FALSE, width = NULL)
                                           ),
                                           column(4,
                                                  disabled(
                                                    actionButton("matchMapsCMPONH", "Report OCT on Compass")
                                                  )
                                           ),
                                           column(3,
                                                  selectInput("segmListCMPONH", label = NULL, choices = "None")
                                           )
                                         ),
                                         
                                         
                                         plotOutput("fundusImageONH", width = "100%", height = "600px",
                                                    hover = "plot_hoverONH",
                                                    click = "plot_clickONH",
                                                    hoverOpts(id = "plot_hoverONH", delay = 0, delayType = c("debounce"))),
                                         
                                         plotOutput("fundusImagecolorbarONH", width = "100%", height = "70px")
                                         
                                         
                                     )
                                   )
                          )
              )
      ),
      
      
      
      ########### Test tab ####################
      tabItem(tabName = "testTab",
              fluidRow(
                box(title = paste("Test settings", sep = ""), width = 6, style = "font-size: 18px",
                    #Three tabs, one for macular scans, one for ONH and one for FAF
                    tabsetPanel(type = "tabs", id = "settingTabs",
                                #General settings
                                tabPanel("General settings", value = "set_Tabs",
                                         fluidRow(
                                           column(6,
                                                  radioButtons("standGrid", label = "Select a standard grid",
                                                               choices = c("None","10-2", "24-2", "30-2", "Upload your perimetric grid"),
                                                               selected = NULL),
                                                  
                                                  fileInput("user_grid", label = NULL)
                                           ),
                                           column(6,
                                                  radioButtons("selectEye", label = "Eye",
                                                               choices = c("Right","Left"),  
                                                               selected = NULL),
                                                  radioButtons("gridCenter", label = "Center of the grid",
                                                               choices = c("Fovea","PRL"),  
                                                               selected = NULL),
                                                  
                                                  checkboxInput("alignON", strong("Align with Fovea - disc angle"), FALSE),
                                                  disabled(
                                                    checkboxInput("STRUCTURE", strong("Use structural data"), FALSE)
                                                  )
                                           )
                                         ),
                                         fluidRow(
                                           box(title = "Grid locations", width = 6,
                                               DT::dataTableOutput('grid_perim')
                                           ),
                                           box(title = "Manual locations", width = 6,
                                               DT::dataTableOutput('manual_add'),
                                               actionButton("clearManual", label = "Clear all manual", width = "100%")
                                           )
                                         ),
                                         br(),
                                         fluidRow(
                                           column(4,
                                                  disabled(
                                                    actionButton("startTest", label = "Start!", width = "100%")
                                                  )
                                           ),
                                           column(4,
                                                  disabled(
                                                    actionButton("pauseTest", label = "Pause", width = "100%")
                                                  )
                                           ),
                                           column(4,
                                                  disabled(
                                                    actionButton("stopTest", label = "Stop", width = "100%")
                                                  )
                                           )
                                         ),
                                         br(),
                                         progressBar(id = "pb1", value = 0, total = 100, title = "   ", striped = TRUE),
                                         dataTableOutput('showResults')
                                         
                                ),
                                tabPanel("Structural strategy (10-2)", value = "struc_Setting",
                                         br(),
                                         fluidRow(
                                           column(3,
                                                  selectInput("segmListStrat", label = NULL, choices = "None")
                                           ),
                                           column(5,
                                                  actionButton("calcStruct", label = "Calculate structural data", width = "100%")
                                           )
                                         ),
                                         br(),
                                         fluidRow(
                                           column(3,
                                                  plotOutput("previewOCTcolorbarStrat")
                                           ),
                                           column(9,
                                                  plotOutput("previewOCTStrat") %>% withSpinner(color="#0dc5c1")
                                           )
                                         ),
                                         fluidRow(
                                           dataTableOutput('showStruct')
                                         )
                                )
                    )
                ),
                box(title = paste("Compass fundus image (60°)", sep = ""), width = 6, style = "font-size: 18px",
                    fluidRow(
                      column(6,
                             selectInput("segmListCMPTest", label = NULL, choices = "None")
                      ),
                      column(6,
                             disabled(
                               checkboxInput("plotDrasdo", label = "RGC Displacement", value = FALSE)
                             ),
                             
                             disabled(
                               checkboxInput("useOCTONH", label = "Use OCT ONH location", value = FALSE)
                             )
                      )
                    ),
                    plotOutput("fundusImageTest", width = "100%", height = "600px",
                               click = "plot_click_add",
                               dblclick = "plot_dblclick_add"),
                    plotOutput("fundusImagecolorbarTest", width = "100%", height = "70px")
                )
              )
      )
      
    ),
    tags$footer("This application is intended for informational, 
                educational and research purposes only. It is not, and is not intended, 
                for the use in diagnosis of disease and other conditions. Health care providers 
                should exercise their own independent clinical judgement when using the application 
                in conjunction with patient care.", align = "center", style = "
                width:100%;
                height:100%;   /* Height of the footer */
                color: white;
                padding: 10px;
                background-color: black;
                z-index: 1000;")
  )
)

