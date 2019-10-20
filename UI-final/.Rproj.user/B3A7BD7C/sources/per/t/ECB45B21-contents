library(shiny)
shinyUI(navbarPage("Navbar",
           tabPanel("Graph Creation",
        
        
                  fluidRow(
            
            
                    titlePanel("File Input"),
                    sidebarLayout(
                      sidebarPanel(
                        fileInput("file","Upload the file"), # fileinput() function is used to get the file upload contorl option
                        helpText("Default max. file size is 5MB"),
                        tags$hr(),
                        h5(helpText("Select the read.table parameters below")),
                        checkboxInput(inputId = 'header', label = 'Header', value = TRUE),
                        br(),
                        radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = ''),
                        br(),br(),br()
                        
                      ),
                       mainPanel(
                              
                                tabsetPanel(
                                    tabPanel("Graph",visNetworkOutput("network")),
                                    tabPanel("Data",numericInput("idInput","Insert the id",0),
                                                    dataTableOutput("dataset"))
                                )
                        # use below code if you want the tabset programming in the main panel. If so, then tabset will appear when the app loads for the first time.
                        #       tabsetPanel(tabPanel("Summary", verbatimTextOutput("sum")),
                        #                   tabPanel("Data", tableOutput("table")))
                      )
                      
                    )
                  ),
                 fluidRow(
                   br(),
                    sliderInput("slider", label = h3("Edge Threshold"), min = 0, 
                               max = 1, value = 0.6),
                   br(),
                    actionButton("graphButton","Create Graph")
                 
                 
                 
               )   
          
          
          ),
          tabPanel("GraphFiltering",
                   fluidRow(
                     column(2,
                        hr("Include"),
                        uiOutput("selectbox1"),
                        uiOutput("textbox1"),
                        actionButton("includeButton","Add to graph")
                        ),
                     column(4,
                        hr("Exclude"),
                        uiOutput("selectbox2"),
                        uiOutput("textbox2"),
                        actionButton("excludeButton","Delete from graph")
                        ),
                     column(6,
                       hr("Filters Added"),
                       tableOutput("filtertable")
                            )
                   ),
                   br(),br(),
                   fluidRow(
                      column(1,
                         actionButton("FilterButton","Filter")),
                      column(2,
                         actionButton("ResetButton","Reset")),
                      column(3,
                         checkboxInput("ReverseButton","Reverse Filters"))
                            ),
                   verbatimTextOutput("Indexes")
                    
                    ),
          
          
          tabPanel("MST",
                   
                   visNetworkOutput("mstnetwork"),
                   selectInput("colormst","Select an attribute for background coloring",choices=c(Default="Default"),selected = "Default"),
                   selectInput("bordermst","Select an attribute for border coloring",choices=c(Default="Default"),selected = "Default"),
                   actionButton("mstButton","Create MST"),
                   verbatimTextOutput("Cendroids"),
                   verbatimTextOutput("Links")
            
            
          ),
          
          
          tabPanel("Centralities",
                   actionButton("centralButton","Calculate Centralities"),
                   br(),br(),
                   tabsetPanel(
                     
                          tabPanel("Centralities",DT::dataTableOutput("Centralities")),
                          tabPanel("Rankings",DT::dataTableOutput("Rankings")),
                          tabPanel("Summary",DT::dataTableOutput("Summary"))
                     
                        
                     
                     
                     
                     
                   )
                   
                   
                   
                   
          )
))
                
              
        
                      
                      
            
            
            
            
          
          
          
          
          
          
      
