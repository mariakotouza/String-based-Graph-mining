
shinyUI(
  
  fluidPage(theme=shinytheme(theme="flatly"),
            
            navbarPage("Biograph",
                     
                       tabPanel("Graph Creation",
                                
                                
                                fluidRow(
                                  add_busy_bar(color = "#FF0000"),
                                  #add_busy_gif(src = "https://media.giphy.com/media/tXL4FHPSnVJ0A/giphy.gif", height = 120, width = 120),
                                  
                                  titlePanel(""),
                                  sidebarLayout(
                                    sidebarPanel(width=2,
                                                 shinyFilesButton("files","Insert a file","Insert a file",multiple=FALSE),
                                                 #fileInput("file","Upload the file"), # fileinput() function is used to get the file upload contorl option
                                                 # helpText("Default max. file size is 5MB"),
                                                 tags$hr(),
                                                 h5(helpText("Select the read.table parameters below")),
                                                 checkboxInput(inputId = 'header', label = 'Header', value = TRUE),
                                                 br(),
                                                 radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = '\t'),
                                                 br(),br(),br(),
                                                 br(),
                                                 selectInput("seqSelect","Select the column of the sequence",choices=sequences_for_distance_calc,selected="Summary.AA.JUNCTION",width=350),
                                                 selectInput("simSelect","Select the similarity metric",choices=list("OSA"="osa",
                                                                                                                     "Levenshtein"= "lv",
                                                                                                                     "Full DL"="dl", 
                                                                                                                     "Hamming"="hamming",
                                                                                                                     "LCS"="lcs",
                                                                                                                     "Q-Gram"="qgram",
                                                                                                                     "Cosine"="cosine",
                                                                                                                     "Jaccard"="jaccard",
                                                                                                                     "Jaro Winklar"="jw"),selected="lv",width=350),
                                                 
                                                 useShinyalert(),
                                                 actionButton("simButton","Calculate distance matrix"),
                                                 
                                                 br(),
                                                 br(),
                                                 selectInput("sparsification","Select the sparsification method",choices=c("knn","enn", "mknn"),selected="enn",width=350),
                                                 
                                                         
                                                 checkboxInput(inputId = "delete_edges_of_diff_classes", label = "Delete edges that connect nodes of different classes", value = FALSE),
                                                 
                                                 conditionalPanel(
                                                   condition = "input.sparsification == 'enn'",
                                                   sliderInput("slider", label = h3("Edge Threshold"), min = 0, max = 1, value = 0.25,step=0.001,width=500)
                                                 ),
                                                 
                                                 conditionalPanel(
                                                   condition = "input.sparsification == 'knn'",
                                                   numericInput("k_knn", "Number of neighbors", 10,  min = 0, max = 100000000, width="140px")
                                                 ),
                                                 
                                                 conditionalPanel(
                                                   condition = "input.sparsification == 'mknn'",
                                                   numericInput("m_k_knn", "Number of neighbors", 10,  min = 0, max = 100000000, width="140px")
                                                 ),
                                                 
                                                 checkboxInput(inputId = "clusterId", label = "Unique sequence-clusterID combination", value = FALSE),
                                                 br(),
                                                 
                                                 actionButton("graphButton","Create Graph"),
                                                 br(),br(),
                                                 selectInput("componentSelect","Components",choices=list("no components"="")),
                                                 actionButton("componentButton","Select a component")),
                                    
                                    
                                    mainPanel(
                                      
                                      tabsetPanel(
                                        tabPanel("Graph",
                                                 column(2,numericInput("nodeSelect","Insert an id",0)),
                                                 column(2,br(),verbatimTextOutput("graphMeasure")),
                                                 
                                                 column(7,
                                                        div(style = "position:absolute;right:2em;", 
                                                            br(),
                                                            column(2,actionButton("getNeigh","Neighbours")),
                                                            column(2,downloadButton("downloadGraph", "Download Graph"),offset=1),
                                                            column(2,offset=2,actionButton("visGraph", "Visualise Graph"))
                                                        )),
                                                 br(),br(),br(),br(),br(),
                                                 visNetworkOutput("network",height = 800,width=1500),
                                                 bsModal("modal1", "Neighbours", "getNeigh", size = "large",
                                                         dataTableOutput("neigh"))
                                                 
                                        ),
                                        tabPanel("Data",
                                                 numericInput("idInput","Insert the id",0), 
                                                 dataTableOutput("dataset")
                                                 #,
                                                 #downloadButton("downloadGraph", "Download Graph")
                                        ),
                                        
                                        tabPanel("Adjacency",
                                                 dataTableOutput("adj")
                                        ),
                                        
                                        tabPanel("Unique Combo",
                                                 dataTableOutput("uni")
                                        )
                                        
                                      )
                                      
                                      # use below code if you want the tabset programming in the main panel. If so, then tabset will appear when the app loads for the first time.
                                      #       tabsetPanel(tabPanel("Summary", verbatimTextOutput("sum")),
                                      #                   tabPanel("Data", tableOutput("table")))
                                    ))),
                                
                                h3("Filters") ,   
                                column(5,
                                       
                                       fluidRow(style="background-color:rgb(238, 239, 239);border-radius:5px;border-color=LightGrey;border-color=Black;border-width=thick;",
                                                
                                                
                                                
                                                
                                                br(),
                                                fluidRow(
                                                  column(4,
                                                         uiOutput("selectbox1"),
                                                         uiOutput("textbox1"),
                                                         br(),
                                                         column(1,actionButton("includeButton","Include"),br(),br(),actionButton("FilterButton","Filter")),
                                                         column(1, actionButton("excludeButton","Exclude"),br(),br(),actionButton("ResetButton","Reset"),offset=2),
                                                         column(1,br(),br(),checkboxInput("ReverseButton","Reverse Filters"),offset=2)
                                                  ),
                                                  
                                                  column(8,
                                                         h4("Filters Used"),
                                                         tableOutput("filtertable")
                                                  )
                                                ),
                                                br(),
                                                
                                                hr("Filtered Indexes"),
                                                verbatimTextOutput("Indexes"))
                                )
                       ),

                       tabPanel("MST",
                                fluidRow(
                                  
                                  #add_busy_gif(src = "https://media.giphy.com/media/tXL4FHPSnVJ0A/giphy.gif", height = 120, width = 120),
                                  
                                  column(10,
                                         div(style = "position:absolute;right:1em;", 
                                             column(2,actionButton("getNeigh2","Neighbours")),
                                             column(2,downloadButton("downloadMST", "Download MST"),offset=1),
                                             column(2,offset=2,actionButton("visMST", "Visualise MST"))
                                             
                                         ),
                                         selectInput("mstAlgoSelect",choices=list("Prim"="Prim","Kruskal"="Kruskal"),selected="Prim",label="Select MST algorithm"),
                                         numericInput("nodeSelect2","Insert an id",0),
                                         visNetworkOutput("mstnetwork",height=1000),
                                         bsModal("modal2", "Neighbours", "getNeigh2", size = "large",
                                                 dataTableOutput("neigh2"))),
                                  column(2,
                                         br(),br(),
                                         plotOutput("mstLegend",width=350),
                                         plotOutput("mstLegend2",width = 350)),  
                                  selectInput("colormst","Select an attribute for background coloring",choices=c(Default="Default"),selected = "Default"),
                                  selectInput("bordermst","Select an attribute for border coloring",choices=c(Default="Default"),selected = "Default"),
                                  checkboxInput("clusterMST","Coloring according to clustering",value=FALSE),
                                  actionButton("mstButton","Create MST"),
                                  hr("Central Nodes"),
                                  verbatimTextOutput("Cendroids"),
                                  br(),
                                  dataTableOutput("adj2")
                                  
                                )
                                
                       ),
                       
                       tabPanel("Centralities",
                                fluidRow(
                                  
                                  
                                  actionButton("centralButton","Calculate Centralities"),
                                  checkboxInput("fastCButton","Only major centralities metrics"),
                                  div(style = "position:absolute;right:1em;", 
                                      downloadButton("downloadCentral", "Download Centralities")),
                                  br(),
                                  
                                  tabsetPanel(
                                    
                                    tabPanel("Centralities",DT::dataTableOutput("Centralities")),
                                    tabPanel("Rankings",DT::dataTableOutput("Rankings")),
                                    tabPanel("Summary",DT::dataTableOutput("Summary")),
                                    tabPanel("Network",
                                             column(10,
                                                    selectInput("centralSelect","Select a centrality type",list('Degree'="Degree.Centrality",
                                                                                                                'Betweenness'="Shortest.Paths.Betweenness.Centrality",
                                                                                                                'Closeness'="Average.Distance",
                                                                                                                'Eigenvector'="eigenvector.centralities"),selected="Shortest.Paths.Betweenness.Centrality"),
                                                    numericInput("nodeSelect3","Insert an id",0),                                                          
                                                    visNetworkOutput("centralnetwork",height=800),
                                                    hr("Note:Triangle indicates the most 'central' code.")),
                                             column(2,
                                                    
                                                    checkboxInput("clusterCentral","Coloring according to clustering",value=FALSE),
                                                    selectInput("centralColor","Select an attribute for nodes color",choice=list("no attr"="")),plotOutput("centralLegend",width=300))
                                             
                                    )           
                                  )
                                )
                       ),
                       
                       
                       tabPanel("Clustering",
                                fluidRow(   
                                  
                                  
                                  tabsetPanel(
                                    tabPanel("Clustering Analysis",
                                             tags$head(
                                               tags$style(
                                                 HTML(".shiny-notification {
                                                      
                                                      position:fixed;
                                                      top: calc(40%);
                                                      left: calc(10%);
                                                      }"
                                                      )
                                                 )
                                                 ),
                                             
                                             selectInput("clusterSelect","Select a clustering algorithm",list("no algorithm"=" ",
                                                                                                              "Louvain"="louvain",
                                                                                                              "Fast Greedy"="fast_greedy",
                                                                                                              "Label Propagation"="label_propagation",
                                                                                                              "Leading Eigenvalue"="leading_eigenvalue",
                                                                                                              "Walktrap"="walktrap",
                                                                                                              "Edge Betweeness"="edge_betweenness",
                                                                                                              "Hierarchical"="hierarchical"),selected = " "),
                                             
                                             fluidRow(
                                               column(8,
                                                      div(style = "position:absolute;right:1em;",
                                                          column(2,downloadButton("downloadCluster","Download Membership")),
                                                          column(2,offset=4,actionButton("visCluster", "Visualise Clusters"))
                                                          
                                                          
                                                      ),   
                                                      numericInput("nodeSelect4","Insert an id",0),
                                                      visNetworkOutput("clusterNetwork",height=800)),
                                               column(2,
                                                      selectInput("clusterColor","Select an attribute for nodes color",choice=list("no attr"="")),plotOutput("clusterLegend",width=300)),
                                               column(2,
                                                      checkboxInput("plotHeat","Plot a heatmap",value=FALSE),
                                                      plotOutput("intraHeatmap"))
                                             ),
                                             
                                             fluidRow(
                                               
                                               column(4,plotOutput("silhouette")),
                                               column(4,br(),br(),br(),verbatimTextOutput("metrics"),offset=2)
                                               
                                               
                                               
                                             )
                                                 ),

                                    tabPanel("Clustering Crossover",
                                             selectInput("clusterSelect1","Select a clustering algorithm",list("Louvain"="louvain",
                                                                                                               "Fast Greedy"="fast_greedy",
                                                                                                               "Label Propagation"="label_propagation",
                                                                                                               "Leading Eigenvalue"="leading_eigenvalue",
                                                                                                               "Walktrap"="walktrap",
                                                                                                               "Edge Betweeness"="edge_betweenness",
                                                                                                               "Hierarchical"="hierarchical"),selected = "louvain"),
                                             selectInput("clusterSelect2","Select a clustering algorithm",list("Louvain"="louvain",
                                                                                                               "Fast Greedy"="fast_greedy",
                                                                                                               "Label Propagation"="label_propagation",
                                                                                                               "Leading Eigenvalue"="leading_eigenvalue",
                                                                                                               "Walktrap"="walktrap",
                                                                                                               "Edge Betweeness"="edge_betweenness",
                                                                                                               "Hierarchical"="hierarchical"),selected = "louvain"),
                                             
                                             fluidRow(
                                               column(2,
                                                      DT::dataTableOutput("confusionMatrix1")),
                                               column(2,
                                                      DT::dataTableOutput("confusionMatrix2"),offset=1),
                                               
                                               column(4,
                                                      verbatimTextOutput("metrics2"),
                                                      br(),br(),br(),br(),
                                                      selectInput("heatSelect",label="",choices=list("Overall"=0,"Conditional under clustering 1"=1,"Conditional under clustering 2"=2),selected=0),
                                                      plotOutput("interHeatmap"),offset=1)
                                               
                                             )
                                             
                                             
                                    )
                                               )
                                  )
                       ),
                       
                       tabPanel("SSL",
                                fluidRow(
                                  
                                  #add_busy_gif(src = "https://media.giphy.com/media/tXL4FHPSnVJ0A/giphy.gif", height = 120, width = 120),
                                  
                                  column(10,
                                         div(style = "position:absolute;right:1em;", 
                                             #column(2,actionButton("getNeigh2","Neighbours")),
                                             column(2,downloadButton("downloadSSL_results", "Download SSL Results"),offset=1)
                                             
                                         ),
                                         selectInput("sslAlgoSelect",choices=list("Label_prop"="Label propagation"),selected="Label_prop",label="Select the SSL algorithm"),
                                         
                                         uiOutput("uiSSLColNames"),
                                         uiOutput("uiSSLRefNames"),
                                         actionButton("runSSLBtn","Run SSL"),
                                         conditionalPanel(
                                           condition = "input.runSSLBtn",
                                           br(),
                                           DT::dataTableOutput("sslResultsMatrix"),
                                           br(),
                                           DT::dataTableOutput("sslEvaluation")
                                         ),
                                         
                                         
                                         #numericInput("nodeSelect2","Insert an id",0),
                                         #visNetworkOutput("mstnetwork",height=1000),
                                         #bsModal("modal2", "Neighbours", "getNeigh2", size = "large",
                                        #         dataTableOutput("neigh2"))),
                                  column(2,
                                         br(),br()#,
                                         #plotOutput("mstLegend",width=350),
                                         #plotOutput("mstLegend2",width = 350)),  
                                  #selectInput("colormst","Select an attribute for background coloring",choices=c(Default="Default"),selected = "Default"),
                                  #selectInput("bordermst","Select an attribute for border coloring",choices=c(Default="Default"),selected = "Default"),
                                  #checkboxInput("clusterMST","Coloring according to clustering",value=FALSE),
                                  #actionButton("mstButton","Create MST"),
                                  #hr("Central Nodes"),
                                  #verbatimTextOutput("Cendroids"),
                                  #br(),
                                  #dataTableOutput("adj2")
                                  
                                )
                                
                       )
                       
  )
  )
  
  )



)
)












