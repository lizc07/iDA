#' @include global.R
############
# ui
#########

#' @import shinydashboard
# header ---------------------------------------------------------------------------------------------------------------------
header  <- dashboardHeader(title = "iDA", disable = F)
# sidebar ---------------------------------------------------------------------------------------------------------------------
sidebar <- dashboardSidebar(
            #sidebarSearchForm(label = "Enter a number", "searchText", "searchButton"),
            sidebarMenu(
              # Setting id makes input$tabs give the tabName of currently-selected tab
              id = "tabs",
              menuItem("Load Data", tabName = "fileInput", icon = icon("dashboard"),badgeLabel = "Step1"),
              menuItem("Process Data", icon = icon("th"), tabName = "dataProcess", badgeLabel = "Step2",
                       badgeColor = "green"),
              menuItem("Advanced Options", icon = icon("bar-chart-o"),badgeLabel = "Step3",
                       badgeColor = "green"),

              fluidRow(column(12,
                              uiOutput("showCategory"),
                              uiOutput("showPlotOption.color"),
                              uiOutput("showPlotOption.shape"),
                              sliderInput(inputId = "sliderInput.plot.height", label = "Plot Height",
                                        min = 0, max = 5000, value = 500, step = 10, round = T, width = "100%")),
                       column(12, align = "center", hr(),downloadButton("pbmc.save", "Save Workspace"))
                      )
            ))
# body ---------------------------------------------------------------------------------------------------------------------
body <- dashboardBody(
  tabItems(
    tabItem(tabName = "fileInput",
            tabBox(title = "Loading Files", id = "tabBox.fileInput", width = 12,
                   tabPanel(title = "CSV Files",
                            fluidRow(
                                    column(12,
                                          column(10,
                                            fileInput("fileInput.csv.umi", "Expression Data (UMI count prefered, Gene x Cell) - Required!", width = "100%",
                                                      accept = c(
                                                      "text/csv",
                                                      "text/comma-separated-values,text/plain",
                                                      ".csv")
                                                      )
                                                ),
                                          br(),column(2, checkboxInput(inputId = "fileInput.csv.umi.transpose", "Do Transpose", value = F))
                                          ),
                                    column(12,
                                          column(10,
                                                 fileInput("fileInput.csv.annot", "Cell Annotation (Cell x Category)", width = "100%",
                                                           accept = c(
                                                             "text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv")
                                                          )
                                                ),
                                          br(),column(2, checkboxInput(inputId = "fileInput.csv.annot.transpose", "Do Transpose", value = F))
                                          )
                                   )
                           ),
                   tabPanel(title = "RData File",
                            fileInput("fileInput.rdata", "Saved RData, including Expression and Annotation Data (Saved by iDA)", width = "100%",
                                      accept = c(".RData")
                                      )
                           )
            ),
            box(actionButton(inputId = "fileInput.action", width = "100%", label = "Submit"), status = "success", width = 12),
            infoBoxOutput("infoBox.umi", width = 6),
            infoBoxOutput("infoBox.annot", width = 6),
            uiOutput("fileInput.error"),
            uiOutput("fileInput.submitted")
    ),

    tabItem(tabName = "dataProcess",
      tabBox(title = "Start Analysis", width = 12,
             #footer = "Produced by Zongcheng Li, using shiny", status = "primary", solidHeader = T, collapsible = F, collapsed = F,
             # QC -------------------------------------------------------------------------
             tabPanel(title = "Quality Control",

                      fluidRow(
                        column(6,
                               h4("Gene Number in each Cell"),
                               column(6,numericInput(inputId = "numericInput.seurat.gene.min",
                                                     label =  "Minimum", value = 500, width = "100%")),
                               column(6,numericInput(inputId = "numericInput.seurat.gene.max",
                                                     label =  "Maximum", value = "Inf", width = "100%"))
                        ),
                        column(6,
                               h4("UMI Number in each Cell"),
                               column(6,numericInput(inputId = "numericInput.seurat.umi.min",
                                                     label =  "Minimum", value = 10000,width = "100%")),
                               column(6,numericInput(inputId = "numericInput.seurat.umi.max",
                                                     label =  "Maximum", value = "Inf",width = "100%"))
                        ),
                        column(12,
                               actionButton(inputId = "dataProcess.qc.action", width = "100%", label = "Apply QC"),

                               column(9, uiOutput("showPlot.qc"))
                        )
                      )),
             # HVGs -------------------------------------------------------------------------
             tabPanel(title = "Identify variable genes",

                      fluidRow(
                        column(6,
                               h4("Expression Cutoff"),
                               column(6,numericInput(inputId = "numericInput.seurat.hvg.x.min",
                                                     label =  "Minimum", value = 1, width = "100%")),
                               column(6,numericInput(inputId = "numericInput.seurat.hvg.x.max",
                                                     label =  "Maximum", value = "Inf", width = "100%"))
                        ),
                        column(6,
                               h4("Dispersion Cutoff"),
                               column(6,numericInput(inputId = "numericInput.seurat.hvg.y.min",
                                                     label =  "Minimum", value = 1,width = "100%")),
                               column(6,numericInput(inputId = "numericInput.seurat.hvg.y.max",
                                                     label =  "Maximum", value = "Inf",width = "100%"))
                        ),
                        column(12,
                               actionButton(inputId = "dataProcess.hvg.action", width = "100%", label = "Get Variable Genes"),
                               uiOutput("showPlot.hvg")
                        )
                      )),
             # PCA tSNE ------------------------------------------------------------------------
             tabPanel(title = "Dimension Reduction",
                      fluidRow(

                        box(width = 6, title = "PCA Analysis",
                               column(6,
                               selectInput(inputId = "selectInput.seurat.pca", label = "Genes to be used", width = "100%",
                                           choices = c("All Genes", "HVGs", "Selected Genes"), selected = "HVGs")),
                               column(6,
                               textInput(inputId = "textInput.seurat.pca", label = "Selected Genes (If Applicable)",
                                         width = "100%", placeholder = "Input selected genes")),
                               actionButton(inputId = "dataProcess.pca.action", width = "100%", label = "Get PCA Plot"),
                               uiOutput("showPlot.pca")
                        ),
                        box(width = 6, title = "tSNE Analysis",
                               column(4,
                               selectInput(inputId = "selectInput.seurat.tsne.method", label = "tSNE Method",width = "100%",
                                               choices = c("Seurat", "Rtsne"), selected = "Seurat")),
                               column(4,
                               selectInput(inputId = "selectInput.seurat.tsne", label = "Features to be used",width = "100%",
                                           choices = c("Selected PCs", "HVGs", "Selected Genes"), selected = "HVGs")),
                               column(4,
                               textInput(inputId = "textInput.seurat.tsne", label = "Select Features:",width = "100%",
                                         placeholder = "Input selected genes or PCs (such as 1,2,4 or 1-3)")),
                               column(4,
                                      sliderInput(inputId = "sliderInput.seurat.tsne.perplexity", label = "perplexity",
                                           min = 0, max = 150, value = 30, round = T, width = "100%")),
                               column(4,
                                      sliderInput(inputId = "sliderInput.seurat.tsne.theta", label = "theta",
                                           min = 0, max = 1, value = 0.5, round = F, width = "100%")),
                               column(4,
                                      sliderInput(inputId = "sliderInput.seurat.tsne.max_iter", label = "iterations",
                                           min = 100, max = 2500, value = 1000, round = T, width = "100%")),
                               actionButton(inputId = "dataProcess.tsne.action", width = "100%", label = "Get tSNE Plot"),
                               uiOutput("showPlot.tsne")
                        )

                      )),
             # Clustering -------------------------------------------------------------------------
             tabPanel(title = "Clustering",
                      fluidRow(
                        box(width = 6, title = "Find Clusters",
                            column(6,
                                   selectInput(inputId = "selectInput.seurat.cluster", label = "Features to be used",width = "100%",
                                               choices = c("Selected PCs", "tSNEs", "HVGs", "Selected Genes"), selected = "HVGs")),
                            column(6,
                                   textInput(inputId = "textInput.seurat.cluster.feature", label = "Select Features:",width = "100%",
                                             placeholder = "Input selected genes or PCs (such as 1,2,4 or 1-3)")),
                            column(4,
                                   sliderInput(inputId = "sliderInput.seurat.cluster.resolution", label = "resolution",
                                               min = 0, max = 10, value = 0.8, step = 0.1, round = F, width = "100%")),
                            column(4,
                                   sliderInput(inputId = "sliderInput.seurat.cluster.k.param", label = "k.param",
                                               min = 0, max = 150, value = 30, round = T, width = "100%")),
                            column(4,
                                   sliderInput(inputId = "sliderInput.seurat.cluster.k.scale", label = "k.scale",
                                               min = 0, max = 150, value = 30, round = T, width = "100%")),

                            actionButton(inputId = "dataProcess.cluster.action", width = "100%", label = "Get Clusters"),
                            uiOutput("showPlot.cluster.tsne")
                        ),
                        tabBox(width = 6,
                               tabPanel(title = "Differential Expression Genes",
                                        column(6,
                                               uiOutput("showCluster.deg.cat")),
                                        # column(6,
                                        #        checkboxInput(inputId = "checkboxInput.seurat.deg.scale", label = "Do Scale",
                                        #                      width = "100%", value =F)),
                                        column(6,
                                               sliderInput(inputId = "sliderInput.seurat.deg.pct", label = "min.pct",
                                                           min = 0, max = 1, value = 0.1, step = 0.01, round = F, width = "100%")),
                                        column(6,
                                               sliderInput(inputId = "sliderInput.seurat.deg.th", label = "thresh.use",
                                                           min = 0, max = 10, value = 1, step = 0.01, round = F, width = "100%")),
                                        column(6,
                                               sliderInput(inputId = "sliderInput.seurat.deg.num", label = "DEGs Num per Group",
                                                           min = 1, max = 100, value = 10, round = T, width = "100%")),

                                        actionButton(inputId = "dataProcess.cluster.deg.action", width = "100%", label = "Get Heatmap"),
                                        uiOutput("showPlot.cluster.heatmap")
                                        ),
                               tabPanel(title = "Lightening Plot",
                                        column(6,
                                               textInput(inputId = "textInput.seurat.light.gene", label = "Input gene(s)",
                                                             width = "100%", placeholder = "Input gene")),
                                        actionButton(inputId = "dataProcess.cluster.light.action", width = "100%", label = "Get LighteningPlot"),
                                        uiOutput("showPlot.cluster.light")
                               )

                        )
                      ))

          )
    ),

    # Cluster -------------------------------------------------------------------------
    tabItem(tabName = "tSNE and Clustering"
    ),
    tabItem(tabName = "pointPlot"
      # DEG
    )

  )
)
# page ---------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(title = "Interactive Data Analysis", skin = "purple",
              header = header, sidebar = sidebar, body = body)
