############
# server.R
############
#' @include global.R


server <- function(input, output){
  ############
  # Load Data
  ############
  # dataInput
  dataInput <- reactive({
    umi <- annot <- color <- para <- NULL
    if(input$tabBox.fileInput == "CSV Files"){
      if(!is.null(input$fileInput.csv.umi)) {
        withProgress(message = "Loading in progress...",
                     tmp <- try(read.csv(file = input$fileInput.csv.umi$datapath[1], header = T, row.names = 1))
        )
        if(class(tmp) == "try-error"){
          output$fileInput.error <- renderUI({box(renderText({tmp}),width = 12, status = "danger",solidHeader=T,
                                                  title = "Error", footer = "Please Contact lizc07@vip.qq.com")
          })
        }else{
          umi <- tmp
          if(input$fileInput.csv.umi.transpose) umi <- as.data.frame(t(umi))
        }


      }
      if(!is.null(input$fileInput.csv.annot)) {
        withProgress(message = "Loading in progress...",
                     tmp <- try(read.csv(file = input$fileInput.csv.annot$datapath[1], header = T, row.names = 1))
        )
        if(class(tmp) == "try-error"){
          output$fileInput.error <- renderUI({box(renderText({tmp}),width = 12, status = "danger",solidHeader=T,
                                                  title = "Error", footer = "Please Contact lizc07@vip.qq.com")
          })
        }else{
          annot <- tmp
          if(input$fileInput.csv.annot.transpose) annot <- as.data.frame(t(annot))
        }
      }
    } else if(input$tabBox.fileInput == "RData File"){
      if(!is.null(input$fileInput.rdata)) {
        tmp <- try(load(file = input$fileInput.rdata$datapath[1]))
        #print(pbmc)
        if(class(tmp) == "try-error"){
          output$fileInput.error <- renderUI({box(renderText({tmp}),width = 8, status = "danger",solidHeader=T,
                                                  title = "Error", footer = "Please Contact lizc07@vip.qq.com")
          })
        }else{
          pbmc <<- pbmc
          umi <- pbmc@raw.data
          annot <- pbmc@meta.data
        }
      }
    }
    return(list(umi = umi, annot = annot, color = color, para = para))
  })

  dataSubmit <- reactive({
    input$fileInput.action
    rm(list = ls())
    if(is.null(input$fileInput.error)){
      isolate(dataInput())
    }
  })

  # output$infoBox
  upload_data <- observeEvent(input$fileInput.action, {
    output$infoBox.umi <- renderInfoBox({
      if(checkdata.showModal(dataSubmit()$umi, "UMI Data")){
        if(is.null(dataSubmit()$annot)){
          text.submitted <- "Warning: Annotation Data is not uploaded!"
          status.submitted <- "warning"
        }else if(!identical(sort(rownames(dataSubmit()$annot)), sort(colnames(dataSubmit()$umi)))){
          text.submitted <- "Warning: Cell Numbler/Names is not identical!"
          status.submitted <- "warning"
        }else{
          text.submitted <- "Please continue to process"
          status.submitted <- "success"
        }
        output$fileInput.submitted <- renderUI({box(text.submitted, width = 12, status = status.submitted, solidHeader=T,
                                                    title = "Submitted", footer = "Produced by Zongcheng Li, using shiny")
        })
        infoBox(title = "Expression", value = paste0("Gene: ",nrow(dataSubmit()$umi),"\n\rCell: ", ncol(dataSubmit()$umi)), fill = F)
      }else{
        infoBox(title = "Expression", value = NULL,fill = T,subtitle="Please Load Data")
      }
    })

    output$infoBox.annot <- renderInfoBox({
      if(!is.null(dataSubmit()$annot)){
        checkdata.showModal(dataSubmit()$annot, "Annotation Data")
        infoBox(title = "Annotation", value = paste0("Cell: ",nrow(dataSubmit()$annot),"\n\rCategory: ", ncol(dataSubmit()$annot)), fill = F)
      }else{
        infoBox(title = "Annotation", value = NULL,fill = T,subtitle="Please Load Data")
      }
    })
    # check box
    output$showCategory <- renderUI({
      if(checkdata.showModal(dataSubmit()$annot, show = F)){
        annot <- formatAnnot(dataSubmit()$annot)
        return(getAnnotUI(annot))
      }else{
        return()
      }
    })

    output$showPlotOption.color <- renderUI({
      selectInput(inputId = "selectInput.plot.color.by", label = "Color By",
                  choices = colnames(dataSubmit()$annot), selected = dataSubmit()$annot[1])
    })
    output$showPlotOption.shape <- renderUI({
      selectInput(inputId = "selectInput.plot.shape.by", label = "Shape By",
                  choices = colnames(dataSubmit()$annot), selected = dataSubmit()$annot[1])
    })
  })




  # Process Data
  output$pbmc.save <- downloadHandler(filename = paste0("iDA", format(Sys.time(),'_%Y%m%d_%H%M%S'),".RData"),
                                        content = function(file) {save(pbmc, file = file)})

  # QC
  seurat_QC <- observeEvent(eventExpr = input$dataProcess.qc.action,
                            handlerExpr = {
                              umi <- dataSubmit()$umi
                              annot <- dataSubmit()$annot
                              if(!checkdata.showModal(umi, message = "Something goes wrong with UMI Data, Please Check UMI Data or Contact lizc07@vip.qq.com")){
                                return()
                              }

                              cell.genenum <- colSums(umi > 0)
                              cell.uminum <- colSums(umi)

                              cell_select <- colnames(umi)[which(cell.genenum > return.numeric(input$numericInput.seurat.gene.min, 0) &
                                                                   cell.genenum < return.numeric(input$numericInput.seurat.gene.max, Inf) &
                                                                   cell.uminum > return.numeric(input$numericInput.seurat.umi.min, 0) &
                                                                   cell.uminum < return.numeric(input$numericInput.seurat.umi.max, Inf))]
                              if(is.null(cell_select) || length(cell_select) < 2) {
                                showModal(modalDialog(title = "Notice!", "Please Modify the QC Standard to include more Cells!", easyClose = T, fade = T))
                                return()
                              }
                              withProgress(message = "Calculation in progress...",
                                           pbmc <<- CreateSeuratObject(project = "iDA", raw.data = umi[,cell_select], meta.data = annot,
                                                                       normalization.method = "LogNormalize", scale.factor = 100000,
                                                                       do.scale = T, do.center = T, save.raw = T))
                              output$qc.vlnplot <- renderPlot(VlnPlot(pbmc, c("nGene", "nUMI"),nCol = 2, group.by = input$selectInput.plot.color.by, x.lab.rot = T, do.return = T))
                              output$qc.geneplot <- renderPlot(GenePlot(pbmc, "nUMI", "nGene"))
                              output$showPlot.qc <- renderUI({
                                list(column(6, plotOutput("qc.vlnplot", height = input$sliderInput.plot.height)),
                                     column(6, plotOutput("qc.geneplot", height = input$sliderInput.plot.height)))
                              })
                            })

  # HVGs
  seurat_HVGs <- observeEvent(eventExpr = input$dataProcess.hvg.action,
                              handlerExpr = {
                                if(!checkdata.showModal(pbmc, message = "Apply QC First!")){
                                  return()
                                }
                                withProgress(message = "Calculation in progress...",
                                             pbmc <<- FindVariableGenes(pbmc, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = F,
                                                                        x.low.cutoff = return.numeric(input$numericInput.seurat.hvg.x.min, 0),
                                                                        x.high.cutoff = return.numeric(input$numericInput.seurat.hvg.x.max, Inf),
                                                                        y.cutoff = return.numeric(input$numericInput.seurat.hvg.y.min, 0),
                                                                        y.high.cutoff = return.numeric(input$numericInput.seurat.hvg.y.max, Inf),
                                                                        num.bin = 20, display.progress = F, sort.results = F))
                                output$hvg.plot <- renderPlot({
                                  VariableGenePlot(pbmc, x.low.cutoff = return.numeric(input$numericInput.seurat.hvg.x.min, 0),
                                                   x.high.cutoff = return.numeric(input$numericInput.seurat.hvg.x.max, Inf),
                                                   y.cutoff = return.numeric(input$numericInput.seurat.hvg.y.min, 0),
                                                   y.high.cutoff = return.numeric(input$numericInput.seurat.hvg.y.max, Inf))
                                })
                                output$hvg.datatable <- DT::renderDataTable({
                                  datatable(data = as.matrix(pbmc@var.genes), class = 'cell-border stripe', colnames = c('ID' = 1), filter = "top",options = list(pageLength = 10))
                                })
                                output$hvg.downloadData <- downloadHandler(filename = suffixFileName("iDA", ".csv"),
                                                                           content = function(file) {
                                                                             write.csv(pbmc@var.genes, file, row.names = FALSE)})
                                output$showPlot.hvg <- renderUI({
                                  list(column(6, align = "center",
                                              plotOutput("hvg.plot", height = input$sliderInput.plot.height),
                                              downloadButton("hvg.downloadData", paste0("Download ", length(pbmc@var.genes)," HVGs"))),
                                       column(6,
                                              hr(),
                                              DT::dataTableOutput("hvg.datatable")
                                       ))
                                })
                              })

  # PCA
  seurat_PCA <- observeEvent(eventExpr = input$dataProcess.pca.action,
                             handlerExpr = {
                               if(is.null(pbmc) || is.null(pbmc@calc.params$FindVariableGenes)){
                                 showModal(modalDialog("Find HVGs First!",title = "Error", easyClose = T, fade = T))
                                 return()
                               }
                               select_gene <- switch(input$selectInput.seurat.pca,
                                                     "All Genes" = rownames(pbmc@data), HVGs = pbmc@var.genes, "Selected Genes" = textInput2genes(input$textInput.seurat.pca))
                               if(is.null(select_gene)){
                                 showModal(modalDialog("Input genes for PCA!",title = "Error", easyClose = T, fade = T))
                                 return()
                               }
                               withProgress(message = "Calculation in progress...",{
                                 pbmc <<- RunPCA(pbmc, pc.genes = select_gene, do.print = F, pcs.print = NULL, genes.print = 0)
                                 pbmc <<- ProjectPCA(pbmc, do.print = F)
                               })
                               output$pca.datatable <- DT::renderDataTable({
                                 datatable(data = pbmc@dr$pca@gene.loadings[,1:4], class = 'cell-border stripe',
                                           colnames = c('ID' = 1), filter = "top",options = list(pageLength = 10)) %>% formatRound(columns = 1:4,digits = 3)
                               })
                               output$pca.plot <- renderPlot({
                                 PCAPlot(pbmc, dim.1 = 1, dim.2 = 2, cells.use = NULL, do.return = F,
                                         group.by = input$selectInput.plot.color.by, pt.shape = input$selectInput.plot.shape.by)
                               })
                               output$pca.elbow.plot <- renderPlot({
                                 PCElbowPlot(pbmc)
                               })
                               output$showPlot.pca <- renderUI({
                                 tabBox(width = 12,
                                        tabPanel(title = "PCA plot",
                                                 plotOutput("pca.plot", height = input$sliderInput.plot.height)),
                                        tabPanel(title = "Elbow Plot",
                                                 plotOutput("pca.elbow.plot", height = input$sliderInput.plot.height)),
                                        tabPanel(title = "PCA Table",
                                                 DT::dataTableOutput("pca.datatable"))
                                 )
                               })
                             })
  # tSNE
  seurat_tSNE <- observeEvent(eventExpr = input$dataProcess.tsne.action,
                              handlerExpr = {
                                if(is.null(pbmc) || is.null(pbmc@calc.params$FindVariableGenes) || is.null(pbmc@calc.params$RunPCA)){
                                  showModal(modalDialog("Perform PCA Analysis First!",title = "Error", easyClose = T, fade = T))
                                  return()
                                }
                                select_feature <- switch(input$selectInput.seurat.tsne,
                                                         "Selected PCs" = textInput2PCs(input$textInput.seurat.tsne),
                                                         HVGs = pbmc@var.genes,
                                                         "Selected Genes" = textInput2genes(input$textInput.seurat.tsne))
                                if(is.null(select_feature)){
                                  showModal(modalDialog("Input Features for tSNE!",title = "Error", easyClose = T, fade = T))
                                  return()
                                }
                                if(is.numeric(select_feature)){
                                  select_pca <- select_feature
                                  select_gene <- NULL
                                }else{
                                  select_pca = NULL
                                  select_gene <- select_feature
                                }
                                #print(input$sliderInput.seurat.tsne.perplexity)
                                if(input$selectInput.seurat.tsne.method == "Seurat"){
                                  withProgress(message = "Calculation in progress...",
                                               tmp <- try(pbmc <<- RunTSNE(pbmc, dims.use = select_pca, genes.use = select_gene, do.fast = T,
                                                                           perplexity = input$sliderInput.seurat.tsne.perplexity,
                                                                           theta = input$sliderInput.seurat.tsne.theta,
                                                                           max_iter= input$sliderInput.seurat.tsne.max_iter, pca = F))
                                  )
                                }else{
                                  withProgress(message = "Calculation in progress...",{
                                    if(is.null(pbmc@dr$tsne)){
                                      pbmc <<- RunTSNE(pbmc, dims.use = 1:2, genes.use = NULL, do.fast = T,
                                                       perplexity = 1,
                                                       theta = 1, max_iter= 100, pca = F)
                                    }
                                    set.seed(666)
                                    d_euclidean <- stats::dist(t(pbmc@data[select_gene,]))
                                    tsne_out <- NULL
                                    tmp <- try(pbmc@dr$tsne@cell.embeddings[] <- Rtsne(d_euclidean, is_distance=TRUE, pca = F,
                                                                                        perplexity=input$sliderInput.seurat.tsne.perplexity,
                                                                                        verbose = TRUE,
                                                                                        max_iter=input$sliderInput.seurat.tsne.max_iter,
                                                                                        theta = input$sliderInput.seurat.tsne.theta)$Y[,c(1,2)]
                                    )
                                  })
                                }
                                if(class(tmp) == "try-error"){
                                  output$showPlot.tsne <- renderUI({box(renderText({tmp}),width = 12, status = "danger",solidHeader=T,
                                                                        title = "Error", footer = "Please Contact lizc07@vip.qq.com")
                                  })
                                  return()
                                }else{
                                  output$tsne.plot <- renderPlot({
                                    TSNEPlot(pbmc, cells.use = NULL, do.return = F, pt.size = 3,
                                             group.by = input$selectInput.plot.color.by, pt.shape = input$selectInput.plot.shape.by)
                                  })
                                  output$showPlot.tsne <- renderUI({
                                    list(column(12, align = "center",
                                                plotOutput("tsne.plot", height = input$sliderInput.plot.height))
                                    )
                                  })
                                }
                              })

  # Clustering - cluster
  seurat_cluster <- observeEvent(eventExpr = input$dataProcess.cluster.action,
                                 handlerExpr = {
                                   if(is.null(pbmc) || is.null(pbmc@calc.params$FindVariableGenes) || is.null(pbmc@calc.params$RunPCA) || is.null(pbmc@calc.params$RunTSNE)){
                                     showModal(modalDialog("Perform tSNE Analysis First!",title = "Error", easyClose = T, fade = T))
                                     return()
                                   }
                                   select_feature <- switch(input$selectInput.seurat.cluster,
                                                            "Selected PCs" = textInput2PCs(input$textInput.seurat.cluster.feature),
                                                            tSNEs = "tsne",
                                                            HVGs = pbmc@var.genes,
                                                            "Selected Genes" = textInput2genes(input$textInput.seurat.cluster.feature))
                                   if(is.null(select_feature)){
                                     showModal(modalDialog("Input Features for tSNE!",title = "Error", easyClose = T, fade = T))
                                     return()
                                   }
                                   if(is.numeric(select_feature)){
                                     reduction.type <- "pca"
                                     select_pca <- select_feature
                                     select_gene <- NULL
                                   }else if(select_feature == "tsne"){
                                     reduction.type <- "tsne"
                                     select_pca <- c(1,2)
                                   }else{
                                     reduction.type <- ""
                                     select_pca = NULL
                                     select_gene <- select_feature
                                   }
                                   withProgress(message = "Calculation in progress...",{
                                     pbmc@calc.params$BuildSNN <- NULL
                                     tmp <- try(pbmc <<- FindClusters(pbmc, genes.use = select_gene, reduction.type = reduction.type, dims.use = select_pca,
                                                                      force.recalc = T, reuse.SNN = F,
                                                                      save.SNN = T, plot.SNN = F, resolution = input$sliderInput.seurat.cluster.resolution,
                                                                      k.param = input$sliderInput.seurat.cluster.k.param,
                                                                      k.scale = input$sliderInput.seurat.cluster.k.scale))
                                   })
                                   if(class(tmp) == "try-error"){
                                     output$showPlot.cluster.tsne <- renderUI({box(renderText({tmp}),width = 12, status = "danger",solidHeader=T,
                                                                                   title = "Error", footer = "Please Contact lizc07@vip.qq.com")
                                     })
                                     return()
                                   }else{
                                     output$showCluster.deg.cat <- renderUI({
                                       selectInput(inputId = "selectInput.seurat.deg.cat", label = "Groups",
                                                   choices = colnames(pbmc@meta.data), selected = "orig.ident", width = "100%")
                                     })

                                     output$cluster.tsne.snn <- renderPlot({
                                       net <- graph.adjacency(adjmatrix = pbmc@snn, mode = "undirected",
                                                              weighted = TRUE, diag = FALSE)
                                       plot.igraph(x = net, layout = as.matrix(x = pbmc@dr$tsne@cell.embeddings),
                                                   edge.width = E(graph = net)$weight, vertex.label = NA,
                                                   vertex.size = 0)
                                     })
                                     output$cluster.tsne.plot <- renderPlot({
                                       TSNEPlot(pbmc, cells.use = NULL, do.return = F, pt.size = 3)
                                     })
                                     output$showPlot.cluster.tsne <- renderUI({
                                       tabBox(width = 12,
                                              tabPanel(title = "tSNE plot",
                                                       plotOutput("cluster.tsne.plot", height = input$sliderInput.plot.height)),
                                              tabPanel(title = "SNN Plot",
                                                       plotOutput("cluster.tsne.snn", height = input$sliderInput.plot.height))
                                       )
                                     })
                                   }
                                 })

  # Clustering - deg
  seurat_deg <- observeEvent(eventExpr = input$dataProcess.cluster.deg.action,
                             handlerExpr = {
                               if(is.null(pbmc) || is.null(pbmc@calc.params$FindVariableGenes) || is.null(pbmc@calc.params$RunPCA) || is.null(pbmc@calc.params$RunTSNE) || !any(grepl("FindClusters.", names(pbmc@calc.params)))){
                                 showModal(modalDialog("Perform Find Clusters First!",title = "Error", easyClose = T, fade = T))
                                 return()
                               }

                               withProgress(message = "Calculation in progress...",{
                                 tmp <- try({
                                   ### diff
                                   pbmc <<- SetAllIdent(pbmc, id = input$selectInput.seurat.deg.cat)
                                   pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = input$sliderInput.seurat.deg.pct,
                                                                   thresh.use = input$sliderInput.seurat.deg.th,
                                                                   test.use = "roc",latent.vars = NULL)
                                   pbmc.markers <- subset(pbmc.markers, myAUC >= 0.7)
                                   pbmc.markers %>% group_by(cluster) %>% top_n(input$sliderInput.seurat.deg.num, power) -> top10
                                 })
                               })
                               if(class(tmp) == "try-error"){
                                 output$showPlot.cluster.heatmap <- renderUI({box(renderText({tmp}),width = 12, status = "danger",solidHeader=T,
                                                                                  title = "Error", footer = "Please Contact lizc07@vip.qq.com")
                                 })
                                 return()
                               }else{
                                 output$cluster.deg.heatmap <- renderPlot({
                                   DoHeatmap(pbmc, genes.use = top10$gene,#col.low = "darkblue", col.mid = "white",col.high = "red",
                                             use.scaled = T, slim.col.label = TRUE, remove.key = TRUE, cex.col = 0.6)
                                 })
                                 output$cluster.deg.heatmap.u <- renderPlot({
                                   pheatmap(pbmc@data[top10$gene, order(pbmc@meta.data[,input$selectInput.seurat.deg.cat])],
                                            cluster_rows = F, cluster_cols = F, annotation_col = pbmc@meta.data[,c(input$selectInput.seurat.deg.cat), drop=F],
                                            show_colnames = F, color = pal_heatmap2, border_color = NA)
                                 })
                                 output$showPlot.cluster.heatmap <- renderUI({
                                   tabBox(width = 12,
                                          tabPanel(title = "Expression Heatmap",
                                                   plotOutput("cluster.deg.heatmap.u", height = input$sliderInput.plot.height)),
                                          tabPanel(title = "Scaled Expression Heatmap",
                                                   plotOutput("cluster.deg.heatmap", height = input$sliderInput.plot.height))
                                   )
                                 })
                               }
                             })

  # Clustering - light
  seurat_light <- observeEvent(eventExpr = input$dataProcess.cluster.light.action,
                               handlerExpr = {
                                 if(is.null(pbmc) || is.null(pbmc@calc.params$FindVariableGenes) || is.null(pbmc@calc.params$RunPCA) || is.null(pbmc@calc.params$RunTSNE)){
                                   showModal(modalDialog("Perform tSNE Analysis First!",title = "Error", easyClose = T, fade = T))
                                   return()
                                 }

                                 withProgress(message = "Calculation in progress...",{
                                   tmp <- try({
                                     myFeaturePlot(pbmc, features.plot = textInput2genes(input$textInput.seurat.light.gene), ncol =4)
                                   })
                                 })
                                 if(class(tmp) == "try-error"){
                                   output$showPlot.cluster.heatmap <- renderUI({box(renderText({tmp}),width = 12, status = "danger",solidHeader=T,
                                                                                    title = "Error", footer = "Please Contact lizc07@vip.qq.com")
                                   })
                                   return()
                                 }else{
                                   output$cluster.light.plot <- renderPlot({
                                     myFeaturePlot(pbmc, features.plot = textInput2genes(input$textInput.seurat.light.gene), ncol = NULL)
                                   })
                                   output$showPlot.cluster.light <- renderUI({
                                     plotOutput("cluster.light.plot", height = input$sliderInput.plot.height)
                                   })
                                 }
                               })
}
