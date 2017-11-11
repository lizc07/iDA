options(stringsAsFactors = F)
rm(list = ls())


#' @import DT
#' @import shiny
#' @import shinydashboard
#' @import Seurat
#' @import Rtsne
#' @importFrom igraph E graph.adjacency plot.igraph
#' @importFrom dplyr %>% group_by top_n
#' @import pheatmap
#' @importFrom gplots colorpanel
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
NULL


myFeaturePlot <- function(object, features.plot, nrow = NULL, ncol = NULL, ...){
  #require(ggplot2)
  #require(gridExtra)
  ggData <- as.data.frame(cbind(object@dr$tsne@cell.embeddings,FetchData(object, features.plot)))
  colnames(ggData) <- c(colnames(object@dr$tsne@cell.embeddings),gsub("-",".",features.plot))
  # print(feature.tmp)
  # ggData[,feature.tmp] <- t(object@data[feature.tmp,])
  ggl <- lapply(features.plot, function(feature){
    ggplot(ggData) + geom_point(mapping = aes_string(x = "tSNE_1", y = "tSNE_2", color = gsub("-",".",feature)), size = 3) +
      scale_color_gradientn(colours = c("grey","yellow","red")) +
      theme(legend.title = element_blank(),axis.title = element_blank()) + ggtitle(feature)
  })
  grid.arrange(grobs = ggl, nrow= nrow,ncol = ncol)
}


textCapitalize <- function(x){
  if(!is.character(x)){
    return(x)
  }
  tmp <- strsplit(as.character(x), split = "")[[1]]
  paste0(toupper(tmp[1]), paste0(tmp[-1], collapse = ""))
}
textInput2PCs <- function(x){
  tmp <- textInput2genes(x)
  x <- NULL
  for(i in tmp){
    if(grepl("-",i)){
      range = as.numeric(strsplit(i, "-")[[1]])
      x <- c(x, seq(range[1], range[2]))
    }else{
      x <- c(x, as.numeric(i))
    }
  }
}
textInput2genes <- function(x){
  x <- gsub("'", "", x)
  x <- gsub('"', "", x)
  x <- unique(unlist(strsplit(x, split = ' |\n|,|\t')))
  sapply(x, textCapitalize)
}

suffixFileName <- function(x, extension){
  paste0(x, format(Sys.time(),'_%Y%m%d_%H%M%S'), extension)
}

return.numeric <- function(x, ret = NULL){
  ifelse(is.numeric(x), x, ret)
}

checkdata.showModal <- function(x, label = NULL, message = NULL, show = T){
  #print(dim(x))
  text.null <- paste0("Please Load ", label)
  text.invalid <- paste0("Please check Your Data! It seems like 0 row/column in the ", label)

  if(!is.null(message)){
    text.null <- text.invalid <- message
  }

  if(is.null(x)){
    if(show) showModal(modalDialog(text.null,title = "Check Data", easyClose = T, fade = T))
    return(F)
  }else if(min(dim(x)) == 0){
    if(show) showModal(modalDialog(text.invalid, title = "Check Data", easyClose = T, fade = T))
    return(F)
  }else{
    return(T)
  }
}

## for guessCategory
guessNumeric <- function(x, exclude = c(NA, "NA", NULL, "NULL", ""," ","null")){
  return(!anyNA(as.numeric(setdiff(x, exclude))))
}
fillNA <- function(x, asNA = c(NA, "NA", NULL, "NULL", ""," ","null")){
  x[x %in% asNA] <- NA
  return(x)
}
formatAnnot <- function(annot, na.fill = T){
  numericCol <- apply(annot, 2, guessNumeric)
  if(na.fill){
    annot <- fillNA(annot, asNA = c(""," "))
  }
  if(any(numericCol)){
    annot[,numericCol] <- apply(annot[,numericCol], 2, as.numeric)
  }
  return(annot)
}
getAnnotUI <- function(annot){
  numericCol <- apply(annot, 2, guessNumeric)
  categoryAnnot <- annot[,!numericCol, drop = F]
  annotUI <- lapply(colnames(categoryAnnot), function(x){
                    cats <- unique(sort(categoryAnnot[,x]))
                    checkboxGroupInput(inputId = paste0("annot.cats.",gsub(pattern = " |-",replacement = ".",x)), label = x,
                                       choices = c("All",cats), selected = c("All",cats), inline = T, width = "100%")
                  })
  #print(annotUI)
  if(any(numericCol)){
    numericAnnot <- annot[,numericCol, drop = F]
    if(ncol(numericAnnot) %% 2 != 0){
      annotUI <- c(annotUI,
                   p("Error: Annotation data includes an ODD number of coordinate data!")
      )
    }else{
      annotUI <- c(annotUI,
                   list(selectInput(inputId = "annot.numeric", label = "Coordinate Annotation",
                               choices = colnames(numericAnnot), selected = colnames(numericAnnot), multiple = T, width = "100%"))
      )
    }
  }
  return(annotUI)
}
