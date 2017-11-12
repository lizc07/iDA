# iDA
R package - iteractive Data Analysis, implemented under shiny
## install R package from Github  
```R
install.packages('devtools')
devtools::install_github('lizc07/iDA')
```
## Dependency
```R
source('http://bioconductor.org/biocLite.R')
biocLite(c('Seurat', 'shiny', 'shinydashboard', 'DT', 'igraph', 
'Rtsne', 'dplyr', 'ggplot2', 'gplots', 'gridExtra', 'pheatmap'))
# if necessary
# chooseCRANmirror()
# chooseBioCmirror()
```
