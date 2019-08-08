# Sparse and Dense Matrix Classes and Methods
library(Matrix)

# Interpreted String Literals
library(glue)

# Create Interactive Web Graphics via 'plotly.js'
library(plotly)

# HTML Widgets for R
library(htmlwidgets)

# wrapper for saveWidget
saveWidgetFix <- function (widget, file, ...) {
  ## A wrapper to saveWidget which compensates for arguable BUG in
  ## saveWidget which requires `file` to be in current working
  ## directory.
  wd <- getwd()
  on.exit(setwd(wd))
  outDir <- dirname(file)
  file <- basename(file)
  setwd(outDir);
  saveWidget(widget, file = file, ...)
}

# exit without additional error message
exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}

# file error message and exit
fileError <- function(a) {
  cat("\nFile error: ", a, " not found or no read/write permission.\nExiting.", sep = "")
  exit()
}

# init
# path to input files
matrix_dir = "/media/hattesoul/data1/scRNAseq/rara/analysis/aggr/v3-1-0_rara_unnorm_aggr/outs/filtered_feature_bc_matrix "
#matrix_dir = "/media/user/data1/scRNAseq/UK-1473/analysis/single_counts/UK-1473_6000/outs/filtered_feature_bc_matrix"

# path to features of interest file
# each line contains one feature, like TP63 or FOXJ1 etc.
foi_file = "./assets/foi.txt"

# path to output file
output_dir = "./output"
#output_dir = "/media/user/data1/scRNAseq/rara/analysis/reanalyze/rara_p2_c_re/outs "

# ID for output file name
plot_file_prefix = "rara_unnorm_aggr"

# limits for UMI count
minUMICount <- 10000
maxUMICount <- 100000 # if maxUMICount = 0 then no limit is assumed

# limits for unique feature count
minUniqueFeatureCount <- 2500
maxUniqueFeatureCount <- 100000 # if maxUniqueFeatureCount = 0 then no limit is assumed

# append "/" to path if necessary
matrix_dir <- gsub("^\\s+|\\s+$", "", matrix_dir)
if(substr(matrix_dir, nchar(matrix_dir), nchar(matrix_dir)) != "/"){
  matrix_dir <- glue("{matrix_dir}/")
}

# append "/" to path if necessary
output_dir <- gsub("^\\s+|\\s+$", "", output_dir)
if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/"){
  output_dir <- glue("{output_dir}/")
}

# check if files are accessible
cat("checking matrix files ...", sep = "")
if(!!file.access(paste0(matrix_dir, "barcodes.tsv.gz"), mode = 4)){
  fileError(paste0(matrix_dir, "barcodes.tsv.gz"))
}
if(!!file.access(paste0(matrix_dir, "features.tsv.gz"), mode = 4)){
  fileError(paste0(matrix_dir, "features.tsv.gz"))
}
if(!!file.access(paste0(matrix_dir, "matrix.mtx.gz"), mode = 4)){
  fileError(paste0(matrix_dir, "matrix.mtx.gz"))
}
cat(" done.\n")

cat("checking write permissions ...", sep = "")
if(!!file.access(".", mode = 2)){
  fileError(paste0(output_dir, plot_file_prefix))
}
cat(" done.\n")

# check if limits are set properly, swap them if not
if(maxUMICount > 0){
  if(maxUMICount < minUMICount){
    cat("WARNING! maxUMICount (", maxUMICount, ") is not greater than minUMICount (", minUMICount, "). Swapping UMI limits ...", sep = "")
    minUMICount <- minUMICount + maxUMICount
    maxUMICount <- minUMICount - maxUMICount
    minUMICount <- minUMICount - maxUMICount
    cat(" done.\n")
  }
}
if(maxUniqueFeatureCount > 0){
  if(maxUniqueFeatureCount < minUniqueFeatureCount){
    cat("WARNING! maxUniqueFeatureCount (", maxUniqueFeatureCount, ") is not greater than minUniqueFeatureCount (", minUniqueFeatureCount, "). Swapping unique feature limits ...", sep = "")
    minUniqueFeatureCount <- minUniqueFeatureCount + maxUniqueFeatureCount
    maxUniqueFeatureCount <- minUniqueFeatureCount - maxUniqueFeatureCount
    minUniqueFeatureCount <- minUniqueFeatureCount - maxUniqueFeatureCount
    cat(" done.\n")
  }
}

# load necessary Cell Ranger files
cat("loading barcodes from ", matrix_dir, " ...", sep = "")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
cat(" done.\n")

cat("loading features from ", matrix_dir, " ...", sep = "")
features.path <- paste0(matrix_dir, "features.tsv.gz")
cat(" done.\n")

cat("loading matrix from ", matrix_dir, " ...", sep = "")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
cat(" done.\n")

# add column and row names
cat("reading feature and barcode names ...")
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2
cat(" done.\n")

# load file with features of interest
cat("loading features of interest from ", foi_file, " ...", sep = "")
foi <- read.table(foi_file, header = FALSE)
cat(" done.\n")
cat("summary:\n  barcodes: ", attr(mat,'Dim')[2], "\n  features: ", attr(mat,'Dim')[1], "\n  features of interest: ", lengths(foi), "\n", sep = "")

# print numerical values rather in fixed notation than in exponential notation
options("scipen" = 10)

# count UMIs
cat("counting UMIs per barcode ...")
UMISums = colSums(mat)
cat(" done.\n")

# apply UMI filter to barcodes (cells)
if(maxUMICount > 0){
  cat("filtering barcodes with minimum ", minUMICount, " and maximum ", maxUMICount, " UMIs ...", sep = "")
} else {
  cat("filtering barcodes with minimum ", minUMICount, " UMIs ...", sep = "")
} 
preValidBarcodes <- vector()
preValidIndices <- vector()
for (i in 1:length(UMISums)) {
  if(UMISums[i] >= minUMICount){
    if(maxUMICount > 0){
      if(UMISums[i] <= maxUMICount){
        preValidBarcodes <- c(preValidBarcodes, UMISums[i])
        preValidIndices <- c(preValidIndices, i)
      }
    } else {
      preValidBarcodes <- c(preValidBarcodes, UMISums[i])
      preValidIndices <- c(preValidIndices, i)
    }
  }
}
cat(" done.\n")
if(maxUMICount > 0){
  cat("summary:\n  barcodes with minimum ", minUMICount, " and maximum ", maxUMICount, " UMIs: ", length(preValidBarcodes), "\n", sep = "")
} else {
  cat("summary:\n  barcodes with minimum ", minUMICount, " UMIs: ", length(preValidBarcodes), "\n", sep = "")
}

# count unique features
cat("counting unique features per barcode ...")
uniqueFeatureSums = colSums(mat[ ,preValidIndices] != 0)
cat(" done.\n")

# apply unique feature filter to barcodes (cells)
if(maxUniqueFeatureCount > 0){
  cat("filtering barcodes with minimum ", minUniqueFeatureCount, " and maximum ", maxUniqueFeatureCount, " unique features ...", sep = "")
} else {
  cat("filtering barcodes with minimum ", minUniqueFeatureCount, " unique features ...", sep = "")
} 
validBarcodes <- vector()
for (i in 1:length(uniqueFeatureSums)) {
  if(uniqueFeatureSums[i] >= minUniqueFeatureCount){
    if(maxUniqueFeatureCount > 0){
      if(uniqueFeatureSums[i] <= maxUniqueFeatureCount){
        validBarcodes <- c(validBarcodes, uniqueFeatureSums[i])
      }
    } else {
      validBarcodes <- c(validBarcodes, uniqueFeatureSums[i])
    }
  }
}
cat(" done.\n")
if(maxUniqueFeatureCount > 0){
  cat("summary:\n  barcodes with minimum ", minUniqueFeatureCount, " and maximum ", maxUniqueFeatureCount, " unique features: ", length(validBarcodes), "\n", sep = "")
} else {
  cat("summary:\n  barcodes with minimum ", minUniqueFeatureCount, " features: ", length(validBarcodes), "\n", sep = "")
}

# set default behavior for prining numerical values
options("scipen" = 0)

# reduce matrix to contain valid barcodes only
cat("filtering by valid barcodes ...")
filtered_mat <- matrix()
codes <- vector()
for (i in 1:length(validBarcodes)) {
  codes <- c(codes, match(attr(validBarcodes, "names")[i], dimnames(mat)[[2]]))
}
filtered_mat <- mat[, codes]
cat(" done.\n")

# get row numbers corresponding to features of interest
cat("filtering by features of interest ...")
hit <- vector()
for (i in 1:lengths(foi)) {
  hit <- c(hit, match(foi[[1]][i], dimnames(filtered_mat)[[1]]))
}
cat(" done.\n")

cat("ploting cell numbers with features of interest ...")
for (i in 1:length(hit)) {
  xlabel <- list(title = "number of gene copies")
  ylabel <- list(title = "number of cells")
  x <- filtered_mat[hit[i], ]
  t <- table(x)
  # get break points
  b <- c(min(strtoi((attr(t, "dimnames")$x))):max(strtoi((attr(t, "dimnames")$x))))
  bb <- c(0, b + 0.5)
  hcum <- h <- hist(x, breaks = bb, plot = FALSE)
  hcum$counts <- cumsum(hcum$counts)
  
  p <- plot_ly(alpha = 0.6) %>%
    add_bars(x = b, y = h$counts, type = "bar", name = "absolute numbers", hovertemplate = paste('<b>%{y}</b> cells have <b>%{x}</b> copies of <b>',foi[[1]][i],'</b>', sep = "")) %>%
    add_bars(x = b, y = hcum$counts, type = "bar", name = "cumulative numbers", hovertemplate = paste('<b>%{y}</b> cells have <b>%{x}</b> or less copies of <b>',foi[[1]][i],'</b>', sep = "")) %>%
    layout(barmode = "overlay", title = glue("absolute gene expression of {foi[[1]][i]}"), xaxis = xlabel, yaxis = ylabel)
  p
  savefile <- glue("{output_dir}{plot_file_prefix}_{foi[[1]][i]}.html")
  saveWidgetFix(as_widget(p), savefile)
}
cat(" done.\n")

#cat("saving plots ...")
#for (i in 1:length(list_foi)) {
#  png(glue("{plot_dir}/{plot_file_prefix}_{foi[[1]][i]}.png"), width = 1200, height = 1200)
#  plot(list_foi[[i]], main = glue("absolute gene expression of {foi[[1]][i]}"), ylab = "cell numbers", xlab = "gene copies")
#  dev.off()
#}
#cat(" done.\n")

####
#densityplot
## Make some sample data
#x <- sample(0:30, 200, replace=T, prob=15 - abs(15 - 0:30))
#x <- filtered_mat[hit[1], ]
#t <- table(x)
#b <- c(min(strtoi((attr(t, "dimnames")$x))):max(strtoi((attr(t, "dimnames")$x))))
#bb <- c(0, b + 0.5)
#b <- c(0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5)
#x <- t
## Calculate and plot the two histograms
#hcum <- h <- hist(x, breaks = c(1:length(x)-1), plot = FALSE)
#hcum <- h <- hist(x, breaks = c(min(strtoi((attr(t, "dimnames")$x))):max(strtoi((attr(t, "dimnames")$x)))), plot = FALSE)

#hcum <- h <- hist(x, breaks = bb, plot = FALSE)
#hcum$counts <- cumsum(hcum$counts)
#plot(hcum, main = "", freq = TRUE)
#plot(h, add = TRUE, col = "grey", freq = TRUE)

## Plot the density and cumulative density
#d <- density(x)
#d <- density(t)
#lines(x = d$x, y = d$y * length(bb) * diff(h$breaks)[1], lwd = 2)
#lines(x = d$x, y = cumsum(d$y)/max(cumsum(d$y)) * length(bb), lwd = 2)
