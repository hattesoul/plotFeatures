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

# path to output file
output_dir = "./output"
#output_dir = "/media/user/data1/scRNAseq/rara/analysis/reanalyze/rara_p2_c_re/outs "

# ID for output file name
plot_file_prefix = "rara_unnorm_aggr"

# limits for UMI count
minUMICount <- 0
maxUMICount <- 0 # if maxUMICount = 0 then no limit is assumed

# limits for unique feature count
minUniqueFeatureCount <- 0
maxUniqueFeatureCount <- 0 # if maxUniqueFeatureCount = 0 then no limit is assumed

# append "/" to path if necessary
matrix_dir <- gsub("^\\s+|\\s+$", "", matrix_dir)
if (substr(matrix_dir, nchar(matrix_dir), nchar(matrix_dir)) != "/") {
  matrix_dir <- glue("{matrix_dir}/")
}

# append "/" to path if necessary
output_dir <- gsub("^\\s+|\\s+$", "", output_dir)
if (substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") {
  output_dir <- glue("{output_dir}/")
}

# check if files are accessible
cat("checking matrix files ...", sep = "")
if (!!file.access(paste0(matrix_dir, "barcodes.tsv.gz"), mode = 4)) {
  fileError(paste0(matrix_dir, "barcodes.tsv.gz"))
}
if (!!file.access(paste0(matrix_dir, "features.tsv.gz"), mode = 4)) {
  fileError(paste0(matrix_dir, "features.tsv.gz"))
}
if (!!file.access(paste0(matrix_dir, "matrix.mtx.gz"), mode = 4)) {
  fileError(paste0(matrix_dir, "matrix.mtx.gz"))
}
cat(" done.\n")

cat("checking write permissions ...", sep = "")
if (!!file.access(".", mode = 2)) {
  fileError(paste0(output_dir, plot_file_prefix))
}
cat(" done.\n")

# check if limits are set properly, swap them if not
if (maxUMICount > 0) {
  if (maxUMICount < minUMICount) {
    cat("WARNING! maxUMICount (", maxUMICount, ") is not greater than minUMICount (", minUMICount, "). Swapping UMI limits ...", sep = "")
    minUMICount <- minUMICount + maxUMICount
    maxUMICount <- minUMICount - maxUMICount
    minUMICount <- minUMICount - maxUMICount
    cat(" done.\n")
  }
}
if (maxUniqueFeatureCount > 0) {
  if (maxUniqueFeatureCount < minUniqueFeatureCount) {
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

# print numerical values rather in fixed notation than in exponential notation
options("scipen" = 10)

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
cat("summary:\n  barcodes: ", attr(mat,'Dim')[2], "\n  features: ", attr(mat,'Dim')[1], "\n", sep = "")

# remove non-expressed features
cat("removing features that are not expressed at all ...")
reducedMat = mat[apply(mat, 1, function(row) any(row !=0 )), ]
cat(" done.\n")
cat("summary:\n  removed features: ", dim(mat)[1] - dim(reducedMat)[1], "\n", sep = "")
cat("summary:\n  features left: ", dim(reducedMat)[1], "\n", sep = "")
#dim(reducedMat)[1]

# count UMIs
cat("counting UMIs per barcode ...")
UMISums = colSums(reducedMat)
cat(" done.\n")

# apply UMI filter to barcodes (cells)
if (maxUMICount > 0) {
  cat("filtering barcodes with minimum ", minUMICount, " and maximum ", maxUMICount, " UMIs ...", sep = "")
} else {
  cat("filtering barcodes with minimum ", minUMICount, " UMIs ...", sep = "")
} 
preValidBarcodes <- vector()
preValidIndices <- vector()
for (i in 1:length(UMISums)) {
  if (UMISums[i] >= minUMICount) {
    if (maxUMICount > 0) {
      if (UMISums[i] <= maxUMICount) {
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
if (maxUMICount > 0) {
  cat("summary:\n  barcodes with minimum ", minUMICount, " and maximum ", maxUMICount, " UMIs: ", length(preValidBarcodes), "\n", sep = "")
} else {
  cat("summary:\n  barcodes with minimum ", minUMICount, " UMIs: ", length(preValidBarcodes), "\n", sep = "")
}

# count unique features
cat("counting unique features per barcode ...")
uniqueFeatureSums = colSums(reducedMat[ ,preValidIndices] != 0)
cat(" done.\n")

# apply unique feature filter to barcodes (cells)
if (maxUniqueFeatureCount > 0) {
  cat("filtering barcodes with minimum ", minUniqueFeatureCount, " and maximum ", maxUniqueFeatureCount, " unique features ...", sep = "")
} else {
  cat("filtering barcodes with minimum ", minUniqueFeatureCount, " unique features ...", sep = "")
} 
validBarcodes <- vector()
for (i in 1:length(uniqueFeatureSums)) {
  if (uniqueFeatureSums[i] >= minUniqueFeatureCount) {
    if (maxUniqueFeatureCount > 0) {
      if (uniqueFeatureSums[i] <= maxUniqueFeatureCount) {
        validBarcodes <- c(validBarcodes, uniqueFeatureSums[i])
      }
    } else {
      validBarcodes <- c(validBarcodes, uniqueFeatureSums[i])
    }
  }
}
cat(" done.\n")
if (maxUniqueFeatureCount > 0) {
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
  codes <- c(codes, match(attr(validBarcodes, "names")[i], dimnames(reducedMat)[[2]]))
}
filtered_mat <- reducedMat[, codes]
cat(" done.\n")

# count UMIs of filtered barcodes
cat("counting UMIs per barcode ...")
validUMISums = colSums(filtered_mat)
cat(" done.\n")

# count unique features of filtered barcodes
cat("counting unique features per barcode ...")
validUniqueFeatureSums = colSums(filtered_mat != 0)
cat(" done.\n")

# manually calculate bin sizes and bin numbers for relative values
cat("ploting unique features/UMIs per barcode ...")

# set bin sizes
binSizeUniqueFeatures = 250
binSizeUMIs = 500

# calculate upper limits
ceiling_dec <- function(x, level = 1) round(x + 5 * 10 ^ (-level - 1), level)
upperLimitUniqueFeatures = ceiling_dec(max(validUniqueFeatureSums), -3)
upperLimitUMIs = ceiling_dec(max(validUMISums), -4)

# calculate number of bins
binNumbersUniqueFeatures = upperLimitUniqueFeatures/binSizeUniqueFeatures
binNumbersUMIs = upperLimitUMIs/binSizeUMIs

# calculate breakpoints for bins
breakPointsUniqueFeatures = c(0, c(1:binNumbersUniqueFeatures) * binSizeUniqueFeatures - 1)
breakPointsUMIs = c(1:binNumbersUMIs) * binSizeUMIs - 1

# calculate values for bins
tempHistUniqueFeatures = hist(validUniqueFeatureSums, breaks = breakPointsUniqueFeatures, plot = FALSE)
tempHistUMIs = hist(validUMISums, breaks = breakPointsUMIs, plot = FALSE)

# calculate relative values for bins
relativeUniqueFeatureBins = round(tempHistUniqueFeatures$counts/sum(tempHistUniqueFeatures$counts) * 100, digits = 3)
relativeUMIBins = round(tempHistUMIs$counts/sum(tempHistUMIs$counts) * 100, digits = 3)

# Plotly seems to skip leading zero in histograms
# workaround: correct bin values by removing leading zeros
while(relativeUniqueFeatureBins[1] == 0) {
  relativeUniqueFeatureBins <- relativeUniqueFeatureBins[2:length(relativeUniqueFeatureBins)]
}
while(relativeUMIBins[1] == 0) {
  relativeUMIBins <- relativeUMIBins[2:length(relativeUMIBins)]
}

# create text objects for plot tooltips
relativeUniqueFeatureText = c(relativeUniqueFeatureBins, rep(NA_character_, length(validUniqueFeatureSums) - length(relativeUniqueFeatureBins)))
relativeUMIText = c(relativeUMIBins, rep(NA_character_, length(validUMISums) - length(relativeUMIBins)))

xlabel <- list(title = "number of unique features/UMIs", rangemode = "tozero")
ylabel <- list(title = "number of cells")
p <- plot_ly(
  alpha = 0.6
) %>%
add_histogram(x = ~(validUniqueFeatureSums), xbins = list(start = 0, size = binSizeUniqueFeatures, end = upperLimitUniqueFeatures), name = glue("unique features
  (bin size: {binSizeUniqueFeatures})"), text = relativeUniqueFeatureText, hovertemplate = paste('<b>%{y} (%{text} %)</b> cells have <b>%{x}</b> unique features', sep = "")) %>%
add_histogram(x = ~(validUMISums), xbins = list(start = 0, size = binSizeUMIs, end = upperLimitUMIs), name = glue("UMIs
  (bin size: {binSizeUMIs})"), text = relativeUMIText, hovertemplate = paste('<b>%{y} (%{text} %)</b> cells have <b>%{x}</b> UMIs', sep = "")) %>%
layout(barmode = "overlay", title = glue("distribution of unique feature/UMI counts"), xaxis = xlabel, yaxis = ylabel, hovermode = "compare")
p
savefile <- glue("{output_dir}{plot_file_prefix}_feat_dist.html")
saveWidgetFix(as_widget(p), savefile)
cat(" done.\n")
