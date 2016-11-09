# Bioconductor Package for Multiple Sample Peak Calling
#
#' @title readPeakFiles
#' @param  peakFolder set bed files to be read as GRanges objects
#' @return GRanges object
#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext
#' @description
#' reading bed format peak files as GRanges object to ease genomic interval manipulation
#' @author  Julaiti Shayiding
#' @example
## myInput <- readPeakFiles(peakFolder = "data/", verbose = FALSE)

readPeakFiles <- function(peakFolder, verbose=FALSE, ...) {
  # input param checking
  if(missing(peakFolder)) {
    stop("input param is missing!")
  }
  if (verbose)
    cat(">> reading all peakfiles from given folder...\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  stopifnot(length(peakFolder)>=1)
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- setNames(
    lapply(files, function(ele_) {
      out <- as(import.bed(ele_), "GRanges")
    }), tools::file_path_sans_ext(basename(files))
  )
  res <- f.read
  return(res)
}
