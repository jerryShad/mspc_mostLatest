# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title readPeakFiles
#' @param  peakFolder set of bed files to be read as GRanges objects
#' @param pvalueBase user can choose the scale of pvalue by custom
#' @param verbose logical that control whether the output is printed or not
#' @return GRanges object
#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext
#' @description
#' reading bed format peak files as GRanges object to ease genomic interval manipulation
#' @author  Julaiti Shayiding

readPeakFiles <- function(peakFolder, pvalueBase = 1L, verbose=FALSE, ...) {
  # input param checking
  if (verbose)
    cat(">> reading all peakfiles from given folder...\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  stopifnot(length(peakFolder)>0)
  stopifnot(is.numeric(pvalueBase))
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- setNames(
    lapply(files, function(ele_) {
      .gr <- as(import.bed(ele_), "GRanges")
      if(is.null(.gr$p.value)) {
        .gr <- .pvalueConversion(.gr, 1L)
      }
    }), tools::file_path_sans_ext(basename(files))
  )
  res <- f.read
  return(res)
}

#' @example
#' inBED <- readPeakFiles(peakFolder = "test/testData/", pvalueBase=1L,verbose = FALSE)
