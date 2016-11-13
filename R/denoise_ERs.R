## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @description implementation for denoise peak interval
##' @title .denoise.ERs
##' @param grs GRangesList with peakInterval
##' @param tau.w threshold value for weakly enriched peak
##' @param .fileName user can name noise peak by custom which exported as BED file
##' @param outDir the folder where noisePeak bed goes
##' @param verbose control whether the output is printed or not
##' @return GRangesList with validated peakInterval
##' @export
##' @importFrom rtracklayer export.bed
##' @author Julaiti Shayiding


.denoise.ERs <- function(grs, tau.w= 1.0E-04, .fileName="", outDir=getwd(), verbose=FALSE, ...) {
  # check input param
  stopifnot(class(grs[[1L]])=="GRanges")
  stopifnot(length(grs)>0)
  stopifnot(is.numeric(tau.w))
  if (verbose) {
    cat(">> filter out all background noise peaks whose pvalue above threshold \t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  if(!dir.exists(outDir)) {
    dir.create(file.path(outDir))
    setwd(file.path(outDir))
  }
  res <- lapply(seq_along(grs), function(x) {
    .gr <- grs[[x]]
    .grNM <- names(grs)[x]
    .drop <- .gr[.gr$p.value > tau.w]
    export.bed(.drop, sprintf("%s/%s.%s.bed", outDir, .fileName, .grNM))
    .keep <- .gr[.gr$p.value <= tau.w]
    return(.keep)
  })
  rslt <- setNames(res, names(grs))
  return(rslt)
}

#' @example
total.ERs <- .denoise.ERs(myData, tau.w = 1.0E-04, .fileName = "noisePeak", outDir = "test/")
