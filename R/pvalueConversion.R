## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .pvalueCoversion
##' @description data conversion
##' In standard BED file, significant value of peak signal defined as Score column,
##' so we need to convert it as p.value (- log(p.value), -10 log(p.value), -100 log(p.value))
##'
##' @param x GRanges objects set of peak Interval
##' @param pvalueBase parameter for user defined pvalue format ( - log(p.value), -10 log(p.value), -100 log(p.value))
##' @return GRanges
##' @export
##' @importFrom rtracklayer score
##' @importFrom rtracklayer mcols
##' @author Julaiti Shayiding
##' @usage
##' .pvalueConversion(gr, 1L)

.pvalueConversion <- function(x, pvalueBase = 1L, ...) {
  # input param checking
  stopifnot(class(x) == "GRanges")
  stopifnot(is.numeric(pvalueBase))
  # explore score of all features
  if(is.null(x$pvalue)){
    x$p.value <- 10^(score(x)/(- pvalueBase))
    colnames(mcols(x))[3] <- "p.value"
  } else {
    x
  }
  res <- x
  return(res)
}
