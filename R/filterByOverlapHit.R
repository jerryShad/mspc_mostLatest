## MSPC Package -- Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .filterByOverlapHit
##' @param .ovHit list of overlap hit index
##' @param peakset set of peakInterval
##' @param replicate.type type of replicate user' input
##' @param isSuffOverlap logical vector that check whether sufficient overlap or not
##' @return
##' @export
##' @importFrom S4Vectors Reduce
##' @importFrom S4Vectors lengths
##' @importFrom XVector extractList
##' @details
##' @author Julaiti Shayiding

.filterByOverlapHit <- function(.ovHit, peakset, replicate.type=c("Biological", "Technical"),
                   isSuffOverlap=TRUE, verbose=FALSE, ...) {
  # check input param
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  stopifnot(missing(isSuffOverlap))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  length(peakset)-1,
                  length(peakset))
  cnt.ovHit <- as.matrix(Reduce('+', lapply(hit.List, lengths)))
  if(isSuffOverlap) {
    keepHit <- lapply(hit.List, function(ele_) {
      keepMe <- sapply(cnt.ovHit, function(x) x >= min.c)
      res <- ele_[keepMe]
    })
    return(keepHit)
  } else {
    dropHit <- lapply(hit.List, function(ele_) {
      droped <- sapply(cnt.ovHit, function(x) x < min.c)
      res <- ele_[droped]
    })
    rslt <- Map(unlist,
                mapply(extractList, peakset, dropHit))
    return(rslt)
  }
}

#' @example
keepList <- func.1(Hit, peakset = total.ERs, replicate.type = "Biological", isSuffOverlap=TRUE)
initDisc.ERs <- func.1(Hit, peakset = total.ERs, replicate.type = "Biological", isSuffOverlap=FALSE)
