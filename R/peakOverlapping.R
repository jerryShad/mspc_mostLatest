## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .peakOverlapping
##' @param peakset peak files
##' @param FUN function selected by user needs when keeping most stringent peak when multiple overlap occur
##' @return IntegerList overlap-hit IntegerList across multiple GRanges simulanously
##' @export
##' @importFrom GenomicRanges findOverlaps
##' @importFrom XVector extractList
##' @importfrom rtracklayer score
##' @importFrom IRanges which.min
##' @author Julaiti Shayiding
##' @example

.peakOverlapping <- function(peakset, FUN=which.max, ...) {
  # input param checking
  stopifnot(inherits(peakset[[1]], "GRanges"))
  res <- list()
  for(i in seq_along(peakset)) {
    que <- peakset[[i]]
    queHit <- as(findOverlaps(que), "List")
    supHit <- lapply(peakset[- i], function(ele_) {
      ans <- as(findOverlaps(que, ele_), "List")
      out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
      out.idx0 <- out.idx0[!is.na(out.idx0),]
      ans <- ans[out.idx0]
      ans
    })
    res[[i]] = DataFrame(c(list(que=queHit), sup=supHit))
    names(res[[i]]) = c(names(peakset[i]),names(peakset[- i]))
  }
  rslt <- lapply(res, function(x) as.matrix(x[names(res[[1]])]))
  rslt <- DataFrame(rbind(rslt[[1]],
                          unique(do.call("rbind", rslt[2: length(rslt)]))))
  rslt <- lapply(rslt, function(x) as(x, "CompressedIntegerList"))
  return(rslt)
}

