## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##' @description
##' all updated function code and workflow of MSPC Packages
##'
##==============================================================================================
##' @details read bed files as GRanges objects from file directory
##' function implementation

readPeakFiles <- function(peakFolder, verbose=FALSE, ...) {
  # input param checking
  require
  if(missing(peakFolder)) {
    stop("input param is missing!")
  }
  stopifnot(length(peakFolder)>=1)
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- setNames(
    lapply(files, function(ele_) {
      res <- as(import.bed(ele_), "GRanges")

    }), tools::file_path_sans_ext(basename(files))
  )
  res <- f.read
  return(res)
}

##' @examples
myData <- readPeakFiles(peakFolder = "test/testData/")

##==============================================================================================
##' @details splitting GRanges
##' function implementation
##' FIXME : S4 subscript error

.denoise.ERs <- function(grs, tau.w=1.0E-04, .fileName="", outDir=getwd(), verbose=FALSE, ...) {
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
  res <- lapply(grs, function(x) {
    if(is.null(x$p.value)) {
      x <- .pvalueConversion(x, pvalueBase = 1L)
      .gr <- grs[[x]]
      .grNM <- names(grs)[x]
      .noiseERs <- .gr[.gr$p.value > tau.w]
      export.bed(.noiseERs, sprintf("%s/%s.%s.bed", outDir, .fileName, .grNM), row.names = FALSE)
      .ERs <- .gr[.gr$p.value <= tau.w]
      return(.ERs)
    }
  })
  rslt <- setNames(res, names(grs))
  return(rslt)
}

#' @example
total.ERs <- .denoise.ERs(myData, tau.w = 1.0E-08, .fileName = "noisePeak", outDir = "test/")

##========================================================

.pvalueConversion <- function(x, pvalueBase = 1L, ...) {
  stopifnot(class(x) == "GRanges")
  stopifnot(is.numeric(pvalueBase))
  # explore score of all features
  if(is.null(x$p.value)){
    x$p.value <- 10^(score(x)/(- pvalueBase))
    colnames(mcols(x))[3] <- "p.value"
  } else {
    x
  }
  return(x)
}

##==============================================================================================
##' @details peakOverlapping
##' function implementation
## I am on the right track, keep coding *_*

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

#' @example
#' Hit <- .peakOverlapping(peakset = total.ERs, FUN = which.max)

##=============================================================================================

MSPC.Analyzer <- function(peakset, ovHit, replicate.type=c("Biological","Technical"), tau.s=1.0E-08, ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  stopifnot(is.numeric(tau.s))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  length(peakset)-1,
                  length(peakset))
  cnt.ovHit <- as.matrix(Reduce('+', lapply(ovHit, lengths)))
  keepMe <- cnt.ovHit >= min.c

  dropList <- lapply(ovHit, function(ele_) ele_[!keep_me])
  init.discardPeaks <- Map(unlist,
                           mapply(extractList, peakset, dropList))

  keepList <- lapply(ovHit, function(x) x[keepMe])

  pval_List <- mapply(.get.pvalue, keepList, peakset)
  .helper.PV <- function(p.list) {
    res <- sapply(p.list, function(x) {
      out <- ifelse(length(x)>0,
                    x,
                    0)
      out
    })
  }
  pval.TB <- as.data.frame(Map(.helper.PV, pval_List))
  #pv.stat <- rowSums( -2*log( pval.TB ), na.rm=T )
  #Npresent <- rowSums( !is.na(pval.TB) )
  #comb.pval <- pchisq( pv.stat, df=2*Npresent, lower=FALSE )
  comb.pval <- suppressWarnings(
    .globSC <- apply(pval.TB, 1, function(ro) {
      require(metap)
      res <- sumlog(ro)$p
    })
  )
  comb.pval <- as.matrix(comb.pval)
  ## FIXME: to make more compatible
  Confirmed_idx <- lapply(keepList, function(elm) {
    saved <- sapply(comb.pval, function(x) x <= tau.s)
    res <- elm[saved]
  })

  # warning : do not remove duplicates, try to find pattern
  confirmed <- Map(extractList, peakset, Confirmed_idx)
  .Confirmed.ERs <- lapply(confirmed, function(ele_) unlist(ele_))
  ## TODO
  Discarded_idx <- lapply(keepList, function(elm) {
    droped <- sapply(comb.pval, function(x) x > tau.s)
    res <- elm[droped]
  })

  .Fisher.discPeaks <- Map(extractList, peakset, Discarded_idx)
  .Fisher.discPeaks <- lapply(.Fisher.discPeaks, function(elm) unlist(elm))
  .Discarded.ERs <- suppressWarnings(mapply(c, .init.discPeaks, .Fisher.discPeaks))

  ##-------------------------------------------------------------------------------
  ## FIXME : optimize me
  confirmed <- lapply(seq_along(all_Confirmed), function(x) {
    if(x==1) {
      unique(all_Confirmed[[x]])
    } else {
      all_Confirmed[[x]]
    }
  })

  discarded. <- lapply(seq_along(all_Discarded), function(x) {
    if(x==1)
      unique(all_Discarded[[x]])
    else
      all_Discarded[[x]]
  })

  ##-------------------------------------------------------------------------------
  .setPurification <- ifelse(replicate.type=="Biological",
                             res <- .Confirmed.ERs,
                             res <- Map(anti_join, .Confirmed.ERs, .Discarded.ERs))


  .create_OUTP <- function(peaks, pAdjustMethod="BH", alpha=0.05, ...) {
    # input param checking
    stopifnot(class(peaks)=="GRanges")
    stopifnot(is.numeric(alpha))
    pAdjustMethod = match.arg(pAdjustMethod)
    if(is.null(peaks$p.value)) {
      stop("required slot is missing")
    } else {
      p <- peaks$p.value
      p.adj <- p.adjust(p, method = pAdjustMethod)
      peaks$p.adj <- p.adj
      rslt <- split(peaks, ifelse(peaks$p.adj <= alpha,
                                  "pass", "fail"))
      return(rslt)
    }
  }
  #' @example
  .BH_output <- Map(.create_OUTP, .setPurification)

  ##---------------------------------------
  rslt <- c(.Confirmed.ERs, .Discarded.ERs)
  ans.Reslt <- split(rslt, sub("_.*", "", names(rslt)))[sub("_.*", "", names(.Confirmed.ERs))]
  return(ans.Reslt)
}

.get.pvalue <- function(hit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, hit)
  return(res)
}

##================================================================================================
##' @description
##' This script is used for how to get stringent peakset, weakpeaks set like Musera

library(data.table)
confirmed_dt <- rbindlist(all_Confirmed,idcol=TRUE)
discarded_dt <- rbindlist(all_Discarded,idcol=TRUE)

DT <- rbind(confirmed_dt, discarded_dt)
DT[,gr := ifelse(p.value <= tau.S,"Stringent","Weak"),.id]
DT[,.id:= gsub("[.](confirmed|discarded)","",.id)]
res <- by(DT,DT$.id,FUN = function(x) split(x,x$gr))


##================================================================================================
##' @description
##' This scripts how to manipulate list of confirmed peaks, discarded peaks for BH correction test
##'
x <- c(all_Confirmed, all_Discarded)
ans.Reslt <- split(x, sub("_.*", "", names(x)))[sub("_.*", "", names(all_Confirmed))]

##================================================================================================
.setPurification <- ifelse(replicate.type=="Biological",
                           res <- .Confirmed.ERs,
                           res <- Map(anti_join, .Confirmed.ERs, .Discarded.ERs))

#-------------------------------------------------------------------------------------------------
# dummy function

func <- function(L1, L2, replicate.type=c("Biological","Technical")) {
  replicate.type = match.arg(replicate.type)
  if(replicate.type=="Biological") {
    res <- L1
  } else {
    res <- Map(function(x,y) anti_join(x,y), L1, L2)
  }
  return(res)
}

L1 <- Map("data.frame", confirmed.)
L2 <- Map("data.frame", discarded.)

# testme:
.setPurif <- func(L1, L2, "Biological")
.setPurif <- lapply(.setPurif, function(x) as(x, "GRanges"))

#---------------------------------------------------------------------------------------------------

.create_OUTP <- function(peaks, pAdjustMethod="BH", alpha=0.05, ...) {
  # input param checking
  stopifnot(class(peaks)=="GRanges")
  stopifnot(is.numeric(alpha))
  pAdjustMethod = match.arg(pAdjustMethod)
  if(is.null(peaks$p.value)) {
    stop("required slot is missing")
  } else {
    p <- peaks$p.value
    p.adj <- p.adjust(p, method = pAdjustMethod)
    peaks$p.adj <- p.adj
    rslt <- split(peaks, ifelse(peaks$p.adj <= alpha,
                                "pass", "fail"))
    return(rslt)
  }
}

#' @example
.BH_output <- Map(.create_OUTP, .setPurif)
.BH_output <- lapply(.BH_output, unique)
##==================================================================================================


