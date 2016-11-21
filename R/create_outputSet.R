## MSPC Project -- Bioconductor Package for Multiple Sample Peak Calling
##
##' @description
##' @title create_outputSet
##' @param peaklist
##' @param peaklist
##' @param tau.s
##' @param outDir
##' @return BED file
##' @export
##' @importFrom rtracklayer export.bed
##' @importFrom tibble rownames_to_column
##' @import magrittr
##' @author Julaiti Shayiding

library(tidyverse)
library(magrittr)

create.outputSet <- function(mList, mList, threshold, outDir=getwd(), ...) {
  # check input param
  if(!dir.exists(outDir)) {
    dir.create(file.path(outDir))
    setwd(file.path(outDir))
  }
  both <- do.call("rbind", c(.Confirmed.ERs, .Discarded.ERs))
  both %<>% rownames_to_column(var = "cn")
  both %<>% separate(cn, c("letters", "seq"), sep = "\\.")
  both %<>% mutate(isPassed = ifelse(.score <= threshold, "Stringent", "Weak"),
                   isDiscard = ifelse(is.na(seq), ".Confirmed", ".Discarded"))

  list_of_dfs <- both %>% split(list(.$letters, .$isPassed, .$isDiscard))
  file_names <- paste0("test/", names(list_of_dfs), ".bed") # change this path
  res <- mapply(export.bed, list_of_dfs, file_names)
  return(res)
}

#' @example
#' output <- create.outputSet(.confirmed.ERs, .discarded.ERs, alpha=0.05, outDir="")
