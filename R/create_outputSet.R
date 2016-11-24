## MSPC Project -- Bioconductor Package for Multiple Sample Peak Calling
##
##' @title create_outputSet
##' @description
##' more detailed description is needed
##'
##' @param peaklist_1 list of confirmed peaks
##' @param peaklist_2 list of discarded peaks
##' @param tau.s threshold value for stringent peaks
##' @param outDir user can control where exported BED file goes
##' @return BED file
##' @export
##' @importFrom rtracklayer export.bed
##' @import magrittr
##' @importFrom tibble rownames_to_column
##' @import magrittr
##' @author Julaiti Shayiding

library(tidyverse)
library(magrittr)
options(scipen = 0)

create_output <- function(output_path, list1, list2, tau.s) {
  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path))
  }

  list1 <- lapply(list1, as.data.frame)
  list2 <- lapply(list2, as.data.frame)

  names(list1) <- paste("confirmed", names(list1), sep = ".")
  names(list2) <- paste("discarded", names(list2), sep = ".")
  both <- do.call(rbind, c(list1, list2))
  both %<>% rownames_to_column(var = "cn")
  both %<>% separate(cn, c("original_list", "letters", "seq"), sep = "\\.")
  both %<>% mutate(peakStringency = ifelse(p.value <= tau.s , "Stringent", "Weak"))

  list_of_dfs <- both %>% split(list(.$letters, .$peakStringency, .$original_list))
  csv_names <- paste0(output_path, names(list_of_dfs), ".bed")
  return(mapply(export.bed, list_of_dfs, csv_names))
}

#' @example
options(scipen = 0)
create_output(list1 = confirmedERs, list2 = discardedERs, tau.s = 1.0E-08, output_path ="test/")

