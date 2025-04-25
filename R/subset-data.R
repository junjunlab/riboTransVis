
#' Subset Data in `serp` and `ribotrans` Objects
#'
#' @description
#' A generic method to subset `serp` and `ribotrans` objects. This function filters internal data
#' slots (`bam_file`, `library`, `summary_info`) based on user-defined criteria.
#'
#' @param object An object of class `serp` or `ribotrans`.
#' @param ... Logical conditions used for subsetting data frames inside the object.
#' These are passed to the standard `subset()` function.
#'
#' @return An updated object of class `serp` or `ribotrans`, with internal slots filtered
#' according to the given criteria.
#'
#' @details
#' For `ribotrans` objects, this method:
#' \itemize{
#'   \item Subsets the `bam_file` data frame based on `...`.
#'   \item Subsets the `library` data frame based on `...`.
#'   \item Filters `summary_info` to include only rows where `sample` exists in the filtered `library`.
#' }
#'
#' For `serp` objects, this method:
#' \itemize{
#'   \item Subsets the `bam_file` data frame based on `...`.
#'   \item Subsets the `library` data frame based on `...`.
#'   \item Filters `summary_info` where the `sample` column matches
#'   `paste(lib$sample, lib$type, sep = "-")`, typically used for compound sample IDs.
#' }
#'
#' @examples
#' \dontrun{
#' # Subset a ribotrans object to include only a specific sample
#' subset_data(ribo, sample == "SampleA")
#'
#' # Subset a serp object based on type in library
#' subset_data(serp_obj, type == "ip")
#' }
#'
#' @export
setGeneric("subset_data",function(object,...) standardGeneric("subset_data"))



#' @describeIn subset_data Method for ribotrans objects
#' @export
setMethod("subset_data",
          signature(object = "ribotrans"),
          function(object,...){
            bm <- subset(object@bam_file, ...)
            object@bam_file <- bm

            lib <- subset(object@library, ...)
            object@library <- lib

            sry <- subset(object@summary_info, sample %in% lib$sample)
            object@summary_info <- sry

            return(object)
          }
)



#' @describeIn subset_data Method for serp objects
#' @export
setMethod("subset_data",
          signature(object = "serp"),
          function(object,...){
            bm <- subset(object@bam_file, ...)
            object@bam_file <- bm

            lib <- subset(object@library, ...)
            object@library <- lib

            sry <- subset(object@summary_info, sample %in% paste(lib$sample, lib$type,sep = "-"))
            object@summary_info <- sry

            return(object)
          }
)
