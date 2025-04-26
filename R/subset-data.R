
#' Subset Data in `serp` and `ribotrans` Objects
#'
#' @description
#' A generic method to subset data contained within `serp` and `ribotrans` objects. This function filters internal data slots,
#' including `bam_file`, `library`, `summary_info`, and if present, `reads_offset_info`, based on user-specified logical conditions.
#'
#' @param object An object of class `serp` or `ribotrans`.
#' @param ... Logical expressions passed to the internal `subset()` functions to filter data frames inside the object.
#'
#' @return An updated `serp` or `ribotrans` object with internal data slots subsetted according to the specified criteria.
#'
#' @details
#' For `ribotrans` objects, the method:
#' \itemize{
#'   \item Subsets the `bam_file` data frame according to the filtering conditions.
#'   \item Subsets the `library` data frame according to the filtering conditions.
#'   \item Filters `summary_info` to retain only rows where the `sample` column is present in the filtered `library`.
#'   \item If the slot `reads_offset_info` contains rows, it is subsetted similarly by matching `sample` values in the filtered `library`.
#' }
#'
#' For `serp` objects, the method:
#' \itemize{
#'   \item Subsets the `bam_file` data frame according to the filtering conditions.
#'   \item Subsets the `library` data frame accordingly.
#'   \item Filters `summary_info` where the `sample` column matches the compound sample IDs formed by concatenating `lib$sample` and `lib$type` with a hyphen (`paste(lib$sample, lib$type, sep = "-")`).
#'   \item If `reads_offset_info` is present and non-empty, it is subsetted by samples in the filtered `library`.
#' }
#'
#' @examples
#' \dontrun{
#' # Subset a ribotrans object to include only samples named "SampleA"
#' subset_data(ribo, sample == "SampleA")
#'
#' # Subset a serp object to keep only 'ip' type samples
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

            if(nrow(object@reads_offset_info) > 0){
              offset <- subset(object@reads_offset_info, sample %in% lib$sample)
              object@reads_offset_info <- offset
            }

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

            if(nrow(object@reads_offset_info) > 0){
              offset <- subset(object@reads_offset_info, sample %in% paste(lib$sample, lib$type,sep = "-"))
              object@reads_offset_info <- offset
            }

            return(object)
          }
)
