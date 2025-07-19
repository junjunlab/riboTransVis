#' Trim Adapters from FASTQ Files
#'
#' This function trims adapters from single-end or paired-end FASTQ files using the rfastp package.
#' It processes multiple FASTQ files in batch and saves the trimmed results to a specified output directory.
#'
#' @param read1 Character vector of paths to Read 1 FASTQ files. Required parameter.
#' @param read2 Character vector of paths to Read 2 FASTQ files for paired-end data.
#'   Default is empty string, which will be converted to a vector of empty strings matching the length of read1.
#' @param output_dir Character string specifying the output directory for trimmed files.
#'   Default is "2.trim-data". Directory will be created if it doesn't exist.
#' @param outputFastq Character vector specifying the output FASTQ file names.
#'   Should match the length of read1. Default is NULL.
#' @param adapterSequenceRead1 Character string specifying the adapter sequence for Read 1.
#'   Default is "auto" for automatic detection.
#' @param adapterSequenceRead2 Character string specifying the adapter sequence for Read 2.
#'   Default is "auto" for automatic detection.
#' @param minReadLength Integer specifying the minimum read length after trimming.
#'   Reads shorter than this will be discarded. Default is 15.
#' @param maxReadLength Integer specifying the maximum read length after trimming.
#'   Default is 0 (no maximum length filtering).
#' @param thread Integer specifying the number of threads to use for processing. Default is 1.
#' @param verbose Logical indicating whether to print verbose output during processing. Default is FALSE.
#' @param ... Additional parameters passed to the rfastp function.
#'
#' @return Invisible NULL. The function is called for its side effects (creating trimmed FASTQ files).
#'   Progress messages are printed to the console for each processed file pair.
#'
#' @details
#' The function creates the output directory if it doesn't exist and processes each FASTQ file
#' (or file pair for paired-end data) using the rfastp function. For single-end data,
#' leave read2 as the default empty string. The function will automatically handle the
#' conversion to process single-end files.
#'
#' @examples
#' \dontrun{
#' # Single-end data
#' trim_adaptors(
#'   read1 = c("sample1.fastq", "sample2.fastq"),
#'   outputFastq = c("sample1_trimmed.fastq", "sample2_trimmed.fastq"),
#'   output_dir = "trimmed_data",
#'   minReadLength = 20,
#'   thread = 4
#' )
#'
#' # Paired-end data
#' trim_adaptors(
#'   read1 = c("sample1_R1.fastq", "sample2_R1.fastq"),
#'   read2 = c("sample1_R2.fastq", "sample2_R2.fastq"),
#'   outputFastq = c("sample1_trimmed", "sample2_trimmed"),
#'   output_dir = "trimmed_data",
#'   minReadLength = 20,
#'   thread = 8,
#'   verbose = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link[Rfastp]{rfastp}} for the underlying trimming function
#'
#' @export
trim_adaptors <- function(read1 = NULL,
                          read2 = "",
                          output_dir = "2.trim-data",
                          outputFastq = NULL,
                          adapterSequenceRead1 = "auto",
                          adapterSequenceRead2 = "auto",
                          minReadLength = 15,
                          maxReadLength = 0,
                          thread = 1,
                          verbose = FALSE,
                          ...){
  # check read2
  if(read2 == ""){
    read2 <- rep("",length(read1))
  }

  # OUTPUT
  dir.create(path = output_dir)

  # loop
  # x = 1
  lapply(seq_along(read1),function(x){
    if (requireNamespace("Rfastp", quietly = TRUE)) {
      Rfastp::rfastp(read1 = read1[x],
                     read2 = read2[x],
                     outputFastq = paste0(output_dir,"/",outputFastq[x]),
                     adapterSequenceRead1 = adapterSequenceRead1,
                     adapterSequenceRead2 = adapterSequenceRead2,
                     minReadLength = minReadLength,
                     maxReadLength = maxReadLength,
                     thread = thread,
                     verbose = verbose,
                     ...)
    } else {
      warning("Package 'Rfastp' is needed for this function to work.")
    }

    return(message(paste(read1[x],"and",read2[x],"have been processed!")))
  })

}
