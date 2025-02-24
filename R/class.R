globalVariables(c(".", ".I", "enrich", "exonlen", "gene", "gene_name", "idnew",
                  "mstart", "mstop", "n", "na.omit", "new","count", "counts", "frame", "qwidth",
                  "pos", "rname", "smooth", "smooth.x", "smooth.y", "transcript_id",
                  "translen", "type", "utr3","utr5", "which_label", "width",
                  "mappped_reads", "normval", "rpm","start","end","gene_id", "seqnames"))


#' ribotrans Class
#'
#' An S4 class to represent ribosome profiling data and associated information.
#'
#' @slot bam_file A `data.frame` containing BAM file data.
#' @slot library A `data.frame` representing the sequencing library information.
#' @slot features A `data.frame` describing transcript features (e.g., exons, genes).
#' @slot summary_info A `data.frame` describing all reads detailed information.
#' @slot gene_name A `character` string specifying the target gene name.
#' @slot ribo_occupancy A `data.frame` storing ribosome occupancy data.
#' @slot ribo.smoothed A `character` representing smoothed ribosome profiles.
#' @slot RNA_coverage A `data.frame` containing RNA coverage data from sequencing.
#' @slot RNA.smoothed A `character` representing smoothed RNA profiles.
#' @slot scaled_occupancy A `data.frame` with normalized ribosome occupancy data.
#'
#' @return An object of class \code{ribotrans}.
#' @export
ribotrans <- setClass("ribotrans",
                      slots = list("bam_file" = "data.frame",
                                   "library" = "data.frame",
                                   "features" = "data.frame",
                                   "summary_info" = "data.frame",
                                   "gene_name" = "character",
                                   "ribo_occupancy" = "data.frame",
                                   "ribo.smoothed" = "character",
                                   "RNA_coverage" = "data.frame",
                                   "RNA.smoothed" = "character",
                                   "scaled_occupancy" = "data.frame")
)



#' Construct a ribotrans Object
#'
#' Creates a `ribotrans` object from RNA and ribosome profiling BAM files along
#' with transcriptome annotation (`GTF` file).
#'
#' @param gtf_file A `character` string indicating the path to the GTF annotation file.
#' @param RNA_bam_file A `character` vector containing paths to RNA-seq BAM files.
#' @param RNA_sample_name A `character` vector representing RNA-seq sample names.
#' @param Ribo_bam_file A `character` vector containing paths to ribosome profiling BAM files.
#' @param Ribo_sample_name A `character` vector representing ribosome profiling sample names.
#' @param choose_longest_trans A logical value indicating whether to select only
#'   the longest transcript for each gene. Default is `FALSE`.
#'
#' @details
#' This function processes the transcriptome annotation and extracts alignment
#' information from RNA and ribosome profiling BAM files to create
#' an S4 object of class \code{ribotrans}.
#'
#' It computes total mapped reads for each BAM file and stores the information
#' in the `library` slot of the `ribotrans` object.
#'
#' @return An object of class \code{ribotrans} containing BAM file metadata,
#' library information, and parsed transcriptomic features.
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom data.table as.data.table
#' @importFrom methods new
#' @export
construct_ribotrans <- function(gtf_file = NULL,
                                RNA_bam_file = NULL,
                                RNA_sample_name = NULL,
                                Ribo_bam_file = NULL,
                                Ribo_sample_name = NULL,
                                choose_longest_trans = FALSE){
  # ============================================================================
  # prepare trans info
  # ============================================================================
  features <- prepareTransInfo(gtf_file = gtf_file)

  # whether select longest transcript
  if(choose_longest_trans == TRUE){
    dt_features <- data.table::as.data.table(features)

    features <- dt_features[
      dt_features[, .I[which.max(translen)], by = gene]$V1
    ]
  }
  # ============================================================================
  # extract library info
  # ============================================================================
  # total mapped reads
  bams <- c(RNA_bam_file, Ribo_bam_file)
  sp <- c(RNA_sample_name, Ribo_sample_name)
  tp <- rep(c("rna", "ribo"), c(length(RNA_bam_file),length(Ribo_bam_file)))

  lapply(seq_along(bams),function(x){
    total_reads <- sum(Rsamtools::idxstatsBam(bams[x])$mapped)

    data.frame(bam = bams[x],
               mappped_reads = total_reads,
               type = tp[x],
               sample = sp[x])
  }) %>% do.call("rbind", .) %>% data.frame() -> library


  # ============================================================================
  # create ribotrans object
  # ============================================================================
  bam_file <- data.frame(bam = bams,
                         type = rep(c("rna", "ribo"),
                                    c(length(RNA_bam_file),length(Ribo_bam_file))),
                         sample = c(RNA_sample_name, Ribo_sample_name)
  )

  res <- methods::new("ribotrans",
                      "bam_file" = bam_file,
                      "library" = library,
                      "features" = features)

  return(res)
}
