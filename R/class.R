globalVariables(c(".", ".I", "enrich", "exonlen", "gene", "gene_name", "idnew","allexp", "ave", "cds", "relexp",
                  "mstart", "mstop", "n", "na.omit", "new","count", "counts", "frame", "qwidth",
                  "pos", "rname", "smooth", "rpm.x", "rpm.y", "transcript_id", "strand", "codon",
                  "translen", "type", "utr3","utr5", "which_label", "width","f_len", "object", "len",
                  "mappped_reads", "normval", "rpm","start","end","gene_id", "seqnames","idx", "idx.mx"))


#' ribotrans Class
#'
#' An S4 class to represent ribosome profiling data and associated information.
#'
#' @slot bam_file A `data.frame` containing BAM file data.
#' @slot library A `data.frame` representing the sequencing library information.
#' @slot assignment_mode A `character` string indicating the methodology used for read assignment.
#' @slot mapping_type A `character` string specifying the mapping type (`"transcriptome"` or `"genome"`).
#' @slot features A `data.frame` describing transcript features (e.g., exons, genes).
#' @slot genome_trans_features A `data.frame` containing genomic transcript features for genome-based mappings.
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
                                   "mapping_type" = "character",
                                   "assignment_mode" = "character",
                                   "features" = "data.frame",
                                   "genome_trans_features" = "data.frame",
                                   "summary_info" = "data.frame",
                                   "gene_name" = "character",
                                   "ribo_occupancy" = "data.frame",
                                   "ribo.smoothed" = "character",
                                   "RNA_coverage" = "data.frame",
                                   "RNA.smoothed" = "character",
                                   "scaled_occupancy" = "data.frame")
)



#' Construct a Ribotrans Object for Transcriptome or Genome Mapping
#'
#' `construct_ribotrans()` generates a structured object with transcript or genomic annotations
#' and mapped RNA-seq / Ribo-seq library information. It supports selecting the longest transcript for each gene.
#'
#' @param gtf_file A `character` string indicating the path to the GTF annotation file.
#' @param mapping_type Character. Either `"transcriptome"` or `"genome"`, indicating whether
#'        RNA-seq/Ribo-seq alignments are transcriptome-based or genome-based. Default is `"transcriptome"`.
#' @param assignment_mode Character. Specifies the read assignment strategy: `"end5"` (5' end) or `"end3"` (3' end).
#' @param extend Logical. Whether to extend sequences upstream and downstream of transcripts.
#'              Default is FALSE.
#' @param extend_upstream Numeric. Number of bases to extend upstream of transcript start.
#'                       Only used when extend = TRUE. Default is 0.
#' @param extend_downstream Numeric. Number of bases to extend downstream of transcript end.
#'                         Only used when extend = TRUE. Default is 0.
#' @param RNA_bam_file A `character` vector containing paths to RNA-seq BAM files.
#' @param RNA_sample_name A `character` vector representing RNA-seq sample names.
#' @param Ribo_bam_file A `character` vector containing paths to ribosome profiling BAM files.
#' @param Ribo_sample_name A `character` vector representing ribosome profiling sample names.
#' @param choose_longest_trans Logical value indicating whether to select the longest transcript
#' for each gene. If TRUE, only the longest transcript (based on CDS and transcript length)
#' will be kept for each gene. This is useful to reduce redundancy when multiple transcript
#' isoforms exist for a gene. If FALSE (default), all transcripts will be retained.
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
#'
#' @examples
#' \dontrun{
#' ribotrans_obj <- construct_ribotrans(
#'   gtf_file = "example.gtf",
#'   mapping_type = "transcriptome",
#'   RNA_bam_file = c("rna1.bam", "rna2.bam"),
#'   RNA_sample_name = c("sample1", "sample2"),
#'   Ribo_bam_file = c("ribo1.bam"),
#'   Ribo_sample_name = c("sample3"),
#'   choose_longest_trans = TRUE
#' )
#' }
#'
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom methods new
#' @importFrom parallel detectCores
#' @export
construct_ribotrans <- function(gtf_file = NULL,
                                mapping_type = c("transcriptome", "genome"),
                                assignment_mode = c("end5", "end3"),
                                extend = FALSE,
                                extend_upstream = 0,
                                extend_downstream = 0,
                                RNA_bam_file = NULL,
                                RNA_sample_name = NULL,
                                Ribo_bam_file = NULL,
                                Ribo_sample_name = NULL,
                                choose_longest_trans = FALSE){
  mapping_type <- match.arg(mapping_type, choices = c("transcriptome", "genome"))
  assignment_mode <- match.arg(assignment_mode, choices = c("end5", "end3"))

  # ============================================================================
  # prepare trans info
  # ============================================================================
  # transcriptome features
  features <- prepareTransInfo(gtf_file = gtf_file)

  # check extend information
  if(extend == TRUE){
    features <- features %>%
      dplyr::mutate(utr5 = utr5 + extend_upstream,
                    utr3 = utr3 + extend_downstream,
                    exonlen = exonlen + extend_upstream + extend_downstream,
                    translen = translen + extend_upstream + extend_downstream,
                    mstart = dplyr::if_else(cds != 0,mstart + extend_upstream,mstart),
                    mstop = dplyr::if_else(cds != 0,mstop + extend_upstream,mstop))
  }

  # whether select longest transcript
  if(choose_longest_trans == TRUE){
    features <- features %>%
      fastplyr::f_group_by(gene) %>%
      fastplyr::f_arrange(cds, translen,.by_group = T,.descending = T) %>%
      fastplyr::f_slice_head(n = 1,keep_order = T)
  }

  # check mapping_type
  if(mapping_type == "genome"){

    if(extend == TRUE){
      gtf_input <- get_transcript_sequence(gtf_file = gtf_file,
                                           extend = extend,
                                           extend_upstream = extend_upstream,
                                           extend_downstream = extend_downstream,
                                           return_extend_region = TRUE)

      features.g <- prepareTransInfo_forGenome(gtf_file = gtf_input)
    }else{
      features.g <- prepareTransInfo_forGenome(gtf_file = gtf_file)
    }

    features.g2 <- subset(features.g, transcript_id %in% features$transcript_id)

    # if (requireNamespace("data.table", quietly = TRUE)) {
    #   data.table::setDT(features.g)
    #   data.table::setDT(features)
    #   features.g2 <- features.g[transcript_id %in% features$transcript_id]
    # } else {
    #   warning("Package 'data.table' is needed for this function to work.")
    # }

  }else{
    features.g <- data.frame()
  }
  # ============================================================================
  # extract library info
  # ============================================================================
  # total mapped reads
  bams <- c(RNA_bam_file, Ribo_bam_file)
  sp <- c(RNA_sample_name, Ribo_sample_name)
  tp <- rep(c("rna", "ribo"), c(length(RNA_bam_file),length(Ribo_bam_file)))

  lapply(seq_along(bams),function(x){
    # check index file for bam
    if(!file.exists(paste(bams[x],"bai",sep = "."))){
      Rsamtools::indexBam(bams[x], nThreads = parallel::detectCores())
    }

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
                      "features" = features,
                      "mapping_type" = mapping_type,
                      "assignment_mode" = assignment_mode,
                      "genome_trans_features" = features.g)

  return(res)
}
