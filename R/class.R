globalVariables(c(".", ".I", "enrich", "exonlen", "gene", "gene_name", "idnew","allexp", "ave", "cds", "relexp",
                  "mstart", "mstop", "n", "na.omit", "new","count", "counts", "frame", "qwidth", "group_id",
                  "pos", "rname", "smooth", "rpm.x", "rpm.y", "transcript_id", "strand", "codon", "rel_pos",
                  "translen", "type", "utr3","utr5", "which_label", "width","f_len", "object", "len", "chr",
                  "mappped_reads", "normval", "rpm","start","end","gene_id", "seqnames","idx", "idx.mx","rel_scale",
                  "all_counts", "periodicity","lower", "upper", "win_ip", "win_total", "xend","codon_pos",
                  "avg", "avg_ct", "cdslen", "norm_avg", "normsm","count.tt", "density", "label", "relsp",
                  "relst", "rpm.tt", "sp", "wi", "x", "y", "yend","reloccup", "value","avg_val",  "dist","rpks",
                  "Abbreviation1", "Abbreviation3", "AminoAcid", "abbrev", "codon_seq", "freq", "occup","rpk",
                  "occurrence", "pause_score", "pep_seq", "ratio", "rel_pause", "tripep_val","motif","nt_pos",
                  "sample_group","sum_pi","counts.y", "periodicity.x","codons","log2FoldChange","Amplitude",
                  "Period","mapped", "pvalue", "seq_type","ccdf","pos3sp", "pos3st", "pos5sp", "pos5st",
                  "sm1", "sm2", "smratio","bin", "nmft", "winstart","ci_high", "ci_low","norm.x","norm.y",
                  "above", "background_mean", "group", "ip", "tt","gene_biotype","max_sum", "min_sum",
                  "sum_1", "sum_2", "total_sum","pause_mean", "pause_sd", "total_counts","lp", "rp","ribo_count_feature"))


#' ribotrans Class
#'
#' An S4 class to represent ribosome profiling data and associated information.
#'
#' @slot bam_file A `data.frame` containing BAM file data.
#' @slot library A `data.frame` representing the sequencing library information.
#' @slot gtf_data A `GRanges` representing the gtf data for gene annotations.
#' @slot gtf_path A `character` string specifying the gtf file path.
#' @slot assignment_mode A `character` string indicating the methodology used for read assignment.
#' @slot mapping_type A `character` string specifying the mapping type (`"transcriptome"` or `"genome"`).
#' @slot features A `data.frame` describing transcript features (e.g., exons, genes).
#' @slot genome_trans_features A `data.frame` containing genomic transcript features for genome-based mappings.
#' @slot summary_info A `data.frame` describing all reads detailed information.
#' @slot reads_offset_info A \code{data.frame} recording information on ribosome P-site offsets.
#' @slot gene_name A `character` string specifying the target gene name.
#' @slot ribo_occupancy A `data.frame` storing ribosome occupancy data.
#' @slot ribo.smoothed A `character` representing smoothed ribosome profiles.
#' @slot RNA_coverage A `data.frame` containing RNA coverage data from sequencing.
#' @slot RNA.smoothed A `character` representing smoothed RNA profiles.
#' @slot scaled_occupancy A `data.frame` with normalized ribosome occupancy data.
#' @slot counts A list for saving RPF and mRNA counts data.
#'
#' @details
#' The `ribotrans` class is designed to facilitate the analysis of ribosome profiling data
#' by organizing relevant information, including BAM file metadata, gene annotations (GTF),
#' and various summary statistics related to ribosome and RNA-seq reads.
#'
#' Key functionalities include:
#' - Storing and accessing essential metadata
#' - Aggregating ribosome profiling signals
#' - Computing P-site offsets
#' - Normalizing ribosome occupancy with RNA-seq coverage
#'
#'
#' @return An object of class \code{ribotrans}.
#' @export
ribotrans <- setClass("ribotrans",
                      slots = list("bam_file" = "data.frame",
                                   "library" = "data.frame",
                                   "gtf_data" = "GRanges",
                                   "gtf_path" = "character",
                                   "mapping_type" = "character",
                                   "assignment_mode" = "character",
                                   "features" = "data.frame",
                                   "genome_trans_features" = "data.frame",
                                   "summary_info" = "data.frame",
                                   "reads_offset_info" = "data.frame",
                                   "gene_name" = "character",
                                   "ribo_occupancy" = "data.frame",
                                   "ribo.smoothed" = "character",
                                   "RNA_coverage" = "data.frame",
                                   "RNA.smoothed" = "character",
                                   "scaled_occupancy" = "data.frame",
                                   "counts" = "list")
)



#' Construct a Ribotrans Object for Transcriptome or Genome Mapping
#'
#' `construct_ribotrans()` generates a structured object with transcript or genomic annotations
#' and mapped RNA-seq / Ribo-seq library information. It supports selecting the longest transcript for each gene.
#'
#' @param genome_file A `character` string indicating the path to the genome sequence file
#' if the extend option is `TRUE`. Default is NULL.
#' @param gtf_file A `character` string indicating the path to the GTF annotation file.
#' @param mapping_type Character. Either `"transcriptome"` or `"genome"`, indicating whether
#'        RNA-seq/Ribo-seq alignments are transcriptome-based or genome-based. Default is `"transcriptome"`.
#'
#' @param assignment_mode Character. Specifies the read assignment strategy. Choices are:
#' \itemize{
#' \item{\code{"end5"}: Assign reads based on the 5' end of each aligned read.}
#' \item{\code{"end3"}: Assign reads based on the 3' end of each aligned read.}
#' }
#'
#'
#'
#' @param extend Logical. Whether to extend sequences upstream and downstream of transcripts.
#'              Default is FALSE.
#' @param extend_upstream Numeric. Number of bases to extend upstream of transcript start.
#'                       Only used when extend = TRUE. Default is 0.
#' @param extend_downstream Numeric. Number of bases to extend downstream of transcript end.
#'                         Only used when extend = TRUE. Default is 0.
#' @param RNA_bam_file A `character` vector containing paths to RNA-seq BAM files.
#' @param RNA_sample_name A `character` vector representing RNA-seq sample names.
#' @param RNA_sample_group Character vector. Group labels or conditions for RNA-seq samples.
#' @param Ribo_bam_file A `character` vector containing paths to ribosome profiling BAM files.
#' @param Ribo_sample_name A `character` vector representing ribosome profiling sample names.
#' @param Ribo_sample_group Character vector. Group labels or conditions for Ribo-seq samples.
#' @param choose_longest_trans Logical value indicating whether to select the longest transcript
#' for each gene. If TRUE, only the longest transcript (based on CDS and transcript length)
#' will be kept for each gene. This is useful to reduce redundancy when multiple transcript
#' isoforms exist for a gene. If FALSE (default), all transcripts will be retained.
#' @param prepareTransInfo_params A list parameters to prepareTransInfo function.
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
#' @importFrom utils modifyList
#' @export
construct_ribotrans <- function(genome_file = NULL,
                                gtf_file = NULL,
                                mapping_type = c("transcriptome", "genome"),
                                assignment_mode = c("end5", "end3"),
                                extend = FALSE,
                                extend_upstream = 0,
                                extend_downstream = 0,
                                RNA_bam_file = NULL,
                                RNA_sample_name = NULL,
                                RNA_sample_group = NULL,
                                Ribo_bam_file = NULL,
                                Ribo_sample_name = NULL,
                                Ribo_sample_group = NULL,
                                choose_longest_trans = FALSE,
                                prepareTransInfo_params = list()){
  options(fastplyr.inform = FALSE)
  options(fastplyr.optimise = FALSE)

  mapping_type <- match.arg(mapping_type, choices = c("transcriptome", "genome"))
  assignment_mode <- match.arg(assignment_mode, choices = c("end5", "end3"))

  # ============================================================================
  # prepare trans info
  # ============================================================================
  # load gtf
  if(is.character(gtf_file)){
    if (requireNamespace("rtracklayer", quietly = TRUE)) {
      gtf <- rtracklayer::import.gff(gtf_file,format = "gtf")
    } else {
      warning("Package 'rtracklayer' is needed for this function to work.")
    }
  }else if(inherits(gtf_file,"data.frame")){
    gtf <- gtf_file
  }else {
    stop("'gtf_file' must be a file path or a data frame.")
  }

  # check gene_name column
  if(!("gene_name" %in% colnames(gtf))){
    gtf$gene_name <- gtf$gene_id
  }

  # transcriptome features
  features <- do.call(prepareTransInfo,
                      modifyList(list(gtf_file = gtf_file),
                                 prepareTransInfo_params))

  # features <- prepareTransInfo(gtf_file = gtf_file)

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
      dplyr::group_by(gene) %>%
      dplyr::arrange(dplyr::desc(cds), dplyr::desc(translen),.by_group = T) %>%
      fastplyr::f_slice_head(n = 1,keep_order = T) %>%
      data.frame()
  }

  # check mapping_type
  if(mapping_type == "genome"){

    if(extend == TRUE){
      gtf_input <- get_transcript_sequence(gtf_file = gtf_file,
                                           genome_file = genome_file,
                                           extend = extend,
                                           extend_upstream = extend_upstream,
                                           extend_downstream = extend_downstream,
                                           return_extend_region = TRUE)

      features.g <- prepareTransInfo_forGenome(gtf_file = gtf_input)
    }else{
      features.g <- prepareTransInfo_forGenome(gtf_file = gtf_file)
    }

    features.g <- subset(features.g, transcript_id %in% features$transcript_id)

  }else{
    features.g <- data.frame()
  }
  # ============================================================================
  # extract library info
  # ============================================================================
  # total mapped reads
  bams <- c(RNA_bam_file, Ribo_bam_file)
  sp <- c(RNA_sample_name, Ribo_sample_name)

  RNA_sample_group <- if (is.null(RNA_sample_group)) RNA_sample_name else RNA_sample_group
  Ribo_sample_group <- if (is.null(Ribo_sample_group)) Ribo_sample_name else Ribo_sample_group
  gp <- c(RNA_sample_group, Ribo_sample_group)

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
               sample = sp[x],
               sample_group = gp[x])
  }) %>% do.call("rbind", .) %>% data.frame() -> library


  # ============================================================================
  # create ribotrans object
  # ============================================================================
  bam_file <- data.frame(bam = bams,
                         type = rep(c("rna", "ribo"),
                                    c(length(RNA_bam_file),length(Ribo_bam_file))),
                         sample = c(RNA_sample_name, Ribo_sample_name),
                         sample_group = gp
  )

  res <- methods::new("ribotrans",
                      "bam_file" = bam_file,
                      "gtf_data" = gtf,
                      "gtf_path" = gtf_file,
                      "library" = library,
                      "features" = features,
                      "mapping_type" = mapping_type,
                      "assignment_mode" = assignment_mode,
                      "genome_trans_features" = features.g)

  return(res)
}



# ==============================================================================
# serp class
# ==============================================================================


#' SERP (Selective Ribosome Profiling) Class
#'
#' @description
#' An S4 class for representing and analyzing Selective Ribosome Profiling data.
#' The class extends the 'ribotrans' class and provides specific slots for storing
#' total occupancy, IP occupancy, and enriched ratio data.
#'
#' @slot total_occupancy A data frame containing ribosome occupancy data from total samples.
#' @slot ip_occupancy A data frame containing ribosome occupancy data from IP (immunoprecipitation) samples.
#' @slot enriched_ratio A data frame containing the calculated enrichment ratios between IP and total samples.
#'
#' @details
#' The SERP class is designed for analyzing selective ribosome profiling experiments,
#' which typically involve comparing IP-enriched ribosome footprints with total ribosome footprints.
#' This class inherits all slots from the 'ribotrans' class and adds specific slots for
#' storing and analyzing the enrichment patterns.
#'
#' @seealso \code{\link{construct_serp}} for creating instances of this class.
#'
#'
#' @export
serp <- setClass("serp",
                 slots = list("total_occupancy" = "data.frame",
                              "ip_occupancy" = "data.frame",
                              "enriched_ratio" = "data.frame"),
                 contains = "ribotrans"
)





# ==============================================================================
# construct_serp function
# ==============================================================================


#' Construct a Ribotrans Object for Transcriptome or Genome Mapping
#'
#' `construct_ribotrans()` generates a structured object with transcript or genomic annotations
#' and mapped RNA-seq / Ribo-seq library information. It supports selecting the longest transcript for each gene.
#'
#' @param genome_file A `character` string indicating the path to the genome sequence file
#' if the extend option is `TRUE`. Default is NULL.
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
#' @param total_bam_file Character vector. Path(s) to BAM file(s) for total ribosome profiling samples.
#' @param total_sample_name Character vector. Sample name(s) for total ribosome profiling samples.
#' @param total_sample_group Character vector. Group labels or conditions for total ribosome profiling samples.
#' @param IP_bam_file Character vector. Path(s) to BAM file(s) for IP ribosome profiling samples.
#' @param IP_sample_name Character vector. Sample name(s) for IP ribosome profiling samples.
#' @param IP_sample_group Character vector. Group labels or conditions for IP ribosome profiling samples.
#' @param choose_longest_trans Logical value indicating whether to select the longest transcript
#' for each gene. If TRUE, only the longest transcript (based on CDS and transcript length)
#' will be kept for each gene. This is useful to reduce redundancy when multiple transcript
#' isoforms exist for a gene. If FALSE (default), all transcripts will be retained.
#' @return An object of class 'serp' containing the processed data.
#'
#' @details
#' The function performs several key steps:
#' 1. Loads and processes the GTF file to extract transcript information
#' 2. Applies any requested extensions to transcript regions
#' 3. Optionally selects the longest transcript per gene
#' 4. Prepares genome-based transcript features if genome mapping is specified
#' 5. Extracts library information from BAM files
#' 6. Creates and returns a new SERP object with all the processed information
#'
#' The function requires the 'rtracklayer' package for GTF file processing and
#' 'Rsamtools' for BAM file handling.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' # Extended usage with transcript extension
#' serp_obj <- construct_serp(
#'   genome_file = "path/to/genome.fa",
#'   gtf_file = "path/to/annotations.gtf",
#'   mapping_type = "genome",
#'   assignment_mode = "end5",
#'   extend = TRUE,
#'   extend_upstream = 50,
#'   extend_downstream = 50,
#'   total_bam_file = c("path/to/total_rep1.bam", "path/to/total_rep2.bam"),
#'   total_sample_name = c("total_rep1", "total_rep2"),
#'   IP_bam_file = c("path/to/IP_rep1.bam", "path/to/IP_rep2.bam"),
#'   IP_sample_name = c("IP_rep1", "IP_rep2"),
#'   choose_longest_trans = TRUE
#' )
#' }
#'
#'
#' @export
construct_serp <- function(genome_file = NULL,
                           gtf_file = NULL,
                           mapping_type = c("transcriptome", "genome"),
                           assignment_mode = c("end5", "end3"),
                           extend = FALSE,
                           extend_upstream = 0,
                           extend_downstream = 0,
                           total_bam_file = NULL,
                           total_sample_name = NULL,
                           total_sample_group = NULL,
                           IP_bam_file = NULL,
                           IP_sample_name = NULL,
                           IP_sample_group = NULL,
                           choose_longest_trans = FALSE){
  options(fastplyr.inform = FALSE)
  options(fastplyr.optimise = FALSE)

  mapping_type <- match.arg(mapping_type, choices = c("transcriptome", "genome"))
  assignment_mode <- match.arg(assignment_mode, choices = c("end5", "end3"))

  # ============================================================================
  # prepare trans info
  # ============================================================================
  # load gtf
  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    gtf <- rtracklayer::import.gff(gtf_file,format = "gtf")
  } else {
    warning("Package 'rtracklayer' is needed for this function to work.")
  }


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
      fastplyr::f_slice_head(n = 1,keep_order = T) %>%
      data.frame()
  }

  # check mapping_type
  if(mapping_type == "genome"){

    if(extend == TRUE){
      gtf_input <- get_transcript_sequence(gtf_file = gtf_file,
                                           genome_file = genome_file,
                                           extend = extend,
                                           extend_upstream = extend_upstream,
                                           extend_downstream = extend_downstream,
                                           return_extend_region = TRUE)

      features.g <- prepareTransInfo_forGenome(gtf_file = gtf_input)
    }else{
      features.g <- prepareTransInfo_forGenome(gtf_file = gtf_file)
    }

    features.g <- subset(features.g, transcript_id %in% features$transcript_id)

  }else{
    features.g <- data.frame()
  }
  # ============================================================================
  # extract library info
  # ============================================================================
  # total mapped reads
  bams <- c(total_bam_file, IP_bam_file)
  sp <- c(total_sample_name, IP_sample_name)

  total_sample_group <- if (is.null(total_sample_group)) total_sample_name else total_sample_group
  IP_sample_group <- if (is.null(IP_sample_group)) IP_sample_name else IP_sample_group
  gp <- c(total_sample_group, IP_sample_group)

  tp <- rep(c("total", "ip"), c(length(total_bam_file),length(IP_bam_file)))

  # summary
  # x = 1
  lapply(seq_along(bams),function(x){
    # check index file for bam
    if(!file.exists(paste(bams[x],"bai",sep = "."))){
      Rsamtools::indexBam(bams[x], nThreads = parallel::detectCores())
    }

    total_reads <- sum(Rsamtools::idxstatsBam(bams[x])$mapped)

    data.frame(bam = bams[x],
               mappped_reads = total_reads,
               type = tp[x],
               sample = sp[x],
               sample_group = gp[x])
  }) %>% do.call("rbind", .) %>% data.frame() -> library


  # ============================================================================
  # create ribotrans object
  # ============================================================================
  bam_file <- data.frame(bam = bams,
                         type = rep(c("total", "ip"),
                                    c(length(total_bam_file),length(IP_bam_file))),
                         sample = c(total_sample_name, IP_sample_name),
                         sample_group = gp
  )

  res <- methods::new("serp",
                      "bam_file" = bam_file,
                      "gtf_data" = gtf,
                      "library" = library,
                      "features" = features,
                      "mapping_type" = mapping_type,
                      "assignment_mode" = assignment_mode,
                      "genome_trans_features" = features.g)

  return(res)
}

