
#' Extract Transcript Information from a GTF File
#'
#' This function loads a GTF file and extracts exon information for each transcript.
#' It computes the total transcript length and determines the exon positions
#' relative to the transcript.
#'
#' @param gtf_file Character. Path to the GTF file.
#'
#' @return A `data.frame` containing transcript position information with the columns:
#' \itemize{
#'   \item `seqnames`: Chromosome name.
#'   \item `start`, `end`: Genomic coordinates of exons.
#'   \item `width`: Length of each exon.
#'   \item `strand`: Strand direction (`+` or `-`).
#'   \item `gene_name`: Gene name.
#'   \item `transcript_id`: Transcript ID.
#'   \item `tx_len`: Cumulative transcript length up to this exon.
#'   \item `f_len`: Total length of the transcript.
#' }
#'
#' @importFrom dplyr filter group_by summarise arrange mutate left_join select
#'
#' @examples
#' \dontrun{
#' gtf_file <- "path/to/your.gtf"
#' transcripts <- prepareTransInfo_forGenome(gtf_file)
#' head(transcripts)
#' }
#'
#' @export
prepareTransInfo_forGenome <- function(gtf_file = NULL){
  # load gtf
  if(is.character(gtf_file)){
    if (requireNamespace("rtracklayer", quietly = TRUE)) {
      gtf <- rtracklayer::import.gff(gtf_file,format = "gtf") %>%
        data.frame()
    } else {
      warning("Package 'rtracklayer' is needed for this function to work.")
    }
  }else if("data.frame" %in% class(gtf_file)){
    gtf <- gtf_file
  }


  # total transcript lenth
  features_len <- gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::group_by(gene_name,transcript_id) %>%
    dplyr::summarise(f_len = sum(width),.groups = "drop")

  # calculate transcript position block for positive strand
  exon_pos <- gtf %>%
    dplyr::filter(type == "exon" & strand == "+") %>%
    dplyr::select(seqnames,start,end,width,strand,gene_name,transcript_id) %>%
    dplyr::group_by(gene_name,transcript_id) %>%
    dplyr::arrange(start,end,.by_group = T) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width") %>%
    dplyr::left_join(y = features_len,by = c("gene_name","transcript_id"))

  # calculate transcript position block for negtive strand
  exon_neg <- gtf %>%
    dplyr::filter(type == "exon" & strand == "-") %>%
    dplyr::select(seqnames,start,end,width,strand,gene_name,transcript_id) %>%
    dplyr::group_by(gene_name,transcript_id) %>%
    dplyr::arrange(dplyr::desc(start),dplyr::desc(end),.by_group = T) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width") %>%
    dplyr::left_join(y = features_len,by = c("gene_name","transcript_id"))

  # combine + and - data
  dist_features <- rbind(exon_pos, exon_neg)

  return(dist_features)
}


# ==============================================================================
# function to calculate ribosome occupancy on genome
# ==============================================================================

#' Compute Read Occupancy on the Genome
#'
#' The `getOccupancyGenome()` function extracts read positions from a BAM file
#' for a given gene, converts them to transcriptome positions, and calculates
#' normalized read occupancy (RPM).
#'
#' @param bam_file Character. Path to the BAM file.
#' @param gene_name Character. The gene name for which read occupancy is extracted.
#' @param features A `data.frame` or `GRanges` object containing transcript annotation.
#' @param total_reads Numeric. The total number of mapped reads (for RPM normalization).
#' @param assignment_mode Character. Specifies the read assignment strategy: `"end5"` (5' end) or `"end3"` (3' end).
#' @param coordinate_to_trans Logical. Whether to convert genomic coordinates to transcript
#'        coordinates. Default is FALSE.
#' @param return_all_position Logical; if \code{TRUE}, extracts data for **all mapped regions**,
#' not just the specified gene (default: \code{FALSE}).
#'
#' @return A `data.frame` with read occupancy information:
#' \itemize{
#'   \item `rname`: Transcript and gene identifier (`"transcript_id|gene_name"`).
#'   \item `pos`: Transcriptomic read position.
#'   \item `count`: Number of reads at the position.
#'   \item `rpm`: Reads per million (RPM) normalized occupancy.
#' }
#'
#' @importFrom dplyr rename mutate group_by summarise filter select case_when
#' @importFrom fastplyr f_filter
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom Rsamtools scanBam ScanBamParam scanBamFlag
#' @importFrom parallel detectCores
#'
#' @examples
#' \dontrun{
#' bam_file <- "example.bam"
#' gene <- "TP53"
#' total_mapped_reads <- 1e7
#' features <- prepareTransInfo_forGenome(gtf_file = "example.gtf")
#'
#' occupancy_data <- getOccupancyGenome(
#'   bam_file = bam_file,
#'   gene_name = gene,
#'   features = features,
#'   total_reads = total_mapped_reads
#' )
#' head(occupancy_data)
#' }
#'
#' @export
getOccupancyGenome <- function(bam_file = NULL,
                               gene_name = NULL,
                               features = NULL,
                               total_reads = NULL,
                               assignment_mode = NULL,
                               coordinate_to_trans = FALSE,
                               return_all_position = FALSE){
  # check return type
  if(return_all_position == TRUE){
    # get position
    bam_data <- Rsamtools::scanBam(file = bam_file,
                                   nThreads = parallel::detectCores(),
                                   param = Rsamtools::ScanBamParam(what = c("rname", "pos", "strand", "qwidth"),
                                                                   flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)))
  }else{
    # filter genes
    query_region <- features %>%
      dplyr::rename(gene = gene_name) %>%
      fastplyr::f_filter(gene == gene_name) %>%
      GenomicRanges::GRanges()

    # get position
    bam_data <- Rsamtools::scanBam(file = bam_file,
                                   nThreads = parallel::detectCores(),
                                   param = Rsamtools::ScanBamParam(what = c("rname", "pos", "strand", "qwidth"),
                                                                   which = query_region,
                                                                   flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
                                                                                                 isSecondaryAlignment = FALSE,
                                                                                                 isDuplicate = FALSE)))

  }

  # list to data frame
  bminfo <- lapply(bam_data,FUN = data.frame) %>% do.call("rbind",.)

  if(assignment_mode == "end5"){
    # check strand to assign 5'end position
    bminfo <- bminfo %>%
      dplyr::mutate(pos = ifelse(strand == "-", pos + qwidth - 1, pos))
  }else{
    # check strand to assign 3'end position
    bminfo <- bminfo %>%
      dplyr::mutate(pos = ifelse(strand == "+", pos + qwidth - 1, pos))
  }

  # summarise
  bminfo <- bminfo %>%
    dplyr::group_by(rname,pos,strand,qwidth) %>%
    dplyr::summarise(count = dplyr::n(),.groups = "drop")

  if(return_all_position == FALSE){
    bminfo <- bminfo %>%
      dplyr::filter(pos >= min(IRanges::start(query_region@ranges)) &
                      pos <= max(IRanges::end(query_region@ranges)))
  }


  # ============================================================================
  # check coordinate transform
  if(coordinate_to_trans == TRUE){
    bminfo <- bminfo %>%
      dplyr::rename(seqnames = rname, start = pos) %>%
      dplyr::mutate(end = start,.after = "start")

    # to GRanges
    bminfo.rg <- GenomicRanges::GRanges(bminfo)

    # transcript annotation regions
    trans_rg <- GenomicRanges::GRanges(features)


    # genomic position to transcriptome position
    ov <- IRanges::findOverlaps(query = bminfo.rg,subject = trans_rg, ignore.strand = TRUE)

    # get overlap data
    if (requireNamespace("S4Vectors", quietly = TRUE)) {
      lo <- cbind(as.data.frame(bminfo.rg[S4Vectors::queryHits(ov)]),
                  as.data.frame(trans_rg[S4Vectors::subjectHits(ov)]))
    } else {
      warning("Package 'S4Vectors' is needed for this function to work.")
    }


    # make unique names
    names(lo) <- make.names(names(lo),unique = T)

    # calculate transcript position
    tgene <- gene_name
    lo <- lo %>%
      dplyr::mutate(pos = dplyr::case_when(strand == "+" & strand.1 == "+" ~ tx_len - abs(end.1 - start),
                                           # strand == "-" & strand.1 == "+" ~ tx_len - abs(end.1 - (start - qwidth + 1)),
                                           strand == "-" & strand.1 == "-" ~ tx_len - abs(start - start.1)
                                           # strand == "+" & strand.1 == "-" ~ tx_len - abs(start - (start + qwidth - 1))
                                           ),
                    rname = paste(transcript_id,gene_name,sep = "|")) %>%
      dplyr::filter(gene_name == tgene)

    # output
    tdf <- lo %>% dplyr::select(rname,pos,qwidth,count)
  }else{
    tdf <- bminfo
  }

  # rpm normalization
  tdf$rpm <- (tdf$count/total_reads)*10^6


  return(tdf)
}



# ==============================================================================
# function to calculate rna coverage on genome
# ==============================================================================

#' Compute RNA-seq Coverage on the Genome
#'
#' The `getCoverageGenome()` function extracts read coverage information from a BAM file
#' for a given gene, converts it to transcriptomic positions, and calculates normalized
#' coverage (RPM).
#'
#' @param bam_file Character. Path to the RNA-seq BAM file.
#' @param gene_name Character. The gene name for which transcript coverage is extracted.
#' @param features A `data.frame` or `GRanges` object containing transcript annotation.
#' @param total_reads Numeric. The total number of mapped reads (for RPM normalization).
#' @param coordinate_to_trans Logical. Whether to convert genomic coordinates to transcript
#'        coordinates. Default is FALSE.
#'
#' @return A `data.frame` with transcript coverage information:
#' \itemize{
#'   \item `rname`: Transcript and gene identifier (`"transcript_id|gene_name"`).
#'   \item `pos`: Transcriptomic position of the read coverage.
#'   \item `count`: Coverage (number of reads covering the position).
#'   \item `rpm`: Reads per million (RPM) normalized coverage.
#' }
#'
#' @importFrom dplyr rename mutate filter select case_when
#' @importFrom fastplyr f_filter
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom Rsamtools pileup BamFile PileupParam ScanBamParam
#'
#' @examples
#' \dontrun{
#' bam_file <- "example.bam"
#' gene <- "TP53"
#' total_mapped_reads <- 1e7
#' features <- prepareTransInfo_forGenome(gtf_file = "example.gtf")
#'
#' coverage_data <- getCoverageGenome(
#'   bam_file = bam_file,
#'   gene_name = gene,
#'   features = features,
#'   total_reads = total_mapped_reads
#' )
#' head(coverage_data)
#' }
#'
#' @export
getCoverageGenome <- function(bam_file = NULL,
                              gene_name = NULL,
                              features = NULL,
                              total_reads = NULL,
                              coordinate_to_trans = FALSE){
  # filter genes
  query_region <- features %>%
    dplyr::rename(gene = gene_name) %>%
    fastplyr::f_filter(gene == gene_name) %>%
    GenomicRanges::GRanges()

  # get coverage
  pileup_result <- Rsamtools::pileup(file = Rsamtools::BamFile(bam_file),
                                     pileupParam = Rsamtools::PileupParam(distinguish_nucleotides = FALSE,
                                                                          distinguish_strands = FALSE,
                                                                          max_depth = 10^6),
                                     scanBamParam = Rsamtools::ScanBamParam(which = query_region)) %>%
    dplyr::select(-which_label)

  colnames(pileup_result)[1] <- "rname"

  # ============================================================================
  # check coordinate transform
  if(coordinate_to_trans == TRUE){
    pileinfo <- pileup_result %>%
      dplyr::rename(seqnames = rname, start = pos) %>%
      dplyr::mutate(end = start,.after = "start")

    # to GRanges
    pileinfo.rg <- GenomicRanges::GRanges(pileinfo)

    # transcript annotation regions
    trans_rg <- GenomicRanges::GRanges(features)


    # genomic position to transcriptome position
    ov <- IRanges::findOverlaps(query = pileinfo.rg,subject = trans_rg, ignore.strand = TRUE)

    # get overlap data
    if (requireNamespace("S4Vectors", quietly = TRUE)) {
      lo <- cbind(as.data.frame(pileinfo.rg[S4Vectors::queryHits(ov)]),
                  as.data.frame(trans_rg[S4Vectors::subjectHits(ov)]))
    } else {
      warning("Package 'S4Vectors' is needed for this function to work.")
    }


    # make unique names
    names(lo) <- make.names(names(lo),unique = T)

    # calculate transcript position
    tgene <- gene_name
    lo <- lo %>%
      dplyr::mutate(pos = dplyr::case_when(strand.1 == "+" ~ tx_len - abs(end.1 - start),
                                           strand.1 == "-" ~ tx_len - abs(start - start.1)),
                    rname = paste(transcript_id,gene_name,sep = "|")) %>%
      dplyr::filter(gene_name == tgene)

    # output
    tdf <- lo %>%
      dplyr::select(rname,pos,count)
  }else{
    tdf <- pileup_result %>%
      dplyr::filter(pos >= min(IRanges::start(query_region@ranges)) &
                      pos <= max(IRanges::end(query_region@ranges)))
  }


  # rpm normalization
  tdf$rpm <- (tdf$count/total_reads)*10^6

  return(tdf)
}




