


#' Process GTF File for Transcript Annotations
#'
#' Extracts transcript features from a GTF file, including exon length, CDS length,
#' and UTR lengths, and returns a processed data frame.
#'
#' @param gtf_file A `character` string specifying the path to a GTF annotation file.
#'
#' @details
#' This function reads a GTF file and extracts relevant transcript features:
#'
#' - Exon lengths (`exonlen`)
#' - 5' UTR lengths (`utr5`)
#' - CDS lengths (`cds`)
#' - 3' UTR lengths (`utr3`)
#'
#' It also calculates total transcript lengths and determines start/stop
#' positions for coding sequences. Transcripts are assigned new unique IDs
#' combining `transcript_id` and `gene_name`.
#'
#' @return A `data.frame` with the following columns:
#' \describe{
#'   \item{transcript_id}{Character: The transcript ID.}
#'   \item{idnew}{Character: Combined ID (`transcript_id|gene_name`).}
#'   \item{utr5}{Integer: Length of the 5' UTR.}
#'   \item{cds}{Integer: Length of the coding sequence (CDS).}
#'   \item{utr3}{Integer: Length of the 3' UTR.}
#'   \item{exonlen}{Integer: Total exon length.}
#'   \item{translen}{Integer: Full transcript length including UTRs and CDS.}
#'   \item{mstart}{Integer: Start position of CDS within the transcript.}
#'   \item{mstop}{Integer: Stop position of CDS within the transcript.}
#'   \item{gene}{Character: Gene name associated with the transcript.}
#' }
#'
#' @importFrom rtracklayer import.gff
#' @importFrom dplyr filter group_by summarise mutate select left_join
#' @importFrom stats na.omit
#' @export
prepareTransInfo <- function(gtf_file = NULL){
  # load gtf
  gtf <- rtracklayer::import.gff(gtf_file,format = "gtf") %>%
    data.frame()

  # filter type
  exon <- gtf %>%
    dplyr::filter(type == "exon") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(exonlen = sum(width))

  st <- gtf %>%
    dplyr::filter(type == "five_prime_utr") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(utr5 = sum(width))

  cds <- gtf %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(cds = sum(width))

  sp <- gtf %>%
    dplyr::filter(type == "three_prime_utr") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(utr3 = sum(width))

  idtran <- gtf %>%
    dplyr::filter(type == "transcript") %>%
    # deal with no gene name symbol with  transcript_id
    dplyr::mutate(gene_name = ifelse(is.na(gene_name),transcript_id,gene_name)) %>%
    dplyr::mutate(idnew = paste(transcript_id,gene_name,sep = "|")) %>%
    dplyr::select(transcript_id,idnew) %>% unique() %>% na.omit()


  # combine info
  all_trans <- idtran %>%
    dplyr::left_join(y = st,by = "transcript_id") %>%
    dplyr::mutate(utr5 = ifelse(is.na(utr5),0,utr5)) %>%
    dplyr::left_join(y = cds,by = "transcript_id") %>%
    dplyr::left_join(y = sp,by = "transcript_id") %>%
    dplyr::left_join(y = exon,by = "transcript_id") %>%
    dplyr::mutate(utr3 = ifelse(is.na(utr3),0,utr3)) %>%
    # na.omit() %>%
    dplyr::mutate(translen = ifelse(is.na(cds),exonlen,utr5 + cds + utr3 + 3),
                  mstart = ifelse(is.na(cds),NA,utr5 + 1),
                  mstop = ifelse(is.na(cds),NA,utr5 + cds))

  all_trans[is.na(all_trans)] <- 0

  # add gene name
  all_trans$gene <- sapply(strsplit(all_trans$idnew,split = "\\|"),"[",2)

  return(all_trans)
}





# ==============================================================================
# smooth data for each position
# ==============================================================================


#' Smooth Data for Each Transcript Position
#'
#' Applies a sliding window smoothing function to mapped read counts (RPM) across transcript positions.
#'
#' @param features A `data.frame` containing transcript annotations.
#' @param posdf A `data.frame` containing positional read count information with columns:
#'   \describe{
#'     \item{rname}{Character: Transcript name in the format `"transcript_id|gene_name"`.}
#'     \item{pos}{Integer: Position along the transcript where reads are mapped.}
#'     \item{sample}{Character: Sample identifier.}
#'     \item{rpm}{Numeric: Reads per million mapped reads at each position.}
#'   }
#' @param slide_window An `integer` defining the window size for the rolling mean (default: `30`).
#'
#' @details
#' This function extracts transcript-specific mapped read counts from `posdf`,
#' aligns them with transcript annotations in `features`, and applies a moving
#' average smoothing function (`zoo::rollmean`) using a **sliding window approach**.
#'
#' It processes each transcript separately, ensuring that all positions (1 to
#' exon length) are covered, filling missing values with `0` before smoothing.
#'
#' @return A `data.frame` containing smoothed transcript positions with columns:
#' \describe{
#'   \item{rname}{Character: Transcript name (`transcript_id|gene_name`).}
#'   \item{pos}{Integer: Position within the transcript.}
#'   \item{sample}{Character: Sample name.}
#'   \item{rpm}{Numeric: Original reads per million (RPM) values.}
#'   \item{smooth}{Numeric: Smoothed RPM values using a rolling mean filter.}
#' }
#'
#' @importFrom dplyr filter left_join
#' @importFrom purrr map_df
#' @importFrom zoo rollmean
#' @export
smoothEachPosition <- function(features = NULL, posdf= NULL,slide_window = 30){
  genesymbol <- sapply(strsplit(as.character(posdf$rname[1]), split = "\\|"),"[",2)
  tanno <- subset(features, gene == genesymbol)

  # loop for each transcript
  purrr::map_df(1:nrow(tanno),function(x){
    tmp <- tanno[x,]

    tmp2 <- posdf %>%
      dplyr::filter(rname == tmp$idnew)

    sp <- unique(tmp2$sample)

    # loop for each sample
    purrr::map_df(seq_along(sp),function(x){
      tmp3 <- tmp2 %>%
        dplyr::filter(sample == sp[x])

      # each pos
      pos.df <- data.frame(rname = tmp3$rname[1],pos = 1:tmp$exonlen)

      # merge with value
      pos.df <- pos.df %>%
        dplyr::left_join(y = tmp3,by = c("rname", "pos"))
      pos.df$sample <- tmp3$sample[1]
      pos.df[is.na(pos.df)] <- 0

      # smooth data
      pos.df$smooth <- zoo::rollmean(pos.df$rpm, k = slide_window, fill = 0)

      return(pos.df)
    }) -> smooth.df

    return(smooth.df)
  }) -> pldf

  return(pldf)
}



# ==============================================================================
# function to calculate coverage
# ==============================================================================

#' Extract Coverage Information from a BAM File
#'
#' Computes read coverage for a given gene from an RNA-seq or ribosome profiling BAM file.
#'
#' @param bam_file A `character` string specifying the path to a BAM file.
#' @param gene_name A `character` string indicating the gene name to extract coverage for.
#' @param features A `data.frame` containing transcript annotations, including:
#'   \describe{
#'     \item{gene}{Character: Gene name.}
#'     \item{idnew}{Character: Combined transcript ID (`transcript_id|gene_name`).}
#'     \item{translen}{Integer: Total transcript length.}
#'   }
#' @param total_reads A `numeric` value representing the total number of mapped reads
#'     in the BAM file (used for RPM normalization).
#'
#' @details
#' This function extracts per-position read coverage for the given `gene_name` by:
#'
#' 1. Filtering the `features` data frame to match the target gene.
#' 2. Defining the transcript region (`GRanges` format).
#' 3. Using `Rsamtools::pileup()` to compute coverage at each position.
#' 4. Normalizing counts to Reads Per Million (RPM).
#'
#' RPM is computed as:
#' \deqn{RPM = \frac{\text{count at position}}{\text{total_reads}} \times 10^6}
#'
#' @return A `data.frame` with per-position read coverage, containing:
#' \describe{
#'   \item{rname}{Character: Transcript name (`transcript_id|gene_name`).}
#'   \item{pos}{Integer: Position along the transcript.}
#'   \item{count}{Integer: Raw read count at this position.}
#'   \item{strand}{Character: Strand information (if applicable).}
#'   \item{rpm}{Numeric: Normalized read count (Reads Per Million, RPM).}
#' }
#'
#' @importFrom Rsamtools pileup BamFile PileupParam ScanBamParam
#' @importFrom dplyr filter select
#' @importFrom GenomicRanges GRanges
#' @export
getCoverage <- function(bam_file = NULL,
                        gene_name = NULL,
                        features = NULL,
                        total_reads = NULL){
  # filter genes
  gene.ft <- features %>%
    dplyr::filter(gene == gene_name)

  # interval
  region <- data.frame(chr = gene.ft$idnew,
                       start = 1,
                       end = gene.ft$translen) %>%
    GenomicRanges::GRanges()

  # get coverage
  pileup_result <- Rsamtools::pileup(file = Rsamtools::BamFile(bam_file),
                                     pileupParam = Rsamtools::PileupParam(distinguish_nucleotides = FALSE,
                                                                          distinguish_strands = FALSE,
                                                                          max_depth = 10^6),
                                     scanBamParam = Rsamtools::ScanBamParam(which = region)) %>%
    dplyr::select(-which_label)

  colnames(pileup_result)[1] <- "rname"

  # total mapped reads
  # total_reads <- sum(Rsamtools::idxstatsBam(bam_file)$mapped)

  # rpm normalization
  pileup_result$rpm <- (pileup_result$count/total_reads)*10^6

  return(pileup_result)
}





# ==============================================================================
# function to calculate ribosome occupancy
# ==============================================================================


#' Compute Read Occupancy from a BAM File
#'
#' Extracts per-position read occupancy from an RNA-seq or ribosome profiling BAM file
#' for a given gene and normalizes the data to Reads Per Million (RPM).
#'
#' @param bam_file A `character` string specifying the path to a BAM file.
#' @param gene_name A `character` string specifying the gene name whose read occupancy should be extracted.
#' @param features A `data.frame` containing transcript annotation information, including:
#'   \describe{
#'     \item{gene}{Character: Gene name.}
#'     \item{idnew}{Character: Combined transcript ID (`transcript_id|gene_name`).}
#'     \item{translen}{Integer: Transcript length.}
#'   }
#' @param total_reads A `numeric` value indicating the total number of mapped reads
#'     in the BAM file (used for RPM normalization).
#'
#' @details
#' This function extracts **read occupancy** at each transcript position by:
#'
#' - Filtering the `features` table to extract transcript information.
#' - Creating a genomic interval (`GRanges`) corresponding to the transcript.
#' - Using `Rsamtools::scanBam()` to extract alignment start positions.
#' - Aggregating read counts per position.
#' - **Normalizing to Reads Per Million (RPM)**:
#'
#' \deqn{RPM = \frac{\text{count at position}}{\text{total_reads}} \times 10^6}
#'
#' @return A `data.frame` containing per-position occupancy data:
#' \describe{
#'   \item{rname}{Character: Transcript name (`transcript_id|gene_name`).}
#'   \item{pos}{Integer: Position along the transcript.}
#'   \item{count}{Integer: Number of reads that start at this position.}
#'   \item{rpm}{Numeric: Normalized read occupancy (Reads Per Million, RPM).}
#' }
#'
#' @importFrom Rsamtools scanBam ScanBamParam scanBamFlag
#' @importFrom dplyr filter group_by summarise
#' @importFrom purrr map_df
#' @importFrom GenomicRanges GRanges
#' @export
getOccupancy <- function(bam_file = NULL,
                         gene_name = NULL,
                         features = NULL,
                         total_reads = NULL){
  # filter genes
  gene.ft <- features %>%
    dplyr::filter(gene == gene_name)

  # interval
  region <- data.frame(chr = gene.ft$idnew,
                       start = 1,
                       end = gene.ft$translen) %>%
    GenomicRanges::GRanges()

  # get position
  bam_data <- Rsamtools::scanBam(file = bam_file,
                                 nThreads = parallel::detectCores(),
                                 param = Rsamtools::ScanBamParam(what = c("rname", "pos"),
                                                                 which = region,
                                                                 flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
  )

  # list to data frame
  tmpdf <- purrr::map_df(seq_along(bam_data),function(x){
    tmp <- data.frame(bam_data[[x]]) %>%
      dplyr::group_by(rname,pos) %>%
      dplyr::summarise(count = dplyr::n(),.groups = "drop")

    return(tmp)
  })

  # total mapped reads
  # total_reads <- sum(Rsamtools::idxstatsBam(bam_file)$mapped)

  # rpm normalization
  tmpdf$rpm <- (tmpdf$count/total_reads)*10^6

  return(tmpdf)
}
