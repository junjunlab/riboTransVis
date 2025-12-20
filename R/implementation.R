


#' Process GTF File for Transcript Annotations
#'
#' Extracts transcript features from a GTF file, including exon length, CDS length,
#' and UTR lengths, and returns a processed data frame.
#'
#' @param gtf_file A `character` string specifying the path to a GTF annotation file,
#'   or a `data.frame` that has already been loaded from a GTF.
#' @param type_cds A `character` string for the feature type of coding sequences in the
#'   GTF file's 3rd column. Default is `"CDS"`.
#' @param type_utr5 A `character` string for the feature type of 5' UTRs.
#'   Default is `"five_prime_utr"`.
#' @param type_utr3 A `character` string for the feature type of 3' UTRs.
#'   Default is `"three_prime_utr"`.
#' @param type_exon A `character` string for the feature type of exons.
#'   Default is `"exon"`.
#' @param type_transcript A `character` string for the feature type of transcripts.
#'   Default is `"transcript"`.
#' @param transcript_attribute A `character` string specifying the attribute name for
#'   transcript identifiers in the GTF's 9th column. Common values are `"transcript_id"`
#'   (default for Gencode/Ensembl) or `"transcript_name"`.
#' @param gene_attribute A `character` string specifying the attribute name for gene
#'   identifiers. Common values are `"gene_name"` (default) or `"gene_id"`.
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
#' @importFrom dplyr filter group_by summarise mutate select left_join
#' @importFrom stats na.omit
#' @export
prepareTransInfo <- function(gtf_file = NULL,
                             type_cds = "CDS",
                             type_utr5 = "five_prime_utr",
                             type_utr3 = "three_prime_utr",
                             type_exon = "exon",
                             type_transcript = "transcript",
                             transcript_attribute = "transcript_id",
                             gene_attribute = "gene_name"){
  # load gtf
  if(is.character(gtf_file)){
    if (requireNamespace("rtracklayer", quietly = TRUE)) {
      gtf <- rtracklayer::import.gff(gtf_file,format = "gtf") %>%
        data.frame()
    } else {
      warning("Package 'rtracklayer' is needed for this function to work.")
    }
  }else if(inherits(gtf_file,"data.frame")){
    # check gene_name column
    if(!("gene_name" %in% colnames(gtf_file))){
      gtf_file$gene_name <- gtf_file$gene_id
    }

    gtf <- gtf_file
  }else {
    stop("'gtf_file' must be a file path or a data frame.")
  }

  # Check if specified attributes exist
  if (!transcript_attribute %in% names(gtf)) {
    stop(paste("Transcript attribute '", transcript_attribute, "' not found in GTF columns.", sep=""))
  }
  if (!gene_attribute %in% names(gtf)) {
    stop(paste("Gene attribute '", gene_attribute, "' not found in GTF columns.", sep=""))
  }

  # Use modern dplyr with .data pronoun for programming with variable names
  exon <- gtf %>%
    dplyr::filter(type == type_exon) %>%
    dplyr::group_by(.data[[transcript_attribute]]) %>%
    dplyr::summarise(exonlen = sum(width))

  st <- gtf %>%
    dplyr::filter(type == type_utr5) %>%
    dplyr::group_by(.data[[transcript_attribute]]) %>%
    dplyr::summarise(utr5 = sum(width))

  cds <- gtf %>%
    dplyr::filter(type == type_cds) %>%
    dplyr::group_by(.data[[transcript_attribute]]) %>%
    dplyr::summarise(cds = sum(width))

  sp <- gtf %>%
    dplyr::filter(type == type_utr3) %>%
    dplyr::group_by(.data[[transcript_attribute]]) %>%
    dplyr::summarise(utr3 = sum(width))

  idtran <- gtf %>%
    dplyr::filter(type == type_transcript) %>%
    # deal with no gene name symbol with  transcript_id
    dplyr::mutate(gene_name = ifelse(is.na(.data[[gene_attribute]]),
                  .data[[transcript_attribute]],.data[[gene_attribute]])) %>%
    dplyr::mutate(idnew = paste(.data[[transcript_attribute]],.data[[gene_attribute]],sep = "|")) %>%
    dplyr::select(.data[[transcript_attribute]],idnew) %>% unique() %>% na.omit()


  # combine info
  all_trans <- idtran %>%
    dplyr::left_join(y = st,by = transcript_attribute) %>%
    dplyr::mutate(utr5 = ifelse(is.na(utr5),0,utr5)) %>%
    dplyr::left_join(y = cds,by = transcript_attribute) %>%
    dplyr::left_join(y = sp,by = transcript_attribute) %>%
    dplyr::left_join(y = exon,by = transcript_attribute) %>%
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
# function to extract transcript sequence
# ==============================================================================

#' Extract Transcript Sequences from Genome and GTF Annotation
#'
#' This function extracts transcript (or other genomic feature) sequences based on genome FASTA and GTF annotations.
#' It also supports optional region extension upstream/downstream, and can output extended coordinates or the final transcript FASTA.
#'
#' @param genome_file Either a \code{BSgenome} object or the path to an uncompressed genome FASTA file (with accompanying \code{.fai} index or it will be generated).
#' @param gtf_file Path to a GTF annotation file containing gene models (required) or Granges format file from gtf data.
#' @param feature Feature type to extract sequences for (Default: "exon"). Could be "exon", "CDS", etc.
#' @param extend Logical; whether to extend the first and last exons per transcript (Default: \code{FALSE}).
#' @param extend_upstream Integer; number of bases to extend upstream of the first exon (Default: \code{0}).
#' @param extend_downstream Integer; number of bases to extend downstream of the last exon (Default: \code{0}).
#' @param return_extend_region Logical; if \code{TRUE}, returns the genomic coordinates (as a \code{data.frame}) of the extracted (possibly extended) regions instead of sequences (Default: \code{FALSE}).
#' @param output_file Output file path to write the transcript sequences in FASTA format. Required if \code{return_extend_region = FALSE}.
#'
#' @return If \code{return_extend_region = TRUE}, returns a \code{data.frame} containing the extended coordinates of the selected feature (e.g., exon). Otherwise, writes transcript sequences to \code{output_file} and returns nothing.
#'
#' @details
#' - The function supports both UCSC-style and Ensembl-style chromosome names, and will automatically adjust "chr" prefix if needed.
#' - Transcript sequences are reconstructed by collapsing all exons (or chosen \code{feature}) for each transcript.
#' - For \code{feature = "CDS"}, stop codons (+3bp) can optionally be included at the 3' end (or 5' end for minus strand).
#' - Requires the following Bioconductor packages: \code{rtracklayer}, \code{Rsamtools}, \code{Biostrings}, \code{txdbmaker}, \code{GenomicFeatures}, \code{AnnotationDbi}
#'
#' @examples
#' \dontrun{
#' get_transcript_sequence(
#'   genome_file = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#'   gtf_file = "Homo_sapiens.GRCh38.110.gtf",
#'   feature = "exon",
#'   extend = TRUE,
#'   extend_upstream = 50,
#'   extend_downstream = 100,
#'   return_extend_region = FALSE,
#'   output_file = "transcripts.fa"
#' )
#' }
#'
#' @importFrom Rsamtools FaFile indexFa
#' @importFrom GenomicFeatures extractTranscriptSeqs makeTxDbFromGRanges
#' @importFrom Biostrings writeXStringSet
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr %>% filter mutate if_else left_join n select
#' @importFrom fastplyr f_filter f_inner_join
#' @export
get_transcript_sequence <- function(genome_file = NULL,
                                    gtf_file = NULL,
                                    feature = "exon",
                                    extend = FALSE,
                                    extend_upstream = 0,
                                    extend_downstream = 0,
                                    return_extend_region = FALSE,
                                    output_file = NULL){
  # ============================================================================
  # load genome and gtf
  # ============================================================================
  if(!is.null(genome_file)){
    # check db
    if(inherits(genome_file,"BSgenome")){
      fa_file <- genome_file

    }else{
      # extact chromosome length
      if (!file.exists(paste0(genome_file, ".fai"))) {
        Rsamtools::indexFa(genome_file)
      }

      fa_file <- Rsamtools::FaFile(genome_file)
    }

    # chromosome length
    if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
      seqinfo <- GenomeInfoDb::seqinfo(fa_file)
      chr_len <- data.frame(seqnames = names(seqinfo),
                            len = as.numeric(GenomeInfoDb::seqlengths(seqinfo))
      )
    } else {
      warning("Package 'GenomeInfoDb' is needed for this function to work.")
    }

  }

  # load gtf file
  if(inherits(gtf_file,"GRanges")){
    gtf <- gtf_file %>% data.frame()
  }else if(is.character(gtf_file)){
    if (requireNamespace("rtracklayer", quietly = TRUE)) {
      gtf <- rtracklayer::import.gff(gtf_file,format = "gtf") %>% data.frame()
    } else {
      warning("Package 'rtracklayer' is needed for this function to work.")
    }
  }


  # ============================================================================
  # prepare info to extract
  # ============================================================================
  exon <- gtf %>% dplyr::filter(type == feature)

  # check chromosome prefix
  if(inherits(genome_file,"BSgenome") & !startsWith(as.character(exon$seqnames[1]),"chr")){
    exon$seqnames <- paste("chr",exon$seqnames,sep = "")
  }

  # check type
  if(feature != "exon"){
    exon$type <- "exon"
  }

  exon <- exon %>%
    # add tid for non gene symbol
    dplyr::mutate(gene_name = dplyr::if_else(is.na(gene_name),transcript_id,gene_name)) %>%
    # add new tid
    dplyr::mutate(idnew = paste(transcript_id,gene_name,sep = "|"))

  # check whether extend
  if(extend == TRUE){
    exon.final <- exon %>%
      # add exon index for each transcript
      dplyr::mutate(idx = dplyr::if_else(strand == "+",1:dplyr::n(),dplyr::n():1),
                    idx.mx = dplyr::n(),.by = transcript_id) %>%
      # extend upstream and downstream
      dplyr::mutate(start = dplyr::if_else(idx == 1 ,start - extend_upstream,start),
                    end = dplyr::if_else(idx == idx.mx, end + extend_downstream, end)) %>%
      # merge with chromosome length info
      fastplyr::f_inner_join(y = chr_len,by = "seqnames") %>%
      # check start and end boundary
      dplyr::mutate(start = dplyr::if_else(start <= 0,1,start),
                    end = dplyr::if_else(end > len,len,end)) %>%
      # redefine region width
      dplyr::mutate(width = abs(end - start) + 1)
  }else{
    exon.final <- exon
  }

  # check return what
  if(return_extend_region == TRUE){
    exon.final <- exon.final %>% dplyr::select(-idnew)
    return(exon.final)
  }else{
    # filter regions in genome
    if(extend == TRUE){
      exon.final <- exon.final %>%
        fastplyr::f_filter(end <= len)
    }else{
      exon.final <- exon.final %>%
        fastplyr::f_filter(seqnames %in% chr_len$seqnames) %>%
        fastplyr::f_inner_join(chr_len, by = "seqnames") %>%
        fastplyr::f_filter(end <= len)
    }


    # include stop codon
    if(feature == "CDS"){
      exon.final <- exon.final %>%
        # add exon index for each transcript
        dplyr::mutate(idx = dplyr::if_else(strand == "+",1:dplyr::n(),dplyr::n():1),
                      idx.mx = dplyr::n(),.by = transcript_id) %>%
        # extend upstream and downstream
        dplyr::mutate(start = dplyr::if_else(strand == "-" & idx == 1 ,start - 3,start),
                      end = dplyr::if_else(strand == "+" & idx == idx.mx, end + 3, end)) %>%
        # merge with chromosome length info
        fastplyr::f_inner_join(y = chr_len,by = "seqnames") %>%
        # redefine region width
        dplyr::mutate(width = abs(end - start) + 1)
    }

    # to txdb format
    if (requireNamespace("txdbmaker", quietly = TRUE)) {
      exon.final <- exon.final %>% dplyr::mutate(transcript_id = idnew)
      txdb.fl <- txdbmaker::makeTxDbFromGRanges(GenomicRanges::GRanges(exon.final))
    } else {
      warning("Package 'txdbmaker' is needed for this function to work.")
    }

    # id
    if (requireNamespace("AnnotationDbi", quietly = TRUE)) {
      tx2gene <- AnnotationDbi::select(txdb.fl,
                                       keys = AnnotationDbi::keys(txdb.fl, "TXNAME"),
                                       columns = c("TXNAME", "GENEID"),
                                       keytype = "TXNAME")
    } else {
      warning("Package 'AnnotationDbi' is needed for this function to work.")
    }


    # extract sequence
    totaltrans.seq <- GenomicFeatures::extractTranscriptSeqs(x = fa_file,transcripts = txdb.fl)

    # assign names
    names(totaltrans.seq) <- tx2gene$TXNAME


    # output
    Biostrings::writeXStringSet(x = totaltrans.seq,filepath = output_file,format = "fasta")
  }

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
      if (requireNamespace("zoo", quietly = TRUE)) {
        pos.df$smooth <- zoo::rollmean(pos.df$rpm, k = slide_window, fill = 0)
      } else {
        warning("Package 'zoo' is needed for this function to work.")
      }


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

  # check data
  if(nrow(pileup_result) == 0){
    pileup_result <- data.frame(seqnames = gene.ft$idnew, pos = 1, count = 0)
  }

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
                                 param = Rsamtools::ScanBamParam(what = c("rname", "pos", "strand", "qwidth"),
                                                                 which = region,
                                                                 flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
                                                                                               isSecondaryAlignment = FALSE,
                                                                                               isDuplicate = FALSE))
  )

  # list to data frame
  tmpdf <- purrr::map_df(seq_along(bam_data),function(x){
    tmp <- data.frame(bam_data[[x]]) %>%
      dplyr::group_by(rname,pos,qwidth,strand) %>%
      dplyr::summarise(count = dplyr::n(),.groups = "drop")

    return(tmp)
  })

  # check data
  if(nrow(tmpdf) == 0){
    tmpdf <- data.frame(seqnames = gene.ft$idnew, pos = 1, count = 0)
  }

  colnames(tmpdf)[1] <- "rname"

  # total mapped reads
  # total_reads <- sum(Rsamtools::idxstatsBam(bam_file)$mapped)

  # rpm normalization
  tmpdf$rpm <- (tmpdf$count/total_reads)*10^6

  return(tmpdf)
}
