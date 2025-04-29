#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


# Function to calculate probability to odds
prob2odds <- function(p){p / ( 1 - p )}



#' Extract sequences of longest transcripts from a GTF file
#'
#' @description
#' This function extracts the sequences of the longest transcript for each gene from a GTF file
#' using a reference genome. It first loads the GTF file, selects the longest transcript for each gene
#' based on CDS and transcript length, and then extracts the corresponding sequences from the genome file.
#'
#' @param gtf_file Path to a GTF annotation file.
#' @param genome_file Path to a reference genome file or BSgenome object .
#' @param output_file Path where the extracted transcript sequences will be saved.
#'
#' @details
#' The function follows these steps:
#' \enumerate{
#'   \item Loads the GTF file using rtracklayer
#'   \item Processes transcript information using prepareTransInfo()
#'   \item For each gene, selects the longest transcript prioritizing CDS length first, then total transcript length
#'   \item Filters the GTF to keep only the selected transcripts
#'   \item Extracts sequences for the selected transcripts from the reference genome
#' }
#'
#' @note
#' This function requires the 'rtracklayer' package for importing GTF files and uses
#' the internal functions 'prepareTransInfo' and 'get_transcript_sequence'.
#'
#' @return
#' Writes the transcript sequences to the specified output file and longest transcript information.
#'
#' @examples
#' \dontrun{
#' # Extract longest transcript sequences and save to output file
#' get_longest_transcript(
#'   gtf_file = "path/to/annotation.gtf",
#'   genome_file = "path/to/genome.fa",
#'   output_file = "path/to/output/longest_transcripts.fa"
#' )
#' }
#'
#'
#' @importFrom fastplyr f_group_by f_arrange f_slice_head
#'
#' @export
get_longest_transcript <- function(gtf_file = NULL,
                                   genome_file = NULL,
                                   output_file = NULL){
  # load gtf
  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    gtf <- rtracklayer::import.gff(gtf_file,format = "gtf")
  } else {
    warning("Package 'rtracklayer' is needed for this function to work.")
  }


  # transcriptome features
  features <- prepareTransInfo(gtf_file = gtf_file)

  # select longest trans
  features.ft <- features %>%
    fastplyr::f_group_by(gene) %>%
    fastplyr::f_arrange(cds, translen,.by_group = T,.descending = T) %>%
    fastplyr::f_slice_head(n = 1,keep_order = T)

  gtf.ft <- gtf[gtf$transcript_id %in% features.ft$transcript_id]

  # extract sequence
  get_transcript_sequence(genome_file = genome_file,
                          gtf_file = gtf.ft,
                          output_file = output_file)

  message("Extraction has been done!")

  return(features.ft)
}





#' Theme for plot
#' @export
theme_ribo <- function(){
  theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.text = element_text(face = "bold",size = rel(1))
    )
}



# boot_ci <- function(x, n = 1000, conf = 0.95){
#   stat_fun <- function(data, indices) {
#     mean(data[indices])
#   }
#
#   if (requireNamespace("boot", quietly = TRUE)) {
#     boot_obj <- boot::boot(data = x, statistic = stat_fun, R = n)
#     bci <- boot::boot.ci(boot_obj,conf = conf, type = "perc")
#   } else {
#     warning("Package 'boot' is needed for this function to work.")
#   }
#
#
#   ci_low <- bci$percent[4]
#   ci_high <- bci$percent[5]
#
#   return(c(ci_low, ci_high))
# }




#' Bootstrap confidence intervals for per-column statistics
#'
#' Calculate bootstrap confidence intervals for position-wise summary statistics
#'
#' @param mat Numeric matrix of dimension `n_samples` × `n_positions`.
#'   Each row is one independent sample (e.g. a gene), each column a position.
#' @param boot_n Integer - number of bootstrap replicates. Defaults to 1000.
#' @param method Character; one of `"median"`, `"mean"` or `"sum"`.
#'   Defaults to `"median"`.
#' @param conf Numeric in (0,1); confidence level. E.g. `0.95` for 95% CIs.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{rel}{Position identifier from column names}
#'   \item{ci_low}{Lower confidence interval bound}
#'   \item{ci_high}{Upper confidence interval bound}
#' }
#'
#' @details
#' The bootstrap procedure works as follows:
#' 1. Generate `boot_n` bootstrap samples by resampling rows with replacement
#' 2. Calculate column-wise statistics for each bootstrap sample
#' 3. Compute quantiles from bootstrap distribution
#'
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(1000), nrow = 100, dimnames = list(NULL, 1:10))
#' do_boot(mat, boot_n = 1000, method = "median")
#' }
#'
#' @export
#' @importFrom stats quantile
#' @import Rcpp
#' @import RcppArmadillo
#'
#' @useDynLib riboTransVis, .registration = TRUE
do_boot <- function(mat = NULL, boot_n = 1000,
                    method = c("median","mean","sum"), conf = 0.95) {
  # Select the requested summary statistic method
  method <- match.arg(method,choices = c("median","mean","sum"))

  # Number of original samples (rows) and positions (columns)
  n_samples <- nrow(mat)
  n_positions <- ncol(mat)

  # Create a vector of length boot_n * n_samples of resampled row-indices
  boot_indices <- matrix(sample(n_samples, boot_n * n_samples, replace = TRUE), nrow = n_samples)
  boot_indices <- matrix(as.integer(boot_indices - 1), nrow = n_samples)

  mat <- as.matrix(mat)

  # call C++:
  boot_results <- boot_stat(mat, boot_indices, method)

  # compute CI
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  ci <- t(apply(boot_results, 2, function(x) {
    quantile(x, probs = probs, na.rm = TRUE)
  }))


  # return a data.frame
  result <- data.frame(
    rel = as.numeric(colnames(mat)),
    ci_low = ci[,1],
    ci_high = ci[,2],
    stringsAsFactors = FALSE
  )

  return(result)
}




#' Generate a Codon to Amino Acid Translation Table
#'
#' This function creates a codon-to-amino acid translation table based on the standard genetic code.
#'
#' @return A data frame containing:
#'   \item{Codon}{The three-letter codon sequence.}
#'   \item{AminoAcid}{The full name of the corresponding amino acid.}
#'   \item{Abbreviation3}{The three-letter abbreviation of the amino acid.}
#'   \item{Abbreviation1}{The one-letter abbreviation of the amino acid.}
#'
#'
#' @importFrom Biostrings GENETIC_CODE
#' @export
get_aa_table <- function(){
  aa_info <- data.frame(
    Abbreviation1 = c("A", "R", "N", "D", "C", "Q", "E", "G", "H",
                      "I", "L", "K", "M", "F", "P", "S", "T", "W",
                      "Y", "V", "*"),
    Abbreviation3 = c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu",
                      "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
                      "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Stp"),
    FullName = c("Alanine", "Arginine", "Asparagine", "Aspartic acid",
                 "Cysteine", "Glutamine", "Glutamic acid", "Glycine",
                 "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",
                 "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan",
                 "Tyrosine", "Valine", "Stop codon")
  )


  codon_table <- data.frame(Codon = names(Biostrings::GENETIC_CODE),
                            Abbreviation1 = as.character(Biostrings::GENETIC_CODE))

  # combine
  final_table <- merge(codon_table, aa_info, by = "Abbreviation1", all.x = TRUE)

  # order
  final_table <- final_table[, c("Codon", "FullName", "Abbreviation3", "Abbreviation1")]
  colnames(final_table) <- c("Codon", "AminoAcid", "Abbreviation3", "Abbreviation1")

  # deal with stop codon
  final_table$AminoAcid[final_table$Abbreviation1 == "*"] <- "Stop codon"
  final_table$Abbreviation3[final_table$Abbreviation1 == "*"] <- "Stp"

  return(final_table)
}




#' Perform Offset Correction for Ribosome Footprints
#'
#' This function corrects **ribosome footprint positions** based on `reads_offset_info`,
#' adjusting **relative positions** of mapped reads.
#'
#' @param object A `ribotrans` object.
#' @param shift Integer defining the **offset shift** applied to read positions. **Default**: `0`.
#'
#' @return A modified `data.frame` with updated **ribosome footprint positions (`pos`)**.
#'
#' @details
#' - If **reads_offset_info is missing**, an error message is displayed.
#' - The offset correction is performed by adjusting `pos = pos - rel_pos + shift`.
#'
#' @examples
#' \dontrun{
#' # Apply offset correction
#' corrected_data <- do_offset_correction(obj, shift = 3)
#' head(corrected_data)
#' }
#'
#'
#' @export
do_offset_correction <- function(object = NULL, shift = 0) {
  sry <- object@summary_info

  if (nrow(object@reads_offset_info) > 0) {
    if (all(colnames(object@reads_offset_info) %in% c("sample", "qwidth", "rel_pos"))) {
      sry_offset <- object@reads_offset_info %>%
        fastplyr::f_inner_join(y = sry, by = c("sample", "qwidth")) %>%
        dplyr::mutate(pos = pos - rel_pos + shift)

      return(sry_offset)
    }
  } else {
    stop("Error: reads_offset_info is missing in the ribotrans object.")
  }
}





#' Generate Unique k-mers from Sequence Files
#'
#' This function extracts unique k-mers from a given set of sequences, which can be amino acid (proteins), DNA, or RNA.
#' It supports reading sequences from FASTA files or directly from in-memory Biostrings objects.
#'
#' - For DNA/RNA sequences, the function optionally translates them to amino acids before k-mer extraction (using the standard genetic code).
#' - Sequences containing ambiguous characters (e.g., N in DNA or RNA) will be removed by default.
#'
#' @param fa_file Character. Path to the input FASTA file. Required if \code{fa_obj} is not provided.
#' @param fa_obj Biostrings object. An optional preloaded AAStringSet, DNAStringSet, or RNAStringSet object. If provided, \code{fa_file} will be ignored.
#' @param fa_type Character. Type of input sequence. Must be one of \code{c("aa", "dna", "rna")}. Default is \code{"aa"}.
#' @param translate Logical. Whether to translate nucleotide sequences to amino acids (only applicable if \code{fa_type} is \code{"dna"} or \code{"rna"}). Default is \code{FALSE}.
#' @param kmer_length Integer. Length of k-mers to extract. Default is \code{15}.
#'
#' @return A character vector containing unique k-mers.
#'
#' @examples
#' \dontrun{
#' # Extract AA 15-mers from protein FASTA
#' kmers <- generate_kmers(fa_file = "proteome.fasta", fa_type = "aa", kmer_length = 15)
#'
#' # Extract translated 15-mer peptides from DNA
#' peptides <- generate_kmers(fa_file = "transcripts.fasta", fa_type = "dna",
#' translate = TRUE, kmer_length = 15)
#'
#' # Providing in-memory AAStringSet
#' library(Biostrings)
#' fa <- readAAStringSet("proteome.fasta")
#' kmers <- generate_kmers(fa_obj = fa, fa_type = "aa", kmer_length = 15)
#' }
#'
#'
#' @export
generate_kmers <- function(fa_file = NULL,
                           fa_obj = NULL,
                           fa_type = c("aa", "dna", "rna"),
                           translate = FALSE,
                           kmer_length = 15){
  fa_type <- match.arg(fa_type,choices = c("aa", "dna", "rna"))
  # ============================================================================
  # check data type
  if(is.null(fa_obj)){
    if(fa_type == "aa"){
      proteome <- Biostrings::readAAStringSet(fa_file)
    }else if(fa_type %in% c("dna","rna")){
      if(fa_type == "dna"){
        xna <- Biostrings::readDNAStringSet(fa_file)

        # remove transcript contains Ns in sequence
        valid_xna <- xna[!grepl("[^ACGTacgt]", as.character(xna))]
      }else{
        xna <- Biostrings::readRNAStringSet(fa_file)
        # remove transcript contains Ns in sequence
        valid_xna <- xna[!grepl("[^ACGUacgu]", as.character(xna))]
      }


      # translate codon to amino acid
      if(translate == TRUE){
        proteome <- Biostrings::translate(x = valid_xna, genetic.code = Biostrings::GENETIC_CODE)
      }else{
        proteome <- valid_xna
      }
    }
  }else{
    proteome <- fa_obj
  }

  # filter
  valid_idx <- which(Biostrings::width(proteome) >= kmer_length)
  filtered_proteome <- proteome[valid_idx]

  # ============================================================================
  # define func
  extract_nmers <- function(seq, n = 15) {
    len <- Biostrings::width(seq)

    # check length
    if(len == n){
      return(as.character(seq))
    }else{
      interval <- seq(1, len - n, by = 1)
      # get peptide
      if (requireNamespace("stringr", quietly = TRUE)) {
        stringr::str_sub_all(seq,start = interval, end = interval + (n - 1))[[1]]
      } else {
        warning("Package 'stringr' is needed for this function to work.")
      }
    }
  }

  # get all kmers
  all_nmers <- unlist(lapply(seq_along(filtered_proteome), function(x){
    extract_nmers(seq = filtered_proteome[x], n = kmer_length)
  }))

  # remove duplicate kmers
  all_nmers_unique <- unique(all_nmers)

  return(all_nmers_unique)
}





#' Convert Aligned Sequences to a Position Frequency Matrix
#'
#' This function converts a set of aligned biological sequences into a position frequency matrix (PFM) using Biostrings::consensusMatrix(). It supports amino acid (AA), DNA, and RNA sequence types and filters the result to include only expected characters for the given type.
#'
#' @param seqs A Biostrings XStringSet object. Must be an instance of AAStringSet (for protein sequences), DNAStringSet (for DNA), or RNAStringSet (for RNA).
#' @param seq_type Character. Type of sequence, must be one of: \code{"aa"}, \code{"dna"}, or \code{"rna"}. Determines the allowed residues/bases for row filtering.
#'
#' @return A matrix where rows correspond to residues or bases, and columns correspond to positions in the aligned sequences. Each value in the matrix represents the count of a particular character at a given position.
#'
#' @details
#' This function uses \code{Biostrings::consensusMatrix()} to calculate a position-wise summary of aligned sequences. The result is filtered to retain only rows corresponding to valid residues/bases for the specified sequence type:
#' \itemize{
#' \item \code{"aa"}: 20 standard amino acids (A, R, N, D, ... V)
#' \item \code{"dna"}: A, C, G, T
#' \item \code{"rna"}: A, C, G, U
#' }
#'
#' @seealso \code{\link[Biostrings]{consensusMatrix}}
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#'
#' # DNA example
#' dna_seqs <- DNAStringSet(c("ATGG", "ATGC", "ATGT"))
#' mat_dna <- to_position_matrix(dna_seqs, seq_type = "dna")
#' print(mat_dna)
#'
#' # Protein example
#' aa_seqs <- AAStringSet(c("ARND", "ARNE", "ARNQ"))
#' mat_aa <- to_position_matrix(aa_seqs, seq_type = "aa")
#' print(mat_aa)
#' }
#'
#' @importFrom Biostrings consensusMatrix
#' @export
to_position_matrix <- function(seqs,seq_type){
  if(seq_type == "aa"){
    base_char <- c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  }else if(seq_type == "dna"){
    base_char <- c("A", "G", "C", "T")
  }else if(seq_type == "rna"){
    base_char <- c("A", "G", "C", "U")
  }

  mat <- Biostrings::consensusMatrix(seqs, as.prob = FALSE)
  mat <- mat[intersect(base_char, rownames(mat)), ]
  return(mat)
}



#' Calculate Log-Ratio Matrix Based on Binomial Significance
#'
#' This function computes a log-ratio matrix for each amino acid at each position
#' based on a foreground count matrix and corresponding background frequencies.
#' The log-ratio quantifies the statistical significance of over- or under-representation
#' of each residue using a binomial model, and results are suitable for pLogo-style visualization.
#'
#' @param fg A numeric matrix representing the counts of residues in the foreground data.
#'           Rows correspond to residue types (e.g., amino acids), and columns correspond to positions.
#'           Row names should be residue single-letter codes.
#' @param bg (Optional) Background matrix. This argument is unused in the current function.
#' @param bg_freq A numeric matrix of same dimension as \code{fg}, representing the expected
#'        background frequencies (probabilities) for each residue at each position.
#'
#' @return A numeric matrix of the same dimensions as \code{fg}, where each element gives
#'         the negative base-10 logarithm of the over/under binomial probability ratio:
#'         \deqn{-log_{10}\left( \frac{P(k \geq K \mid N, p)}{P(k \leq K \mid N, p)} + \varepsilon \right)}
#'         If the observed count \code{K} is zero, the corresponding log-ratio is set to 0.
#'
#' @details
#' For each residue \code{i} and position \code{j}, this function calculates:
#' \itemize{
#'   \item \code{K}: Observed count of residue i at position j in the foreground.
#'   \item \code{N}: Total count of residues at position j in the foreground.
#'   \item \code{p}: Background frequency for residue i at position j.
#' }
#' Then computes:
#' \itemize{
#'   \item \code{pr_over = P(k >= K | N, p)}: Binomial upper-tail probability
#'   \item \code{pr_under = P(k <= K | N, p)}: Binomial lower-tail probability
#'   \item \code{ratio = pr_over / pr_under + 1e-10}
#'   \item The final score is \code{-log10(ratio)}
#' }
#' This formulation captures the log-ratio between over- and underrepresentation probabilities,
#' following a heuristic similar to pLogo residue height calculations.
#'
#'
#' @export
get_log_ratio_mat <- function(fg = NULL,bg = NULL,bg_freq = NULL) {
  aa <- rownames(fg)
  n_pos <- ncol(fg)
  logmat <- matrix(0, nrow = length(aa), ncol = n_pos)
  rownames(logmat) <- aa

  # add pvalue
  for (j in 1:n_pos) {
    for (i in 1:length(aa)) {
      K <- fg[i, j]
      N <- sum(fg[, j])
      p <- bg_freq[i, j]

      # calculate binomial pvalue
      if (K == 0){
        logmat[i, j] <- 0
      }else{
        pr_over <- 1 - stats::pbinom(q = K - 1, size = N, prob = p)
        pr_under <- stats::pbinom(q = K, size = N, prob = p)
        ratio <- pr_over / pr_under + 1e-10

        logmat[i, j] <- -log10(ratio)
      }
    }
  }
  logmat
}
