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





#' Generate Unique k-mer Peptides from Coding or Peptide Sequences
#'
#' This function extracts all unique amino acid k-mers (subsequences of fixed length) from
#' either DNA coding sequences (CDS) or peptide sequences. It supports automatic translation
#' from CDS to protein sequences using the standard genetic code and filters out invalid input
#' such as sequences with ambiguous nucleotides.
#'
#' @param cds_fa Character. File path to a FASTA file containing CDS (DNA sequences). If this is
#'     provided, the sequences will be translated to peptides using \code{Biostrings::translate()}.
#' @param peptide_fa Character. File path to a FASTA file containing amino acid sequences. Used
#'     only if \code{cds_fa} is NULL.
#' @param kmer_length Integer. Length of the peptide k-mers to be extracted. Default is 15.
#'
#' @return A character vector containing all unique amino acid k-mer subsequences
#'     (length = \code{kmer_length}) from the translated or provided protein sequences.
#'
#' @details
#' This function performs the following operations:
#' \itemize{
#'   \item Reads CDS (DNA) or peptide sequences from a FASTA file.
#'   \item Filters out DNA sequences containing ambiguous nucleotides (non-ACGT).
#'   \item Translates CDS into amino acid sequences using the standard genetic code.
#'   \item Selects sequences with length >= \code{kmer_length}.
#'   \item Extracts overlapping k-mers from each sequence using a sliding window (step = 1).
#'   \item Removes duplicate k-mers and returns a unique list.
#' }
#' Requires the \pkg{Biostrings} and \pkg{stringr} packages.
#'
#' @importFrom Biostrings readDNAStringSet translate readAAStringSet width
#'
#' @examples
#' \dontrun{
#' # Example: Generate 15-mers from CDS
#' kmers <- generate_kmers(cds_fa = "transcripts.fasta", kmer_length = 15)
#' head(kmers)
#'
#' # Example: Generate 10-mers from amino acid sequences
#' kmers <- generate_kmers(peptide_fa = "proteins.fasta", kmer_length = 10)
#' head(kmers)
#' }
#'
#' @export
generate_kmers <- function(cds_fa = NULL,
                           peptide_fa = NULL,
                           kmer_length = 15){
  # check data type
  if(!is.null(cds_fa)){
    cds <- Biostrings::readDNAStringSet(cds_fa)

    # remove transcript contains Ns in sequence
    valid_cds <- cds[!grepl("[^ACGTacgt]", as.character(cds))]

    # translate codon to amino acid
    proteome <- Biostrings::translate(x = valid_cds,genetic.code = Biostrings::GENETIC_CODE)

  }else{
    proteome <- Biostrings::readAAStringSet(peptide_fa)
  }

  # filter
  valid_idx <- which(Biostrings::width(proteome) >= kmer_length)
  filtered_proteome <- proteome[valid_idx]

  # define func
  extract_nmers <- function(seq, n = 15) {
    len <- nchar(seq)
    # sapply(1:(len - n), function(i) substr(seq, i, i + n - 1))
    interval <- seq(1, len - n, by = 1)
    # get peptide
    if (requireNamespace("stringr", quietly = TRUE)) {
      stringr::str_sub_all(seq,start = interval, end = interval + (n - 1))[[1]]
    } else {
      warning("Package 'stringr' is needed for this function to work.")
    }

  }

  # get all kmers
  all_nmers <- unlist(lapply(filtered_proteome, function(x) extract_nmers(x, n = kmer_length)))

  # remove duplicate kmers
  all_nmers_unique <- unique(all_nmers)

  return(all_nmers_unique)
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
