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
