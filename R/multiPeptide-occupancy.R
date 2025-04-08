# ==============================================================================
# method
# ==============================================================================

#' @title Compute Multi-Peptide Ribosome Occupancy Scores
#'
#' @description
#' This function calculates the ribosome density at the tri-peptide level from
#' Ribo-seq data, summarizing their pause scores and occurrences.
#'
#' @details
#' The function processes reads mapped to coding regions, filters low-coverage
#' transcripts, computes pause densities, and translates nucleotide sequences into
#' peptides. It subsequently aggregates pause scores for tri-peptides and filters
#' based on their occurrence.
#'
#' @param object A `ribotrans` object containing **Ribo-seq** data.
#' @param merge_rep Logical. Whether to merge replicate samples by \code{sample_group}. Default is \code{FALSE}.
#' @param cds_fa A **FASTA file path** containing CDS sequences for translation.
#' @param do_offset_correct Logical. If `TRUE`, performs **offset correction**
#' using `do_offset_correction()`. **Default**: `FALSE`.
#' @param position_shift Integer defining how much to adjust **ribosome footprint positions**
#' during offset correction. **Default**: `0`.
#' @param exclude_length A numeric vector of length 2 (default: `c(100,100)`)
#'   specifying the number of nucleotides to exclude from the start and end of CDS.
#' @param min_counts An integer (default: `64`), the minimum read count threshold
#'   for filtering low-expression CDS.
#' @param peptide_length An integer (default: `3`), defining the length of the
#'   peptide fragment (e.g., **tri-peptides**).
#' @param peptide_occurrence An integer (default: `100`), the minimum occurrence
#'   threshold for tri-peptide inclusion in the final dataset.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named `list` containing two data frames:
#' \itemize{
#'   \item `long_format`: A tidy-format (long) `data.frame`, listing tri-peptides,
#'      their pause scores, and occurrences across samples.
#'   \item `wider_format`: A wide-format (wide) `data.frame`, where pause scores
#'      are spread across different samples as columns.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Load a ribotrans object with ribosome footprint data
#'   ribo_obj <- load_ribotrans_data("riboseq_data.rds")
#'
#'   # Specify a FASTA file containing CDS sequences
#'   cds_fa_path <- "cds_sequences.fasta"
#'
#'   # Compute tri-peptide ribosome occupancy
#'   result <- multi_peptide_occupancy(ribo_obj, cds_fa = cds_fa_path)
#'
#'   # Extract long and wide format results
#'   long_format <- result$long_format
#'   wide_format <- result$wider_format
#' }
#'
#' @import dplyr
#' @importFrom Biostrings readDNAStringSet translate GENETIC_CODE width
setGeneric("multi_peptide_occupancy",function(object,...) standardGeneric("multi_peptide_occupancy"))




#' @rdname multi_peptide_occupancy
#' @export
setMethod("multi_peptide_occupancy",
          signature(object = "ribotrans"),
          function(object,
                   merge_rep = FALSE,
                   cds_fa = NULL,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   exclude_length = c(100,100),
                   min_counts = 64,
                   peptide_length = 3,
                   peptide_occurrence = 100){
            # ==================================================================
            # filter coding gene and get peptide seqs
            features <- object@features %>%
              dplyr::filter(cds > 0)

            cds <- Biostrings::readDNAStringSet(cds_fa)

            # remove transcript contains Ns in sequence
            valid_cds <- cds[!grepl("[^ACGTacgt]", as.character(cds))]

            # translate codon to amino acid
            aa <- Biostrings::translate(x = valid_cds,genetic.code = Biostrings::GENETIC_CODE)

            ids_retain <- intersect(features$idnew,names(aa))
            aa <- aa[ids_retain]

            # Filter sequences whose length > peptide_length
            aa <- aa[Biostrings::width(aa) > peptide_length]

            # ==================================================================
            # whether do reads offset correction
            if(do_offset_correct == TRUE){
              sry <- do_offset_correction(object = object,shift = position_shift)
            }else{
              sry <- object@summary_info
            }

            # exclude cds first and last nts
            sry <- sry %>%
              fastplyr::f_filter(mstart != 0 | mstop != 0) %>%
              fastplyr::f_select(-qwidth,-translen) %>%
              dplyr::mutate(relst = pos - mstart, relsp = pos - mstop) %>%
              fastplyr::f_filter(relst > exclude_length[1] & relsp < -exclude_length[2])

            # ==================================================================
            # expect reads per position
            avg.ct <- sry %>%
              dplyr::mutate(cdslen = mstop - mstart + 1) %>%
              fastplyr::f_group_by(sample,rname,cdslen) %>%
              fastplyr::f_summarise(counts = sum(count)) %>%
              dplyr::mutate(avg_ct = counts/cdslen) %>%
              fastplyr::f_select(sample,rname,counts,avg_ct)

            # average reads
            density.tt <- sry %>%
              dplyr::inner_join(y = avg.ct,by = c("sample", "rname")) %>%
              # filter low counts
              fastplyr::f_filter(counts > min_counts) %>%
              dplyr::mutate(norm = count/avg_ct) %>%
              fastplyr::f_group_by(sample,sample_group,rname,relst) %>%
              fastplyr::f_summarise(normsm = sum(norm)) %>%
              # codon position
              dplyr::mutate(rel = (relst %/% 3) + 1) %>%
              fastplyr::f_group_by(sample,sample_group,rname,rel) %>%
              fastplyr::f_summarise(value = mean(normsm)) %>%
              fastplyr::f_filter(rname %in% features$idnew)

            # whether aggregate replicates
            if(merge_rep == TRUE){
              density.tt <- density.tt %>%
                fastplyr::f_group_by(sample_group,rname,rel) %>%
                fastplyr::f_summarise(value = mean(value)) %>%
                dplyr::rename(sample = sample_group)
            }else{
              density.tt <- density.tt %>% dplyr::select(-sample_group)
            }
            # ==================================================================
            # count codon numbers
            anno <- dplyr::filter(features, idnew %in% names(aa)) %>%
              dplyr::mutate(codons = cds / 3 + 1)

            # make codon index
            idfull <- tidyr::uncount(anno, weights = codons) %>%
              dplyr::group_by(idnew) %>%
              dplyr::mutate(rel = dplyr::row_number()) %>%
              dplyr::ungroup() %>%
              dplyr::select(rname = idnew, rel)

            sp <- unique(density.tt$sample)

            # add samples
            purrr::map_df(seq_along(sp),function(x){
              idfull$sample <- sp[x]
              return(idfull)
            }) -> idfull

            # merge with density
            fullanno <- idfull %>%
              fastplyr::f_left_join(y = density.tt,by = c("sample","rname","rel"))

            fullanno$value[is.na(fullanno$value)] <- 0

            # ==================================================================
            # prepare tri-peptide info for all transcripts
            # x = 1
            lapply(seq_along(aa),function(x){
              tmp <- aa[x]
              cdslen <- Biostrings::width(tmp)

              # get tri-peptide
              interval <- seq(1, cdslen - (peptide_length - 1), by = 1)

              # get peptide
              seqs <- stringr::str_sub_all(tmp,start = interval,
                                           end = interval + (peptide_length - 1))[[1]]

              pepdf <- data.frame(rname = names(tmp),rel = 1:length(interval),pep_seq = seqs)

              return(pepdf)
            }) %>% do.call("rbind",.) %>% data.frame() -> peptide_info

            # ==================================================================
            # calculate total density for peptide
            if (requireNamespace("zoo", quietly = TRUE)) {
              tripep_pause_score <- fullanno %>%
                dplyr::mutate(tripep_val = zoo::rollsum(value, k = peptide_length, fill = NA),
                              .by = c(sample,rname)) %>%
                na.omit() %>%
                fastplyr::f_group_by(sample,rname) %>%
                # re-assign positions
                dplyr::mutate(rel = 1:dplyr::n()) %>%
                # merge with multi-peptide info
                fastplyr::f_inner_join(y = peptide_info,by = c("rname","rel")) %>%
                fastplyr::f_group_by(sample,pep_seq) %>%
                fastplyr::f_summarise(pause_score = sum(tripep_val),
                                      occurrence = dplyr::n()) %>%
                dplyr::mutate(pause_score = ifelse(pause_score < 1e-10, 0, pause_score)) %>%
                fastplyr::f_filter(occurrence > peptide_occurrence) %>%
                # calculate relative density
                dplyr::mutate(rel_pause = pause_score/occurrence)

            } else {
              warning("Package 'zoo' is needed for this function to work.")
            }


            # long format to wider format
            if (requireNamespace("tidyr", quietly = TRUE)) {
              tripep_pause_score_wide <- tripep_pause_score %>%
                dplyr::select(sample,pep_seq,rel_pause) %>%
                tidyr::pivot_wider(names_from = sample, values_from = rel_pause)
            } else {
              warning("Package 'tidyr' is needed for this function to work.")
            }


            # return
            res <- list(long_format = tripep_pause_score,
                        wider_format = tripep_pause_score_wide)
          }
)





#' @title Scatter Plot of Peptide Pause Scores
#'
#' @description
#' This function generates a scatter plot comparing ribosome pause scores
#' between two sample conditions, highlighting the **top N peptides** with
#' the highest ratio difference.
#'
#' @details
#' The function plots each peptide's ribosome occupancy in two conditions
#' (specified as `x` and `y`). It highlights the `top_motif` peptides with
#' the **largest ratio difference (`y / x`)**, labeling them on the plot.
#' A **dashed y = x line** is included as a reference.
#'
#' @param data A `data.frame` containing peptide sequences and their pause scores
#'   in different sample conditions.
#' @param x A character string specifying the column name for the **X-axis condition**.
#' @param y A character string specifying the column name for the **Y-axis condition**.
#' @param color A string defining the **scatter color** (default: `"#003366"`).
#' @param label_size A numeric value controlling the labeled peptide text size
#'   (default: `3`).
#' @param hjust A numeric value specifying horizontal justification for the labels. The default is `1.2`.
#' @param vjust A numeric value specifying vertical justification for the labels. The default is `1.2`.
#' @param top_motif An integer specifying the number of top peptides to highlight based
#'   on **highest fold-change (`y / x`)** (default: `20`).
#'
#' @return A **ggplot2 scatter plot** comparing `x` vs. `y` peptide density values.
#'   The top `top_motif` peptides with **largest fold-change** are labeled.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Generate scatter plot comparing wt_rep1 vs sgeIF5A_rep1
#'   peptide_scatter_plot(
#'     data = tripep_pause_score_wide,
#'     x = "wt_rep1",
#'     y = "sgeIF5A_rep1",
#'     color = "#009933",
#'     top_motif = 5
#'   )
#' }
#'
peptide_scatter_plot <- function(data = NULL,
                                 x = NULL,y = NULL,
                                 color = "#003366",
                                 label_size = 3,
                                 hjust = 1.2,vjust = 1.2,
                                 top_motif = 20){

  # get ratio
  label_df <- data[,c("pep_seq",x,y)]
  label_df$ratio <- label_df[,3]/label_df[,2]
  label_df <- label_df %>% dplyr::arrange(dplyr::desc(ratio)) %>%
    dplyr::slice_max(order_by = ratio,n = top_motif)

  if (requireNamespace("grDevices", quietly = TRUE)) {
    xlims <- grDevices::extendrange(data[,2:ncol(data)],f = 0.1)[2]
  } else {
    warning("Package 'grDevices' is needed for this function to work.")
  }


  # plot
  ggplot(data) +
    geom_point(aes(x = get(x),y = get(y)),color = color) +
    geom_text(data = label_df,
              aes(x = get(x), y = get(y), label = pep_seq),
              hjust = hjust, vjust = vjust, size = label_size, check_overlap = TRUE) +
    geom_abline(slope = 1,intercept = 0,lty = "dashed", color = "grey40") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black")) +
    coord_equal() +
    xlab(x) + ylab(y) +
    xlim(0, xlims) +
    ylim(0, xlims)
}
