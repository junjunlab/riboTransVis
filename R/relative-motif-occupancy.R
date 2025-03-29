# ==============================================================================
# method
# ==============================================================================

#' Calculate Relative Motif Occupancy in Ribo-seq Data
#'
#' This function calculates the relative ribosome occupancy around specific
#' motifs in coding sequences (CDS) using Ribo-seq data.
#'
#' @param object An object of class `ribotrans` that contains experimental data, including
#' ribosome profiling information.
#' @param cds_fa A character string specifying the file path to a FASTA file containing
#' the nucleotide sequences of coding regions (CDS).
#' @param do_offset_correct Logical. If `TRUE`, performs **offset correction**
#' using `do_offset_correction()`. **Default**: `FALSE`.
#' @param position_shift Integer defining how much to adjust **ribosome footprint positions**
#' during offset correction. **Default**: `0`.
#' @param search_type Search mode, either `"amino"` (default) to search **protein motifs**
#' or `"codon"` to search **specific codon triplets**.
#' @param exclude_length A numeric vector of length two specifying the number of nucleotides
#' to exclude from the start and end of the CDS to avoid boundary effects. Default is `c(100, 100)`.
#' @param min_counts A numeric value indicating the minimum read counts required for a transcript
#' to be included in the analysis. Default is `64`.
#' @param motif_pattern A character vector specifying one or more amino acid motifs to be analyzed.
#' @param motif_upstream An integer specifying the number of nucleotides upstream of the motif
#' to consider for analysis. Default is `50`.
#' @param motif_downstream An integer specifying the number of nucleotides downstream of the motif
#' to consider for analysis. Default is `50`.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `data.frame` containing the relative occupancy of ribosomes around each motif.
#' The columns include:
#' \itemize{
#'   \item `sample` - The sample name.
#'   \item `motif` - The identified motif.
#'   \item `dist` - The relative nucleotide position from the motif.
#'   \item `value` - The unnormalized count of ribosome footprints.
#'   \item `avg_val` - The normalized translation occupancy (relative to transcript averages).
#' }
#'
#'
#' @details
#' - This function first **filters and extracts coding sequences from input FASTA (`cds_fa`)**
#' and translates sequences into **protein sequences** using `Biostrings::translate()`.
#' - If `search_type = "amino"`, the function searches **protein motifs** (e.g., `"PPP"` for Proline triplets).
#' - If `search_type = "codon"`, the function searches **specific codons** (e.g., `"CCG"` for Proline in codon format).
#' - Only **motifs occurring in the correct reading frame** are counted.
#' - **Ribosome occupancy** is calculated in a **fixed window (`motif_upstream` & `motif_downstream`)**
#' around the motif center.
#' - Transcripts with **low read counts (< `min_counts`) are excluded** to reduce noise.
#'
#' @examples
#' \dontrun{
#' result <- relative_motif_occupancy(
#'   object = ribo_data,
#'   search_type = "amino",
#'   cds_fa = "cds_sequences.fasta",
#'   motif_pattern = c("PPP", "PPAP"),
#'   motif_upstream = 100,
#'   motif_downstream = 100,
#'   min_counts = 50
#' )
#' head(result)
#' }
#'
#' @importFrom Biostrings readDNAStringSet translate vmatchPattern
#' @importFrom dplyr filter mutate select group_by summarise left_join
#' @importFrom fastplyr f_filter f_select f_group_by f_summarise
#' @export
setGeneric("relative_motif_occupancy",function(object,...) standardGeneric("relative_motif_occupancy"))





#' @rdname relative_motif_occupancy
#' @export
setMethod("relative_motif_occupancy",
          signature(object = "ribotrans"),
          function(object,
                   cds_fa = NULL,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   search_type = c("amino", "codon"),
                   exclude_length = c(100,100),
                   min_counts = 64,
                   motif_pattern = "PPP",
                   motif_upstream = 50,
                   motif_downstream = 50){
            search_type <- match.arg(search_type,choices = c("amino", "codon"))
            # ==================================================================
            # filter coding gene and get peptide seqs
            features <- object@features %>% dplyr::filter(cds > 0)

            cds <- Biostrings::readDNAStringSet(cds_fa)

            # remove transcript contains Ns in sequence
            valid_cds <- cds[!grepl("[^ACGTacgt]", as.character(cds))]

            # translate codon to amino acid
            aa <- Biostrings::translate(x = valid_cds,genetic.code = Biostrings::GENETIC_CODE)

            ids_retain <- intersect(features$idnew,names(aa))
            aa <- aa[ids_retain]

            # ==================================================================
            # find specfic motif on each transcript

            # loop for each motif
            purrr::map_df(seq_along(motif_pattern),function(x){
              # check search_type
              if(search_type == "codon"){
                motif_matches <- Biostrings::vmatchPattern(motif_pattern[x], cds)
              }else{
                motif_matches <- Biostrings::vmatchPattern(motif_pattern[x], aa)
              }


              # to dataframe
              if (requireNamespace("data.table", quietly = TRUE)) {
                motif.pos <- data.table::rbindlist(lapply(seq_along(motif_matches), function(i) {
                  pos <- Biostrings::start(motif_matches[[i]])
                  if (length(pos) > 0) {
                    data.table::data.table(motif = motif_pattern[x], rname = names(motif_matches)[i], pos = pos)
                  }
                }), fill = TRUE)
              } else {
                warning("Package 'data.table' is needed for this function to work.")
              }

              # filter in-frame codon position
              if(search_type == "codon"){
                motif.pos <- motif.pos %>% fastplyr::f_filter(pos %% 3 == 1)
              }

              return(motif.pos)
            }) -> motif.pos

            # condon position to nucleotide position
            if(search_type == "amino"){
              motif.pos <- motif.pos %>% dplyr::mutate(nt_pos = pos*3 - 2)
            }else{
              motif.pos <- motif.pos %>% dplyr::mutate(nt_pos = pos)
            }


            # ==================================================================
            # load density

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
              fastplyr::f_group_by(sample,rname,relst) %>%
              fastplyr::f_summarise(normsm = sum(norm))

            # filter transcript in motifs found
            density.tt <- density.tt %>%
              fastplyr::f_filter(rname %in% unique(motif.pos$rname))

            # ==================================================================
            # calculate relative occupancy
            tid <- unique(density.tt$rname)

            rel2motif.df <- density.tt %>%
              dplyr::inner_join(y = motif.pos,by = "rname",relationship = "many-to-many") %>%
              dplyr::mutate(dist = relst - nt_pos) %>%
              fastplyr::f_filter(dist >= -motif_upstream & dist <= motif_downstream)


            # average occupancy in fixed region for each sample
            sp <- unique(rel2motif.df$sample)

            purrr::map_df(seq_along(sp),function(x){
              tmp <- subset(rel2motif.df, sample == sp[x]) %>%
                fastplyr::f_group_by(sample,motif,dist) %>%
                fastplyr::f_summarise(value = sum(normsm))

              avg.exp <- sum(tmp$value)/nrow(tmp)
              tmp$avg_val <- tmp$value/avg.exp

              return(tmp)
            }) -> rel2motif.avg.df

            # final results
            return(rel2motif.avg.df)
          }
)




#' Plot Relative Motif Occupancy from Ribo-seq Data
#'
#' This function generates a plot of ribosome occupancy relative to specific
#' amino acid motifs identified in coding sequences (CDS) using Ribo-seq data.
#'
#' @param data A `data.frame` containing the processed ribosome occupancy data.
#' It must include the following columns:
#' \itemize{
#'   \item `dist` - The relative nucleotide distance from the motif.
#'   \item `avg_val` - The average normalized ribosome occupancy.
#'   \item `sample` - The sample name.
#'   \item `motif` - The amino acid motif.
#' }
#'
#' @param mode Character string specifying the unit for the x-axis:
#'   - `"nt"`: Nucleotide level resolution (default).
#'   - `"codon"`: Codon-level resolution.
#'
#' @return A `ggplot2` object visualizing ribosome occupancy relative to motifs.
#' The plot includes:
#' \itemize{
#'   \item A **line plot** (`geom_path`) of ribosome occupancy (`avg_val`) across nucleotide positions (`dist`).
#'   \item A **facet grid** (`facet_grid(motif ~ sample)`) to display results separately for each motif/sample.
#'   \item Custom **theme** settings for improved visualization.
#' }
#'
#' @details
#' This function takes as input the output from `relative_motif_occupancy()`, which
#' contains motif-centered ribosome occupancy profiles, and visualizes them in a grid format
#' where each row represents a different motif and each column represents a different sample.
#'
#' @examples
#' \dontrun{
#' result <- relative_motif_occupancy(object = ribo_data, cds_fa = "cds_sequences.fasta",
#'                                    motif_pattern = c("PPP", "PPAP"), motif_upstream = 100,
#'                                    motif_downstream = 100, min_counts = 50)
#'
#' # Generate the occupancy plot
#' relative_motif_plot(result)
#' }
#'
#
#' @export
relative_motif_plot <- function(data = NULL, mode = c("nt", "codon")){
  mode <- match.arg(mode,choices = c("nt", "codon"))

  data <- data %>% dplyr::mutate(codon_dist = dist %/% 3)

  # check mode
  if(mode == "codon"){
    xaes <- "codon_dist"
  }else{
    xaes <- "dist"
  }

  # plot
  ggplot(data) +
    geom_path(aes(x = get(xaes),y = avg_val)) +
    # facet_wrap(~sample) +
    facet_grid(motif~sample,scales = "free_y") +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text = element_text(colour = "black")) +
    xlab(paste("Distance of E/P/A site to motif"," (",mode,")",sep = "")) +
    ylab("Average ribosome occupancy")

}
