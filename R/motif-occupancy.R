
#' Motif Occupancy Analysis
#'
#' This function calculates the occupancy of specified motifs (either amino acid string or codon pattern) within coding sequence (CDS) regions using ribosome profiling data stored in a ribotrans object.
#' It optionally performs read offset correction, aggregates biological replicates, and visualizes the relative occupancy as cumulative distribution curves (ECDFs).
#'
#' @param object A \code{ribotrans} object containing processed ribosome profiling data.
#' @param merge_rep Logical. Whether to merge replicates (\code{sample_group}) during analysis. Default: \code{FALSE}.
#' @param cds_fa Character. Path to the CDS FASTA file. The sequences are filtered and translated depending on \code{search_type}.
#' @param do_offset_correct Logical. Whether to perform read offset correction before summarizing coverage. Default: \code{FALSE}.
#' @param position_shift Numeric. Integer value to apply as position shift during offset correction. Default: \code{0}.
#' @param search_type Character. Specifies motif space used for searching: either \code{"amino"} (peptide sequence) or \code{"codon"} (nucleotide codon sequences). Default: \code{"amino"}.
#' @param exclude_length A numeric vector of length two specifying the number of nucleotides
#' to exclude from the start and end of the CDS to avoid boundary effects. Default is `c(100, 100)`.
#' @param min_counts Integer. Minimum number of total reads required for a transcript to be included in occupancy calculation. Default: 64.
#' @param motif_pattern Character vector. Motif string(s) to search in amino acid or codon space, e.g., \code{c("PPP", "KKK")}.
#' @param facet_layer ggplot2 facet call. Faceting layer for ECDF plot. Default: \code{facet_wrap(~motif)}.
#' @param return_data Logical. If \code{TRUE}, returns the data.frame of motif-level occupancy; otherwise returns a ggplot object. Default: \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return If \code{return_data = FALSE} (default), returns a \code{ggplot} object showing ECDF curves of log2-transformed codon/amino acid occupancies at motif locations.
#' If \code{return_data = TRUE}, returns a \code{data.frame} containing normalized occupancy values for each motif instance.
#'
#' @details This function performs:
#' \enumerate{
#' \item CDS filtering and translation
#' \item Identification of motif positions across transcripts
#' \item Read normalization (per-length scaling)
#' \item Relative density calculation
#' \item Visualization of cumulative codon/amino acid occupancy at motif positions
#' }
#'
#' @import ggplot2
#' @importFrom Biostrings readDNAStringSet translate vmatchPattern start
#' @importFrom dplyr filter mutate inner_join rename select summarise group_by
#' @importFrom fastplyr f_filter f_group_by f_select f_summarise
#'
#' @examples
#' \dontrun{
#' motif_occupancy(
#' object = ribo_obj,
#' cds_fa = "coding_sequences.fasta",
#' search_type = "amino",
#' motif_pattern = c("PPP", "KKK"),
#' do_offset_correct = TRUE,
#' merge_rep = TRUE
#' )
#' }
#'
#' @export
#' @seealso \code{\link{do_offset_correction}}
setGeneric("motif_occupancy",function(object,...) standardGeneric("motif_occupancy"))





#' @rdname motif_occupancy
#' @export
setMethod("motif_occupancy",
          signature(object = "ribotrans"),
          function(object,
                   merge_rep = FALSE,
                   cds_fa = NULL,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   search_type = c("amino", "codon"),
                   exclude_length = c(100,100),
                   min_counts = 64,
                   motif_pattern = "PPP",
                   facet_layer = facet_wrap(~motif),
                   return_data = FALSE){
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
                motif.pos <- motif.pos %>% fastplyr::f_filter(pos %% 3 == 0) %>%
                  dplyr::mutate(pos = pos %/% 3)
              }

              return(motif.pos)
            }) -> motif.pos

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
              fastplyr::f_group_by(sample,sample_group,rname,relst) %>%
              fastplyr::f_summarise(normsm = sum(norm)) %>%
              # codon position
              dplyr::mutate(codon_pos = (relst %/% 3) + 1) %>%
              fastplyr::f_group_by(sample,sample_group,rname,codon_pos) %>%
              fastplyr::f_summarise(value = mean(normsm))

            # whether aggregate replicates
            if(merge_rep == TRUE){
              density.tt <- density.tt %>%
                fastplyr::f_group_by(sample_group,rname,codon_pos) %>%
                fastplyr::f_summarise(value = mean(value)) %>%
                dplyr::rename(sample = sample_group)
            }

            # ==================================================================
            # calculate relative occupancy
            rel2motif.df <- density.tt %>%
              dplyr::inner_join(y = motif.pos,by = "rname",relationship = "many-to-many") %>%
              dplyr::filter(codon_pos == pos)

            # ==================================================================
            # plot
            p <-
              ggplot(rel2motif.df) +
              stat_ecdf(mapping = aes(x = log2(value),color = sample)) +
              theme_bw(base_size = 12) +
              facet_layer +
              theme(panel.grid = element_blank(),
                    # strip.background = element_blank(),
                    strip.text = element_text(face = "bold",size = rel(1)),
                    axis.text = element_text(colour = "black")) +
              ylab("Cumulative fracton") +
              xlab("Log2 (codon/amino occupancy)")


            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(rel2motif.df)
            }
          }
)
