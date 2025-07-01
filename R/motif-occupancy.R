
#' Calculate motif occupancy from ribosome profiling data
#'
#' @description
#' This method calculates the occupancy of specific amino acid or codon motifs
#' from ribosome profiling (Ribo-seq) data. It normalizes ribosome footprint
#' counts to identify positions with high occupancy, which typically indicate
#' slower translation or ribosome pausing sites.
#'
#' @details
#' The method performs several key steps:
#' \enumerate{
#'   \item Optionally applies offset correction to ribosome footprint positions
#'   \item Filters transcripts based on read counts and length criteria
#'   \item Calculates normalized occupancy values relative to transcript-specific
#'         average coverage
#'   \item Searches for specified motifs in either amino acid or codon sequences
#'   \item Generates statistical comparisons and visualizations
#' }
#'
#' The occupancy calculation uses the formula:
#' \code{Normalized_Occupancy = (Position_Count / Transcript_Average_Count)}
#'
#' Higher occupancy values indicate slower ribosome movement, potentially due to:
#' \itemize{
#'   \item Rare codons with limited tRNA availability
#'   \item Secondary structure formation
#'   \item Regulatory pausing sites
#'   \item Co-translational folding events
#' }
#'
#' @param object A \code{ribotrans} object containing ribosome profiling data
#' @param merge_rep Logical. Whether to merge biological replicates by averaging
#'   their occupancy values. Default is \code{FALSE}.
#' @param cds_fa Character. Path to the CDS (Coding DNA Sequence) FASTA file
#'   required for sequence analysis. Must be provided for motif searching.
#' @param do_offset_correct Logical. Whether to apply P-site offset correction
#'   to account for the distance between the 5' end of ribosome footprints and
#'   the ribosome P-site. Default is \code{FALSE}.
#' @param position_shift Numeric. Additional position shift to apply when
#'   \code{do_offset_correct = TRUE}. Default is \code{0}.
#' @param search_type Character. Type of motif search to perform:
#'   \itemize{
#'     \item \code{"amino"}: Search for amino acid motifs in translated sequences
#'     \item \code{"codon"}: Search for codon motifs in DNA sequences
#'   }
#'   Default is \code{"amino"}.
#' @param exclude_length Numeric vector of length 2. Number of codons to exclude
#'   from the start and end of CDS regions \code{c(start_exclude, end_exclude)}
#'   to avoid edge effects. Default is \code{c(45, 45)}.
#' @param min_counts Numeric. Minimum total read counts required per transcript
#'   for inclusion in analysis. Transcripts with fewer reads are filtered out.
#'   Default is \code{64}.
#' @param motif_pattern Character vector. The motif pattern(s) to search for.
#'   Can be amino acid codes (e.g., \code{"PPP"} for proline repeats) or codon
#'   sequences (e.g., \code{c("CCC", "CCG")} for proline codons).
#'   Default is \code{"PPP"}.
#' @param ref_group Character. Reference group for statistical comparisons.
#'   If \code{NULL}, all pairwise comparisons are performed.
#' @param insert_box Logical. Whether to insert a boxplot into the main
#'   cumulative distribution plot. Default is \code{TRUE}.
#' @param box_width Numeric. Width of the boxplot error bars.
#'   Default is \code{0.25}.
#' @param error_bar_width Numeric. Line width of the boxplot error bars.
#'   Default is \code{0.3}.
#' @param box_pos Numeric vector of length 2. Position of the inset boxplot
#'   as \code{c(x, y)} coordinates in normalized plot coordinates (0-1).
#'   Default is \code{c(0.99, 0.5)}.
#' @param box_width_height Numeric vector of length 2. Width and height of
#'   the inset boxplot as \code{c(width, height)} in normalized coordinates.
#'   Default is \code{c(0.3, 0.8)}.
#' @param return_data Logical. If \code{TRUE}, returns the processed data frame
#'   instead of plots and statistics. Useful for custom downstream analysis.
#'   Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @return
#' The return value depends on the \code{return_data} parameter:
#' \itemize{
#'   \item If \code{return_data = FALSE} (default): A list containing:
#'     \describe{
#'       \item{\code{plot}}{A combined plot showing cumulative distribution
#'         functions of motif occupancy, optionally with inset boxplots}
#'       \item{\code{statistics}}{A data frame with statistical test results
#'         from pairwise comparisons between sample groups}
#'     }
#'   \item If \code{return_data = TRUE}: A data frame with columns:
#'     \describe{
#'       \item{\code{sample}}{Sample identifier}
#'       \item{\code{sample_group}}{Sample group identifier}
#'       \item{\code{rname}}{Transcript name}
#'       \item{\code{rel}}{Relative codon position}
#'       \item{\code{codon_seq}}{Motif sequence}
#'       \item{\code{value}}{Normalized occupancy value}
#'       \item{\code{freq}}{Frequency of the motif across transcripts}
#'     }
#' }
#'
#' @seealso
#' \code{\link{do_offset_correction}} for P-site offset correction methods
#'
#' @examples
#' \dontrun{
#' # Basic amino acid motif analysis
#' result <- motif_occupancy(
#'   object = ribo_data,
#'   cds_fa = "path/to/cds_sequences.fa",
#'   motif_pattern = "PPP",
#'   search_type = "amino"
#' )
#'
#' # View the plot
#' result$plot
#'
#' # Check statistical results
#' result$statistics
#'
#' # Codon-based analysis with multiple motifs
#' codon_result <- motif_occupancy(
#'   object = ribo_data,
#'   cds_fa = "path/to/cds_sequences.fa",
#'   motif_pattern = c("CCC", "CCG", "CCA", "CCT"),
#'   search_type = "codon",
#'   min_counts = 100,
#'   exclude_length = c(30, 30)
#' )
#'
#' # Return raw data for custom analysis
#' raw_data <- motif_occupancy(
#'   object = ribo_data,
#'   cds_fa = "path/to/cds_sequences.fa",
#'   motif_pattern = "KKK",
#'   return_data = TRUE
#' )
#'
#' # Analysis with offset correction
#' corrected_result <- motif_occupancy(
#'   object = ribo_data,
#'   cds_fa = "path/to/cds_sequences.fa",
#'   motif_pattern = "DDD",
#'   do_offset_correct = TRUE,
#'   position_shift = -1
#' )
#' }
#'
#' @export
#' @rdname motif_occupancy
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
                   exclude_length = c(45,45),
                   min_counts = 64,
                   motif_pattern = "PPP",
                   ref_group = NULL,
                   insert_box = TRUE,
                   box_width = 0.25,
                   error_bar_width = 0.3,
                   box_pos = c(0.99,0.5),
                   box_width_height = c(0.3,0.8),
                   return_data = FALSE){
            search_type <- match.arg(search_type,choices = c("amino", "codon"))
            # ==================================================================
            # whether do reads offset correction
            if(do_offset_correct == TRUE){
              sry <- do_offset_correction(object = object,shift = position_shift)
            }else{
              sry <- object@summary_info
            }

            sry <- sry %>% fastplyr::f_filter(mstart != 0 | mstop != 0)


            # expect reads per position
            avg.ct <- sry %>%
              dplyr::mutate(relst = pos - mstart, relsp = pos - mstop) %>%
              fastplyr::f_filter(relst > exclude_length[1] & relsp < -exclude_length[2]) %>%
              dplyr::mutate(cdslen = mstop - mstart + 1) %>%
              dplyr::group_by(sample,rname,cdslen) %>%
              dplyr::summarise(counts = sum(count)) %>%
              dplyr::mutate(avg_ct = counts/cdslen) %>%
              dplyr::select(sample,rname,counts,avg_ct)

            # average reads
            pltdf <- sry %>%
              dplyr::inner_join(y = avg.ct,by = c("sample", "rname")) %>%
              # filter low counts
              fastplyr::f_filter(counts > min_counts) %>%
              dplyr::mutate(rel = pos - mstart, norm = count/avg_ct) %>%
              fastplyr::f_filter(rel >= 0 & rel <= (mstop - mstart + 1)) %>%
              fastplyr::f_group_by(sample,sample_group,rname,rel) %>%
              fastplyr::f_summarise(normsm = sum(norm)) %>%
              # codon position
              dplyr::mutate(rel = (rel %/% 3) + 1) %>%
              fastplyr::f_group_by(sample,sample_group,rname,rel) %>%
              fastplyr::f_summarise(value = mean(normsm))

            # ==================================================================
            # check search_type
            if(search_type == "amino"){
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
                    dplyr::mutate(pos = (pos - 1) %/% 3)
                }

                return(motif.pos)
              }) -> motif.pos

              # ==================================================================
              # calculate relative occupancy
              pltdf3 <- pltdf %>%
                dplyr::inner_join(y = motif.pos,by = "rname",relationship = "many-to-many") %>%
                dplyr::filter(rel == pos) %>%
                dplyr::rename(codon_seq = motif)

            }else{
              # filter cds sequence
              # load cds fasta
              cds <- Biostrings::readDNAStringSet(cds_fa)

              # remove transcript contains Ns in sequence
              valid_cds <- cds[!grepl("[^ACGTacgt]", as.character(cds))]

              ids_retain <- intersect(unique(sry$rname),names(valid_cds))
              valid_cds <- valid_cds[ids_retain]

              # Filter sequences whose length is a multiple of 3
              valid_cds <- valid_cds[Biostrings::width(valid_cds) %% 3 == 0]

              # loop to extract codon seqeunce
              lapply(seq_along(valid_cds),function(x){
                tmp <- valid_cds[x]
                cdslen <- Biostrings::width(tmp)
                interval <- seq(1, cdslen, by = 3)

                seqs <- stringr::str_sub_all(tmp,start = interval, end = interval + 2)[[1]]

                data.frame(rname = names(tmp),rel = 1:(cdslen/3),codon_seq = seqs)
              }) %>% do.call("rbind",.) %>%
                data.frame() %>%
                dplyr::add_count(codon_seq, name = "freq") -> codon_info

              # ==================================================================

              # merge with occupancy
              pltdf2 <- pltdf %>%
                fastplyr::f_inner_join(y = codon_info,by = c("rname", "rel")) %>%
                fastplyr::f_group_by(sample,sample_group,codon_seq,freq)

              # ==================================================================
              # amino acid annotation
              aa_info <- get_aa_table()

              pltdf3 <- pltdf2 %>%
                fastplyr::f_inner_join(y = aa_info,by = c("codon_seq" = "Codon")) %>%
                dplyr::filter(codon_seq %in% motif_pattern)
            }

            # ==================================================================

            # whether aggregate replicates
            if(merge_rep == TRUE){
              pltdf3 <- pltdf3 %>%
                fastplyr::f_group_by(sample_group,rname,rel,codon_seq) %>%
                fastplyr::f_summarise(value = mean(value)) %>%
                dplyr::rename(sample = sample_group)
            }

            # ===========================================================================
            # stastics
            if (!requireNamespace("ggpubr", quietly = TRUE)) {
              stop("Package 'ggpubr' is required. Please install it.")
            }

            stcs <- ggpubr::compare_means(data = pltdf3,
                                          formula = value ~ sample,
                                          group.by = "codon_seq",
                                          ref.group = ref_group,
                                          method = "wilcox.test")

            # ===========================================================================
            # plot

            # loop for each motif
            mt <- unique(pltdf3$codon_seq)

            # x = 1
            lapply(seq_along(mt),function(x){
              pdf <- subset(pltdf3,codon_seq == mt[x])

              pbg <-
                ggplot(pdf) +
                stat_ecdf(mapping = aes(x = log2(value),color = sample)) +
                theme_bw(base_size = 12) +
                facet_wrap(~codon_seq) +
                theme(panel.grid = element_blank(),
                      # strip.background = element_blank(),
                      strip.text = element_text(face = "bold",size = rel(1)),
                      axis.text = element_text(colour = "black")) +
                ylab("Cumulative fracton") +
                xlab("Log2 (codon/amino occupancy)")

              # whether inset boxplot
              if(insert_box == TRUE){
                pin <-
                  ggplot(pdf,
                         aes(x = sample,y = log2(value),fill = sample)) +
                  stat_boxplot(geom = "errorbar", width = box_width,linewidth = error_bar_width) +
                  geom_boxplot(outliers = F,show.legend = F,width = 0.75) +
                  theme_bw() +
                  theme(panel.grid = element_blank(),
                        plot.background = element_blank(),
                        panel.background = element_blank(),
                        axis.text = element_text(colour = "black"),
                        axis.line.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.ticks.x = element_blank()) +
                  ylab("Log2 (occupancy)")


                df <- tibble::tibble(x = box_pos[1], y = box_pos[2], plot = list(pin))

                if (requireNamespace("ggpp", quietly = TRUE)) {
                  pcb <- pbg +
                    expand_limits(x = 0, y = 0) +
                    ggpp::geom_plot_npc(data = df, aes(npcx = x, npcy = y, label = plot),
                                        vp.width = box_width_height[1], vp.height = box_width_height[2])
                } else {
                  warning("Package 'ggpp' is needed for this function to work.")
                }

              }else{
                pcb <- pbg
              }

              return(pcb)
            }) -> plist

            if (requireNamespace("cowplot", quietly = TRUE)) {
              cmb <- cowplot::plot_grid(plotlist = plist)
            } else {
              warning("Package 'cowplot' is needed for this function to work.")
            }


            # return
            if(return_data == FALSE){
              return(list(plot = cmb,statistics = stcs))
            }else{
              return(pltdf3)
            }
          }
)
