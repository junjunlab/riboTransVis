
#' Analyze Codon or Amino Acid Motif Occupancy from Ribosome Profiling Data
#'
#' This function calculates and visualizes motif-level translational occupancy from a \code{ribotrans} object, assessing motif-specific ribosome densities at codon or amino acid resolution. It provides options for offset correction, motif searching, replicate merging, and boxplot insets.
#'
#' @param object A \code{ribotrans} object containing preprocessed ribosome profiling data with corresponding annotations in \code{@features} and positional read density in \code{@summary_info}.
#' @param merge_rep Logical. If \code{TRUE}, replicates from the same \code{sample_group} will be averaged before visualization. Default is \code{FALSE}.
#' @param cds_fa Path to the CDS fasta file with transcript nucleotide sequences.
#' @param do_offset_correct Logical. If \code{TRUE}, perform read offset correction using \code{do_offset_correction()}. Default is \code{FALSE}.
#' @param position_shift Integer. Number of nucleotides used to correct P-site position. Default is 0.
#' @param search_type Search level: "amino" (default) or "codon". Determines whether the motif is searched at the amino acid or codon sequence level.
#' @param exclude_length Integer vector of length 2. Number of codons to exclude from the start and end of CDS. Default is \code{c(45, 45)}.
#' @param min_counts Integer. Minimum read count per transcript to include for analysis. Default is 64.
#' @param motif_pattern Character vector. Motif(s) (amino acids or nucleotides) to search for. Default is "PPP".
#' @param insert_box Logical. Whether to include an inset boxplot in ECDF plots. Default is \code{TRUE}.
#' @param box_pos Numeric vector of length 2. The relative position (npc coordinates) of the inset box. Default is \code{c(0.99, 0.5)}.
#' @param box_width_height Numeric vector of length 2. Absolute width and height of the inset box. Default is \code{c(0.3, 0.8)}.
#' @param return_data Logical. If \code{TRUE}, return data frame instead of plotting. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This function determines transcript positions matching specified motifs (on codon or amino acid level) and evaluates their normalized ribosome occupancy. Occupancy is calculated as read density normalized by transcript average. Plots display the distribution of log2 occupancy for each motif using ECDF curves. Optionally, inset boxplots illustrate sample-wise distributions.
#'
#' The analysis includes CDS filtering, motif location extraction, optional offset correction, normalization of read counts, and final plotting with optional replicate merging.
#'
#' Requires packages: \code{Biostrings}, \code{fastplyr}, \code{data.table}, \code{ggplot2}, \code{cowplot}, and optionally \code{ggpp}.
#'
#' @return If \code{return_data = FALSE}, returns a \code{ggplot2} object with motif occupancy plots. Otherwise returns a data frame of motif positions and relative occupancy values.
#'
#'
#' @examples
#' \dontrun{
#' motif_occupancy(object = ribo_obj,
#' cds_fa = "cds.fa",
#' motif_pattern = c("PPP", "GGG"),
#' search_type = "amino")
#' }
#'
#' @export
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
                   insert_box = TRUE,
                   box_pos = c(0.99,0.5),
                   box_width_height = c(0.3,0.8),
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
                  dplyr::mutate(pos = (pos - 1) %/% 3)
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
            # p <-

            # loop for each motif
            mt <- unique(rel2motif.df$motif)

            lapply(seq_along(mt),function(x){
              pdf <- subset(rel2motif.df,motif == mt[x])

              pbg <-
                ggplot(pdf) +
                stat_ecdf(mapping = aes(x = log2(value),color = sample)) +
                theme_bw(base_size = 12) +
                facet_wrap(~motif) +
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
                  stat_boxplot(geom = "errorbar", width = 0.5,linewidth = 0.25) +
                  geom_boxplot(outliers = F,show.legend = F,width = 0.75) +
                  theme_classic() +
                  theme(panel.grid = element_blank(),
                        plot.background = element_blank(),
                        panel.background = element_blank(),
                        axis.line.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.ticks.x = element_blank()) +
                  ylab("Log2 (codon/amino occupancy)")


                df <- tibble(x = box_pos[1], y = box_pos[2], plot = list(pin))

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
              return(cmb)
            }else{
              return(rel2motif.df)
            }
          }
)
