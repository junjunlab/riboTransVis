# ==============================================================================
# enrichment plot for single gene
# ==============================================================================

#' @title Enrichment Profile Plot for SERP Object
#' @description
#' Visualizes the ribosome profiling signal enrichment (IP / total) from a \code{serp} object.
#' For each transcript, the function plots the enrichment curve and uncertainty shading,
#' optionally with customized signal ranges, CDS boundary annotation, and faceting layout.
#' Requires \code{get_enrichment} to be run beforehand.
#'
#' @param object An object of class \code{serp}. Must contain enrichment data calculated via \code{get_enrichment}.
#' @param cds_line A vector specifying the color and line width of the CDS start marker. Default: \code{c("#CC6600", 4)}.
#' @param line_col Color of enrichment line and the fill of confidence interval rectangles. Default: \code{"#377CC0"}.
#' @param shadow_alpha A numeric vector of length 2 for the alpha transparency range (min, max) of the enrichment CI shading. Default: \code{c(0.1, 0.3)}.
#' @param enrich_threshold A numeric scalar defining the horizontal dashed line threshold for enrichment (e.g., 1 = equal signal). Default: \code{1}.
#' @param new_signal_range Logical, whether to annotate a textual enrichment signal range label (from 0 to max upper). Default: \code{FALSE}.
#' @param range_x X position of the signal range label, in NPC units (0–1). Used only if \code{new_signal_range = TRUE}. Default: \code{0.99}.
#' @param range_y Y position of the signal range label, in NPC units (0–1). Default: \code{0.98}.
#' @param range_size Text size of the signal range annotation. Default: \code{4}.
#' @param range_digit Number of digits to round in signal range label. Default: \code{1}.
#' @param facet A ggplot2 faceting specification, e.g., \code{facet_grid(sample~rname)}. Default: \code{facet_grid(sample~rname, switch = "y")}.
#' @param nrow,ncol Number of rows and columns in multi-figure layout when plotting multiple transcripts. Default: \code{NULL} (automatic).
#' @param ... Additional arguments (currently unused).
#'
#' @return A patchwork of enrichment line plots (using cowplot::plot_grid), one per transcript.
#' Additionally shows enrichment variability via shaded area and optionally range labeling and CDS marks.
#'
#'
#' @importFrom ggside geom_xsidesegment scale_xsidey_continuous ggside
#' @importFrom fastplyr f_group_by f_summarise
#'
#'
#' @seealso \code{\link{get_enrichment}} for enrichment computation.
#'
#' @examples
#' \dontrun{
#' # Assuming `x` is a serp object with enrichment already computed
#' enrichment_plot(x)
#' enrichment_plot(x, cds_line = c("red", 2), enrich_threshold = 1.5, new_signal_range = TRUE)
#' }
#'
#' @export
setGeneric("enrichment_plot",function(object,...) standardGeneric("enrichment_plot"))





#' @rdname enrichment_plot
#' @export
setMethod("enrichment_plot",
          signature(object = "serp"),
          function(object,
                   cds_line = c("#CC6600",4),
                   line_col = "#377CC0",
                   shadow_alpha = c(0.1,0.3),
                   enrich_threshold = 1,
                   new_signal_range = F,
                   range_x = 0.99,
                   range_y = 0.98,
                   range_size = 4,
                   range_digit = 1,
                   facet = facet_grid(sample~rname,switch = "y"),
                   nrow = NULL, ncol = NULL){
            # check data
            if(nrow(object@enriched_ratio) == 0){
              stop("Please run `get_enrichment` first!")
            }

            plotdf <- object@enriched_ratio

            # check data
            feature <- subset(object@features, gene == object@gene_name)

            # window size
            bin <- plotdf$window_size[1]

            tid <- feature$idnew

            # loop for each transcript_id
            # x = 1
            lapply(seq_along(tid),function(x){
              tmp_feature <- subset(feature, idnew == tid[x])
              tmp <- subset(plotdf, rname == tid[x])

              # ================================================================
              # whether add new signal range
              if(new_signal_range == TRUE){
                axis.text.y = element_blank()
                axis.ticks.y = element_blank()

                range <- paste("[0-",round(max(tmp$upper),digits = range_digit),"]",sep = "")

                if (requireNamespace("ggpp", quietly = TRUE)) {
                  range_label <- ggpp::annotate(geom = "text_npc",
                                                npcx = range_x,npcy = range_y,
                                                label = range,
                                                size = range_size)
                } else {
                  warning("Package 'ggpp' is needed for this function to work.")
                }

              }else{
                axis.text.y = NULL
                axis.ticks.y = NULL
                range_label <- NULL
              }

              # ================================================================
              # plot
              ggplot(tmp) +
                geom_rect(aes(xmin = (pos+bin/2)/3-0.5,
                              xmax = (pos+bin/2)/3+0.5,
                              ymin = lower,
                              ymax = upper,
                              alpha = 1/((1/win_total+1/win_ip)*bin) ),
                          fill = line_col,show.legend = F) +
                geom_line(aes(x = (pos+bin/2)/3, y = enrich), col = line_col) +
                geom_hline(yintercept = enrich_threshold,lty = "dashed",color = "black") +
                range_label +
                # scale_y_continuous(trans = "log2", oob = scales::squish, expand = c(0,0)) +
                scale_alpha_continuous(trans = "sqrt", name = "prec", range = shadow_alpha) +
                facet +
                theme_bw() +
                theme(panel.grid = element_blank(),
                      axis.text = element_text(colour = "black"),
                      strip.text.y.left = element_text(angle = 0),
                      strip.background = element_blank(),
                      strip.text = element_text(face = "bold"),
                      strip.placement = "outside",
                      axis.text.y = axis.text.y,
                      axis.ticks.y = axis.ticks.y,
                      ggside.panel.background = element_blank(),
                      ggside.panel.border = element_blank()) +
                xlab("Ribosome position\n(codons / amino acids)") +
                ylab("Mean enrichment\n(IP / total)") +
                ggside::geom_xsidesegment(data = data.frame(x = 1,xend = tmp_feature$cds/3),
                                          aes(x = x,xend = xend,
                                              y = 0.5,yend = 0.5),
                                          inherit.aes = F,
                                          linewidth = as.numeric(cds_line[2]),color = cds_line[1]) +
                ggside::scale_xsidey_continuous(breaks = NULL) +
                # scale_x_continuous(expand = c(0,0)) +
                ggside::ggside(collapse = "x")


            }) -> plist

            # combine plot list
            if (requireNamespace("cowplot", quietly = TRUE)) {
              cowplot::plot_grid(plotlist = plist,nrow = nrow,ncol = ncol)
            } else {
              warning("Package 'cowplot' is needed for this function to work.")
            }
          }
)



# ==============================================================================
# enrichment plot2 for single gene
# ==============================================================================

#' Plot Enrichment Profiles for Target Gene
#'
#' This function visualizes the enrichment signal (IP RPM / Total RPM) across transcripts
#' in a \code{serp} object. It supports sample replication handling, signal smoothing,
#' display in codon or nucleotide resolution, and graphical customization.
#'
#' @param object A \code{serp} object generated and processed by previous analysis steps (must include \code{total_occupancy}, \code{ip_occupancy}, and \code{features}).
#' @param merge_rep Logical; whether to merge biological replicates (sample groups) for plotting. Default is \code{FALSE}.
#' @param var_alpha Alpha transparency level for shaded error region (standard deviation) when \code{merge_rep = TRUE}. Default is 0.5.
#' @param smooth Logical; whether to apply rolling mean smoothing over enrichment signal. Default is \code{TRUE}.
#' @param window_size Integer; sliding window size for smoothing (number of nucleotides or codons depending on \code{mode}). Default is 45.
#' @param retain_cds Logical; if \code{TRUE}, restricts plot to CDS (coding sequence) regions of transcripts. Default is \code{TRUE}.
#' @param mode Character; plotting resolution, either \code{"codon"} (per 3 nt) or \code{"nt"} (single nucleotide). Default is \code{"codon"}.
#' @param new_signal_range Logical; whether to display signal range annotation on the right-top corner of plot. Default is \code{FALSE}.
#' @param enrich_threshold Numeric; add a horizontal dashed line indicating enrichment threshold (useful to visually detect high/low enrichment). Default is \code{1}.
#' @param range_x,range_y Coordinates used for positioning signal range label when \code{new_signal_range = TRUE}. Default to \code{0.99}, \code{0.98}.
#' @param range_size Size of signal range annotation label. Default is \code{4}.
#' @param range_digit Integer; number of digits to round displayed max enrichment range. Default is \code{1}.
#' @param exon_line A vector specifying color and size (e.g., \code{c("grey", 3)}) for exon structure segments if \code{mode = "nt"}. Default is \code{c("grey", 3)}.
#' @param cds_line A vector specifying color and size (e.g., \code{c("grey30", 7)}) for CDS structure line. Default is \code{c("grey30", 7)}.
#' @param facet A \code{ggplot2} facet specification, default is \code{facet_grid(sample ~ rname, switch = "y")}.
#' @param nrow,ncol Set layout for subplot arrangement when multiple transcripts or samples are included (used in cowplot::plot_grid).
#' @param return_data Logical; if \code{TRUE}, return the processed enrichment data (data.frame) used for plotting. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @details
#' The function uses annotation in the \code{serp} object to generate one or more plots showing enrichment (IP/Total) across exons or CDS regions.
#'
#' When \code{smooth = TRUE}, rolling mean is computed using \code{zoo::rollsum} with user-defined \code{window_size}.
#'
#' Samples can be grouped by \code{sample_group}, and replicate merging with SD shading is available.
#'
#' If \code{return_data = TRUE}, returns a data.frame with values including smoothed counts and enrichment ratio.
#'
#' @return A \code{ggplot} object showing the enrichment profile plot, or a data.frame when \code{return_data=TRUE}.
#'
#' @importFrom dplyr select left_join mutate filter rename group_by summarise
#' @importFrom purrr map_df
#' @importFrom ggplot2 ggplot geom_line geom_ribbon geom_hline facet_grid facet_wrap element_blank element_text theme theme_bw xlab ylab aes
#' @importFrom ggside ggside geom_xsidesegment scale_xsidey_continuous
#' @importFrom stats sd
#'
#'
#' @examples
#' \dontrun{
#' # Visualize enrichment with smoothing and signal range label
#' enrichment_plot2(serp_obj, smooth = TRUE, new_signal_range = TRUE)
#'
#' # Only plot codon positions on CDS without smoothing
#' enrichment_plot2(serp_obj, smooth = FALSE, retain_cds = TRUE, mode = "codon")
#'
#' # Return processed data
#' data <- enrichment_plot2(serp_obj, return_data = TRUE)
#' }
#'
#' @export
setGeneric("enrichment_plot2",function(object,...) standardGeneric("enrichment_plot2"))





#' @rdname enrichment_plot2
#' @export
setMethod("enrichment_plot2",
          signature(object = "serp"),
          function(object,
                   merge_rep = FALSE,
                   var_alpha = 0.5,
                   smooth = TRUE,
                   window_size = 45,
                   retain_cds = TRUE,
                   mode = c("codon", "nt"),
                   new_signal_range = FALSE,
                   enrich_threshold = 1,
                   range_x = 0.99,
                   range_y = 0.98,
                   range_size = 4,
                   range_digit = 1,
                   exon_line = c("grey",3),
                   cds_line = c("grey30",7),
                   facet = facet_grid(sample~rname,switch = "y"),
                   nrow = NULL, ncol = NULL,
                   return_data = FALSE){
            mode <- match.arg(mode,choices = c("codon", "nt"))
            # ==================================================================
            # ribo rpm/rna rpm for each position

            # check rna coverage
            if(is.null(object@total_occupancy)){
              stop("Please supply total_occupancy data!")
            }

            if(is.null(object@ip_occupancy)){
              stop("Please supply ip_occupancy data!")
            }

            # run
            mb <- object@total_occupancy %>%
              dplyr::select(sample,sample_group,rname,pos,count,rpm) %>%
              dplyr::left_join(y = object@ip_occupancy[,c("sample","sample_group","rname","pos","count","rpm")],
                               by = c("sample","sample_group","rname","pos"))

            mb[is.na(mb)] <- 0

            # loop for each sample and transcript to smooth
            tanno <- subset(object@features, gene == object@gene_name)

            # check whether smooth data
            # x = 1
            purrr::map_df(1:nrow(tanno),function(x){
              tmp <- tanno[x,]

              tmp2 <- mb %>%
                dplyr::filter(rname == tmp$idnew)

              sp <- unique(tmp2$sample)

              # loop for each sample
              # s = 1
              purrr::map_df(seq_along(sp),function(s){
                tmp3 <- tmp2 %>%
                  dplyr::filter(sample == sp[s])

                # each pos
                pos.df <- data.frame(rname = tmp3$rname[1],pos = 1:tmp$exonlen)

                # merge with value
                pos.df <- pos.df %>% dplyr::left_join(y = tmp3,by = c("rname", "pos"))
                pos.df$sample <- tmp3$sample[1]
                pos.df$sample_group <- tmp3$sample_group[1]
                pos.df[is.na(pos.df)] <- 0

                # retain CDS region
                if(retain_cds == TRUE){
                  pos.df <- pos.df %>%
                    dplyr::filter(pos >= tmp$mstart & pos <= tmp$mstop) %>%
                    dplyr::mutate(pos = pos - tmp$mstart + 1)

                  # codon pos
                  if(mode == "codon"){
                    pos.df <- pos.df %>%
                      dplyr::mutate(pos = (pos - 1) %/% 3 + 1) %>%
                      fastplyr::f_group_by(rname, sample, sample_group,pos) %>%
                      fastplyr::f_summarise(rpm.x = sum(rpm.x), rpm.y = sum(rpm.y))
                  }
                }

                # soomth data for each position
                if(smooth == TRUE){
                  if (requireNamespace("zoo", quietly = TRUE)) {
                    pos.df$sm1 <- zoo::rollsum(pos.df$rpm.x,k = window_size,fill = 0,align = "center")
                    pos.df$sm2 <- zoo::rollsum(pos.df$rpm.y,k = window_size,fill = 0,align = "center")
                  } else {
                    warning("Package 'zoo' is needed for this function to work.")
                  }


                  pos.df$smratio <- pos.df$sm2/pos.df$sm1
                }else{
                  pos.df$smratio <- pos.df$rpm.y/pos.df$rpm.x
                }

                pos.df <- pos.df %>% tidyr::replace_na(list(smratio = 0))

                return(pos.df)
              }) -> smooth.df

              return(smooth.df)
            }) -> mb.smoothed

            # ==================================================================
            # plot
            if(merge_rep == TRUE){
              mb.smoothed <- mb.smoothed %>%
                fastplyr::f_group_by(sample_group,rname,pos) %>%
                fastplyr::f_summarise(rpm.x = mean(rpm.x),rpm.y = mean(rpm.y),
                                      sm1 = mean(sm1),sm2 = mean(sm2),
                                      sd = stats::sd(smratio),
                                      smratio = mean(smratio)
                ) %>%
                dplyr::rename(sample = sample_group)

              shadow <- geom_ribbon(aes(ymin = smratio - sd,
                                        ymax = smratio + sd,
                                        x = pos,y = smratio,
                                        fill = sample), alpha = var_alpha)
            }else{
              shadow <- NULL
            }
            # ================================================================
            tid <- tanno$idnew

            # loop for each transcript_id
            # x = 1
            lapply(seq_along(tid),function(x){
              tmp_feature <- subset(tanno, idnew == tid[x])
              tmp <- subset(mb.smoothed, rname == tid[x])

              # whether add new signal range
              if(new_signal_range == TRUE){
                axis.text.y = element_blank()
                axis.ticks.y = element_blank()

                range <- paste("[0-",round(max(tmp$smratio),digits = range_digit),"]",sep = "")

                if (requireNamespace("ggpp", quietly = TRUE)) {
                  range_label <- ggpp::annotate(geom = "text_npc",
                                                npcx = range_x,npcy = range_y,
                                                label = range,
                                                size = range_size)
                } else {
                  warning("Package 'ggpp' is needed for this function to work.")
                }

              }else{
                axis.text.y = NULL
                axis.ticks.y = NULL
                range_label <- NULL
              }

              # structure
              if(mode == "codon"){
                struc_cds <- data.frame(x = 1,xend = tmp_feature$cds/3)

                sly2 <- ggside::geom_xsidesegment(data = struc_cds,
                                                  aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                  inherit.aes = F,
                                                  linewidth = as.numeric(cds_line[2]),color = cds_line[1])

                sly1 <- NULL
              }else{
                struc_exon <- data.frame(x = 1,xend = tmp_feature$exonlen)

                sly1 <- ggside::geom_xsidesegment(data = struc_exon,
                                                  aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                  inherit.aes = F,
                                                  linewidth = as.numeric(exon_line[2]),color = exon_line[1])

                struc_cds <- data.frame(x = tmp_feature$mstart,xend = tmp_feature$mstop)

                sly2 <- ggside::geom_xsidesegment(data = struc_cds,
                                                  aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                  inherit.aes = F,
                                                  linewidth = as.numeric(cds_line[2]),color = cds_line[1])
              }

              # xlabel
              if(mode == "codon"){
                xlabel <- "Ribosome position (codons / amino acids)"
              }else{
                xlabel <- "Ribosome position (nucleotides)"
              }

              # ================================================================
              p <-
                ggplot(tmp) +
                shadow +
                geom_line(aes(x = pos,y = smratio, color = sample)) +
                geom_hline(yintercept = enrich_threshold,lty = "dashed",color = "black") +
                range_label +
                facet +
                theme_bw() +
                theme(panel.grid = element_blank(),
                      axis.text = element_text(colour = "black"),
                      strip.text.y.left = element_text(angle = 0, hjust = 1),
                      strip.background = element_blank(),
                      strip.text = element_text(face = "bold"),
                      strip.placement = "outside",
                      axis.text.y = axis.text.y,
                      axis.ticks.y = axis.ticks.y,
                      ggside.panel.background = element_blank(),
                      ggside.panel.border = element_blank()) +
                xlab(xlabel) +
                ylab("Mean enrichment (IP / total)[AU]") +
                # add structure
                sly1 + sly2 +
                ggside::scale_xsidey_continuous(breaks = NULL) +
                # scale_x_continuous(expand = c(0,0)) +
                ggside::ggside(collapse = "x")

              return(p)
            }) -> plist

            # combine plot list
            if (requireNamespace("cowplot", quietly = TRUE)) {
              cmb <- cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol)
            } else {
              warning("Package 'cowplot' is needed for this function to work.")
            }

            # return
            if(return_data == FALSE){
              return(cmb)
            }else{
              return(mb.smoothed)
            }
          }
)





# ==============================================================================
# transcript for single gene
# ==============================================================================


#' Plot translation profiles
#'
#' @description
#' This is a generic function for plotting translation-related profiles (such as
#' ribosome occupancy or enrichment signals) for different types of transcriptome data.
#'
#' @param object An object of class \code{serp} or \code{ribotrans}, containing necessary data for plotting.
#' @param ... Additional arguments passed to the specific method.
#'
#' @return A plot, typically a ggplot2 object.
#'
#' @export
setGeneric("trans_plot",function(object,...) standardGeneric("trans_plot"))




#' Plot Ribosome density Across Transcripts
#'
#' @description
#' Visualize interactive ribosome profiling (IP) or total ribosome (total) signal across
#' the transcript body for one or more transcripts from a \code{serp} object. This function
#' supports codon- or nucleotide-resolution visualization, replicate merging, structural
#' annotations (e.g., exon and CDS information), and optional summary statistics.
#'
#' @param object A \code{serp} object. Must contain the slots \code{total_occupancy}, \code{ip_occupancy}, and \code{features}.
#' @param merge_rep Logical. If \code{TRUE}, biological replicates are merged (by \code{sample_group}) and plotted with mean values and shaded standard deviation ribbons. Default: \code{FALSE}.
#' @param var_alpha Numeric. Transparency level for the shaded ribbon when \code{merge_rep = TRUE}. Default: 0.5.
#' @param retain_cds Logical. If \code{TRUE}, restricts plotting to coding sequence (CDS) regions only. Default: \code{TRUE}.
#' @param mode Character. Resolution of x-axis: either \code{"codon"} (default) for amino acid-level or \code{"nt"} for nucleotide-level.
#' @param new_signal_range Logical. Whether to annotate the plot with the observed signal range. Default: \code{FALSE}.
#' @param enrich_threshold Numeric. Not used in current implementation but included for compatibility with other functions. Default: 1.
#' @param range_x X position (normalized 0–1) for the signal range label. Default: 0.98.
#' @param range_y Y position (normalized 0–1) for the signal range label. Default: 0.98.
#' @param range_size Numeric. Font size of the signal range label. Default: 4.
#' @param range_digit Integer. Number of digits to round the signal range annotation. Default: 1.
#' @param exon_line Character vector of length 2 (e.g., \code{c("grey", 3)}). Specifies the color and line width of the exon annotation in \code{"nt"} mode.
#' @param cds_line Character vector of length 2 (e.g., \code{c("grey30", 7)}). Specifies the color and line width of the CDS annotation.
#' @param facet ggplot2 facet. Default is \code{facet_grid(sample ~ rname, switch = "y")}.
#' @param nrow Number of rows in the final layout when multiple transcripts are plotted. Passed to \code{cowplot::plot_grid()}.
#' @param ncol Number of columns in the final layout when multiple transcripts are plotted. Passed to \code{cowplot::plot_grid()}.
#' @param return_data Logical. If \code{TRUE}, return the underlying processed data instead of plot. Default: \code{FALSE}.
#'
#'
#' @details
#' This function aggregates and visualizes transcript-level signal for both IP (immunoprecipitation)
#' and total RNA data stored in a \code{serp} object. It supports CDS-only rendering, codon resolution,
#' and shows structure annotations such as exon and coding regions using \code{ggside} extensions.
#'
#' If \code{merge_rep = TRUE}, signals for replicate samples belonging to the same \code{sample_group}
#' are averaged, and standard deviation is visualized using \code{geom_ribbon()}.
#'
#' The object returned is either a \code{ggplot} object or a \code{data.frame} if \code{return_data = TRUE}.
#'
#' Requires pre-computed smooth signal (column named \code{smooth}) within the IP and total input data.
#'
#' @return A \code{ggplot} object or a \code{data.frame} if \code{return_data = TRUE}.
#'
#' @importFrom dplyr filter mutate left_join rename
#' @importFrom purrr map_df
#' @importFrom fastplyr f_group_by f_summarise
#' @importFrom ggplot2 ggplot geom_line geom_ribbon facet_grid theme element_blank element_text theme_bw xlab ylab aes
#' @importFrom ggside ggside geom_xsidesegment scale_xsidey_continuous
#'
#'
#'
#' @examples
#' \dontrun{
#' # Plot ribosome and total signal at codon resolution
#' trans_plot(serp_obj, mode = "codon")
#'
#' # Merge replicates and plot with standard deviation shading
#' trans_plot(serp_obj, merge_rep = TRUE)
#'
#' # Return processed data without plotting
#' data <- trans_plot(serp_obj, return_data = TRUE)
#' }
#' @export
setMethod("trans_plot",
          signature(object = "serp"),
          function(object,
                   merge_rep = FALSE,
                   var_alpha = 0.5,
                   retain_cds = TRUE,
                   mode = c("codon", "nt"),
                   new_signal_range = FALSE,
                   enrich_threshold = 1,
                   range_x = 0.99,
                   range_y = 0.98,
                   range_size = 4,
                   range_digit = 1,
                   exon_line = c("grey",3),
                   cds_line = c("grey30",7),
                   facet = facet_grid(sample~rname,switch = "y"),
                   nrow = NULL, ncol = NULL,
                   return_data = FALSE){
            mode <- match.arg(mode,choices = c("codon", "nt"))
            # ==================================================================
            # ribo rpm/rna rpm for each position

            # check rna coverage
            if(is.null(object@total_occupancy)){
              stop("Please supply total_occupancy data!")
            }else{
              tt <- object@total_occupancy
              tt$type <- "Total"
            }

            if(is.null(object@ip_occupancy)){
              stop("Please supply ip_occupancy data!")
            }else{
              ip <- object@ip_occupancy
              ip$type <- "IP"
            }

            # combine

            mb <- rbind(tt, ip)

            # loop for each sample and transcript to smooth
            tanno <- subset(object@features, gene == object@gene_name)

            # check whether smooth data
            # x = 1
            purrr::map_df(1:nrow(tanno),function(x){
              tmp <- tanno[x,]

              tmp2 <- mb %>%
                dplyr::filter(rname == tmp$idnew)

              sp <- unique(tmp2$sample)

              # loop for each sample
              # s = 1
              purrr::map_df(seq_along(sp),function(s){
                tmp3 <- tmp2 %>%
                  dplyr::filter(sample == sp[s])

                # loop for type
                tp <- unique(tmp3$type)

                # t = 1
                purrr::map_df(seq_along(tp),function(t){
                  tmp4 <- tmp3 %>% dplyr::filter(type == tp[t])

                  # each pos
                  pos.df <- data.frame(rname = tmp4$rname[1],pos = 1:tmp$exonlen,type = tp[t])

                  # merge with value
                  pos.df <- pos.df %>% dplyr::left_join(y = tmp4,by = c("rname", "pos", "type"))
                  pos.df$sample <- tmp4$sample[1]
                  pos.df$sample_group <- tmp4$sample_group[1]
                  pos.df[is.na(pos.df)] <- 0

                  # retain CDS region
                  if(retain_cds == TRUE){
                    pos.df <- pos.df %>%
                      dplyr::filter(pos >= tmp$mstart & pos <= tmp$mstop) %>%
                      dplyr::mutate(pos = pos - tmp$mstart + 1)

                    # codon pos
                    if(mode == "codon"){
                      pos.df <- pos.df %>%
                        dplyr::mutate(pos = (pos - 1) %/% 3 + 1) %>%
                        fastplyr::f_group_by(rname, sample, sample_group,type, pos) %>%
                        fastplyr::f_summarise(rpm = sum(rpm), smooth = sum(smooth))
                    }
                  }
                  return(pos.df)
                }) -> tmpdf

                return(tmpdf)
              }) -> tmpdf2

              return(tmpdf2)
            }) -> pldf


            # ==================================================================
            # plot
            if(merge_rep == TRUE){
              pldf <- pldf %>%
                fastplyr::f_group_by(sample_group,type,rname,pos) %>%
                fastplyr::f_summarise(rpm = mean(rpm),
                                      sd = sd(smooth),
                                      smooth = mean(smooth)
                ) %>%
                dplyr::rename(sample = sample_group)

              shadow <- geom_ribbon(aes(ymin = smooth - sd,
                                        ymax = smooth + sd,
                                        x = pos,y = smooth,
                                        fill = type), alpha = var_alpha)
            }else{
              shadow <- NULL
            }
            # ================================================================
            tid <- tanno$idnew

            # loop for each transcript_id
            # x = 1
            lapply(seq_along(tid),function(x){
              tmp_feature <- subset(tanno, idnew == tid[x])
              tmp <- subset(pldf, rname == tid[x])

              tmp$sample <- paste(tmp$sample," (",tmp$type,")",sep = "")

              # whether add new signal range
              if(new_signal_range == TRUE){
                axis.text.y = element_blank()
                axis.ticks.y = element_blank()

                range <- paste("[0-",round(max(tmp$smooth),digits = range_digit),"]",sep = "")

                if (requireNamespace("ggpp", quietly = TRUE)) {
                  range_label <- ggpp::annotate(geom = "text_npc",
                                                npcx = range_x,npcy = range_y,
                                                label = range,
                                                size = range_size)
                } else {
                  warning("Package 'ggpp' is needed for this function to work.")
                }

              }else{
                axis.text.y = NULL
                axis.ticks.y = NULL
                range_label <- NULL
              }

              # structure
              if(mode == "codon"){
                struc_cds <- data.frame(x = 1,xend = tmp_feature$cds/3)

                sly2 <- ggside::geom_xsidesegment(data = struc_cds,
                                                  aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                  inherit.aes = F,
                                                  linewidth = as.numeric(cds_line[2]),color = cds_line[1])

                sly1 <- NULL
              }else{
                struc_exon <- data.frame(x = 1,xend = tmp_feature$exonlen)

                sly1 <- ggside::geom_xsidesegment(data = struc_exon,
                                                  aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                  inherit.aes = F,
                                                  linewidth = as.numeric(exon_line[2]),color = exon_line[1])

                struc_cds <- data.frame(x = tmp_feature$mstart,xend = tmp_feature$mstop)

                sly2 <- ggside::geom_xsidesegment(data = struc_cds,
                                                  aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                  inherit.aes = F,
                                                  linewidth = as.numeric(cds_line[2]),color = cds_line[1])
              }

              # xlabel
              if(mode == "codon"){
                xlabel <- "Ribosome position (codons / amino acids)"
              }else{
                xlabel <- "Ribosome position (nucleotides)"
              }

              # ================================================================
              p <-
                ggplot(tmp) +
                shadow +
                geom_line(aes(x = pos,y = smooth, color = type)) +
                range_label +
                facet +
                theme_bw() +
                theme(panel.grid = element_blank(),
                      axis.text = element_text(colour = "black"),
                      strip.text.y.left = element_text(angle = 0, hjust = 1),
                      strip.background = element_blank(),
                      strip.text = element_text(face = "bold"),
                      strip.placement = "outside",
                      axis.text.y = axis.text.y,
                      axis.ticks.y = axis.ticks.y,
                      ggh4x.facet.nestline = element_line(colour = "black"),
                      ggside.panel.background = element_blank(),
                      ggside.panel.border = element_blank()) +
                xlab(xlabel) +
                ylab("Ribosome density [AU]") +
                # add structure
                sly1 + sly2  +
                ggside::scale_xsidey_continuous(breaks = NULL) +
                ggside::ggside(collapse = "x")

              return(p)
            }) -> plist

            # combine plot list
            if (requireNamespace("cowplot", quietly = TRUE)) {
              cmb <- cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol)
            } else {
              warning("Package 'cowplot' is needed for this function to work.")
            }

            # return
            if(return_data == FALSE){
              return(cmb)
            }else{
              return(pldf)
            }
          }
)
