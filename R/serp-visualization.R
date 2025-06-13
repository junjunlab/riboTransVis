# ==============================================================================
# enrichment plot1 for single gene
# ==============================================================================

#' Plot ribosome profiling enrichment with smoothing and confidence intervals
#'
#' @description
#' `enrichment_plot1` creates a detailed visualization of ribosome profiling data, showing the enrichment of IP (immunoprecipitation) over total RNA across transcripts. It applies sliding window smoothing and calculates binomial confidence intervals using the Agresti–Coull method. Users can customize axis scaling, whether to focus on just the CDS region, and whether to display nucleotide (nt) or codon (amino acid) units on the x-axis.
#'
#' @param object An object of class `serp` containing:
#'   - `@total_occupancy`: Total counts per position
#'   - `@ip_occupancy`: IP counts per position
#'   - `@features`: Transcript feature information
#'   - `@library`: Mapping statistics
#'   - `@gene_name`: The gene of interest
#' @param shadow_alpha A numeric vector of length two specifying the transparency range for the confidence interval shading. Default is `c(0.1, 0.3)`.
#' @param window_size Integer. Size of the sliding window for smoothing. Default `45`.
#' @param retain_cds Logical. If `TRUE`, plot only the CDS region; if `FALSE`, include UTRs as well. Default `TRUE`.
#' @param log2_transform Logical. If `TRUE`, apply log2 scale to the y-axis (enrichment values). Default `FALSE`.
#' @param mode Character. Should be `"codon"` (default) to display x-axis in codons/amino acids, or `"nt"` to display in nucleotides.
#' @param new_signal_range Logical. Whether to add an annotation indicating the signal range (max smoothed ratio). Default `FALSE`.
#' @param enrich_threshold Numeric. A threshold line is drawn at this enrichment value. Default is `1`.
#' @param range_x Numeric. Position of the signal range label on the plot (only used if `new_signal_range = TRUE`).
#' @param range_y Numeric. Position of the signal range label on the plot (only used if `new_signal_range = TRUE`).
#' @param range_size Numeric. Font size for the signal range annotation label.
#' @param range_digit Integer. Number of digits used when displaying the maximum enrichment value.
#' @param exon_line A character vector specifying color and size for exon structure lines, e.g., `c("grey", 3)`.
#' @param cds_line A character vector specifying color and size for CDS structure lines, e.g., `c("grey30", 7)`.
#' @param sample_col Named character vector. Manual color mapping for samples. If `NULL`, default ggplot2 colors are used.
#' @param facet A ggplot2 facetting function, e.g., `facet_grid(sample ~ rname, switch = "y")`. Controls how plots are arranged.
#' @param nrow Integers. Number of rows and columns when combining multiple transcript plots.
#' @param ncol Integers. Number of rows and columns when combining multiple transcript plots.
#' @param return_data Logical. If `TRUE`, returns the smoothed enrichment calculation results as a dataframe instead of plotting. Default `FALSE`.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' - If `return_data = FALSE`, returns a `ggplot` or `cowplot` object containing the plotted results.
#' - If `return_data = TRUE`, returns a `data.frame` containing:
#'   * `sample`, `rname` - sample and transcript ID
#'   * `winstart` - start position of the sliding window
#'   * `win_ip`, `win_total` - counts in IP and Total within each window
#'   * `enrich`, `lower`, `upper` - mean enrichment and its 95% confidence interval
#'   * `bin` - window size
#'   * `nmft` - normalization factor (IP library size/Total library size)
#'
#' @details
#' The function first smooths the raw counts using a sliding window approach. For each window, it computes the proportion enriched (IP over IP+Total), and uses Agresti-Coull binomial confidence intervals via the `binom` package. The enrichment values are normalized using library sizes.
#'
#' Depending on the mode:
#' - `mode = "codon"`: x-axis position is rescaled by 3 to show codon units.
#' - `mode = "nt"`: x-axis position shows nucleotides.
#'
#' If `retain_cds = FALSE`, the full transcript including UTRs is plotted. Structural annotations visually mark exons and CDS regions separately.
#'
#' Users can opt to display enrichment on a log2 scale for better dynamic range visualization. Multiple samples and transcripts can be plotted together leveraging ggplot2 faceting.
#'
#' Requires the following R packages: `dplyr`, `purrr`, `ggplot2`, `ggside`, `cowplot`, `binom`, and optionally `ggpp`.
#'
#' @section Notes:
#' - Plots generated with `ggside` append a horizontal bar displaying transcript structure along the x-axis.
#' - If packages `binom`, `ggpp`, or `cowplot` are missing, a warning will be issued and some features may be disabled.
#'
#' @examples
#' \dontrun{
#' # Load a `serp` object
#' data(serp_example)
#'
#' # Basic plot, codon unit, focus on CDS
#' enrichment_plot1(serp_example)
#'
#' # Plot using nucleotide unit including UTRs
#' enrichment_plot1(serp_example, mode = "nt", retain_cds = FALSE)
#'
#' # Return smoothed data instead of plot
#' df <- enrichment_plot1(serp_example, return_data = TRUE)
#'
#' # Customize sample colors
#' enrichment_plot1(serp_example, sample_col = c(Sample1 = "blue", Sample2 = "red"))
#'
#' # Display enrichment on log2 scale
#' enrichment_plot1(serp_example, log2_transform = TRUE)
#' }
#'
#'
#'
#' @name enrichment_plot1
#' @rdname enrichment_plot1
#' @export
setGeneric("enrichment_plot1",function(object,...) standardGeneric("enrichment_plot1"))





#' @rdname enrichment_plot1
#' @export
setMethod("enrichment_plot1",
          signature(object = "serp"),
          function(object,
                   shadow_alpha = c(0.1,0.3),
                   window_size = 45,
                   retain_cds = TRUE,
                   log2_transform = FALSE,
                   mode = c("codon", "nt"),
                   new_signal_range = FALSE,
                   enrich_threshold = 1,
                   range_x = 0.99,
                   range_y = 0.98,
                   range_size = 4,
                   range_digit = 1,
                   exon_line = c("grey",3),
                   cds_line = c("grey30",7),
                   sample_col = NULL,
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
              dplyr::select(sample,sample_group,rname,pos,count) %>%
              dplyr::left_join(y = object@ip_occupancy[,c("sample","sample_group","rname","pos","count")],
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

                if(retain_cds == TRUE){
                  pos.df <- pos.df %>%
                    dplyr::filter(pos >= tmp$mstart & pos <= tmp$mstop) %>%
                    dplyr::mutate(pos = pos - tmp$mstart + 1)
                }

                # ==============================================================
                # smooth data for each position
                ip_total <- subset(object@library, sample == sp[s] & type == "ip")$mappped_reads
                total_total <- subset(object@library, sample == sp[s] & type == "total")$mappped_reads

                probnorm <- ip_total / total_total

                len_gene <- tmp$cds
                ip <- pos.df$count.y
                total <- pos.df$count.x

                win_starts <- 1:(len_gene - window_size + 1)
                win_ends <- window_size:len_gene

                cumsum_ip <- c(0, cumsum(ip))
                win_ip <- cumsum_ip[win_ends + 1] - cumsum_ip[win_starts]

                cumsum_total <- c(0, cumsum(total))
                win_total <- cumsum_total[win_ends+1] - cumsum_total[win_starts]

                if (requireNamespace("binom", quietly = TRUE)) {
                  cidf <- binom::binom.agresti.coull(win_ip, win_ip + win_total, conf.level = 0.95)
                } else {
                  warning("Package 'binom' is needed for this function to work.")
                }

                lower <- prob2odds(pmax(0, cidf$lower))
                lower <- ifelse(is.nan(lower), 0, lower)

                upper <- prob2odds(pmin(1, cidf$upper))
                upper <- ifelse(is.nan(upper), Inf, upper)

                mean <- prob2odds(cidf$mean)
                mean <- ifelse(is.nan(mean), 0, mean)

                # output
                df <- data.frame(sample = sp[s],
                                 rname = pos.df$rname[1],
                                 winstart = win_starts,
                                 win_ip = win_ip,
                                 win_total = win_total,
                                 enrich = mean/probnorm,
                                 lower = lower/probnorm,
                                 upper = upper/probnorm,
                                 bin = window_size,
                                 nmft = probnorm)

                return(df)
              }) -> smooth.df

              return(smooth.df)
            }) -> mb.smoothed

            # ==================================================================
            # plot
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
                if(retain_cds == TRUE){
                  struc_cds <- data.frame(x = 1,xend = tmp_feature$cds/3)

                  sly2 <- ggside::geom_xsidesegment(data = struc_cds,
                                                    aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                    inherit.aes = F,
                                                    linewidth = as.numeric(cds_line[2]),color = cds_line[1])

                  sly1 <- NULL
                }else{
                  stop('mode == "codon" is only used when retain_cds == TRUE!')
                }

              }else{
                if(retain_cds == TRUE){
                  struc_exon <- data.frame(x = 1,xend = tmp_feature$cds)

                  sly1 <- ggside::geom_xsidesegment(data = struc_exon,
                                                    aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                    inherit.aes = F,
                                                    linewidth = as.numeric(cds_line[2]),color = cds_line[1])

                  sly2 <- NULL
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


              }

              # xlabel
              if(mode == "codon"){
                xlabel <- "Ribosome position (codons / amino acids)"
              }else{
                xlabel <- "Ribosome position (nucleotides)"
              }

              # check mode
              if(mode == "codon"){
                if(retain_cds == TRUE){
                  layer.rect <- geom_rect(aes(xmin = (winstart + bin/2)/3 - 0.5,
                                              xmax = (winstart + bin/2)/3 + 0.5,
                                              ymin = lower/nmft,
                                              ymax = upper/nmft,
                                              alpha = 1/((1/win_total + 1/win_ip)*bin),
                                              fill = sample),
                                          show.legend = F)

                  layer.line <- geom_line(aes(x = (winstart + bin/2)/3, y = enrich/nmft, color = sample))
                }else{
                  stop('mode == "codon" is only used when retain_cds == TRUE!')
                }

              }else{
                layer.rect <- geom_rect(aes(xmin = (winstart + bin/2) - 0.5,
                                            xmax = (winstart + bin/2) + 0.5,
                                            ymin = lower/nmft,
                                            ymax = upper/nmft,
                                            alpha = 1/((1/win_total + 1/win_ip)*bin),
                                            fill = sample),
                                        show.legend = F)

                layer.line <- geom_line(aes(x = (winstart + bin/2), y = enrich/nmft, color = sample))
              }

              # y log2 scale
              if(log2_transform == TRUE){
                layer.y <- scale_y_continuous(trans = "log2",
                                              limits = 2^c(-5,sqrt(ceiling(max(tmp$enrich/tmp$nmft))) + 1),
                                              oob = scales::squish)
              }else{
                layer.y <- scale_y_continuous(trans = "identity", oob = scales::squish)
              }

              # sample color
              if(!is.null(sample_col)){
                scol <- scale_color_manual(values = sample_col)
                sfil <- scale_fill_manual(values = sample_col)
              }else{
                scol <- NULL
                sfil <- NULL
              }
              # ================================================================
              # draw

              p <-
                ggplot(tmp) +
                layer.rect + layer.line +
                geom_hline(yintercept = enrich_threshold,lty = "dashed",color = "black") +
                range_label +
                layer.y +
                scale_alpha_continuous(trans = "sqrt", name = "prec", range = shadow_alpha) +
                facet +
                scol + sfil +
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
                xlab(xlabel) +
                ylab("Mean enrichment (IP / total)[AU]") +
                # add structure
                sly1 + sly2 +
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
              return(mb.smoothed)
            }
          }
)




# ==============================================================================
# enrichment plot2 for single gene
# ==============================================================================

#' Plot smoothed enrichment ratio (IP/Total) for ribosome profiling
#'
#' @description
#' `enrichment_plot2` generates enrichment plots showing the ratio between IP and Total RNA signals across transcripts using rpm (Reads Per Million) normalized values. It optionally smooths the profile using a sliding window and can merge biological replicates by plotting the mean and standard deviation shading.
#'
#' @param object An object of class `serp` containing:
#'   - `@total_occupancy`: total RPM per position
#'   - `@ip_occupancy`: IP RPM per position
#'   - `@features`: transcript annotation (exons, CDS)
#'   - `@library`: library sizes and sample info
#'   - `@gene_name`: the gene to plot
#' @param merge_rep Logical. If `TRUE`, merge biological replicates and plot mean ± SD shading. Default `FALSE`.
#' @param var_alpha Numeric. Transparency value for SD error ribbon if `merge_rep = TRUE`.
#' @param smooth Logical. Whether to apply rolling mean smoothing (`zoo::rollsum`) to RPM profiles. Default `TRUE`.
#' @param window_size Integer. Size of the rolling window for smoothing. Default `45`.
#' @param retain_cds Logical. If `TRUE`, restrict analysis to CDS region only; otherwise include UTRs.
#' @param mode Character. Units on x-axis: `"codon"` for codon positions or `"nt"` for nucleotide positions.
#' @param new_signal_range Logical. Whether to annotate the maximum signal range on the plot.
#' @param enrich_threshold Numeric. Threshold line to highlight enrichment level, default `1`.
#' @param range_x Numeric. Positions (npc units) to place the signal range annotation if enabled.
#' @param range_y Numeric. Positions (npc units) to place the signal range annotation if enabled.
#' @param range_size Numeric. Font size of the signal range annotation.
#' @param range_digit Integer. Number of digits to display in the signal range maximum.
#' @param exon_line A character vector specifying color and thickness for exon structure annotations.
#' @param cds_line A character vector specifying color and thickness for CDS structure annotations.
#' @param facet A `ggplot2` faceting function, e.g., `facet_grid(sample ~ rname, switch = "y")`.
#' @param nrow Integers. Number of rows or columns if plotting multiple transcripts at once.
#' @param ncol Integers. Number of rows or columns if plotting multiple transcripts at once.
#' @param return_data Logical. If `TRUE`, return processed data instead of a plot.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' - If `return_data = FALSE`, returns a `ggplot` (possibly combined with `cowplot::plot_grid`) showing the enrichment ratio profiles.
#' - If `return_data = TRUE`, returns a `data.frame` containing:
#'   * `rname`, `sample`, `pos` - position, sample, transcript ID
#'   * `rpm.x` (total) and `rpm.y` (IP) RPM values
#'   * `sm1`, `sm2` - smoothed values (if smoothing is enabled)
#'   * `smratio` - smoothed ratio of IP/Total
#'   * `sd` - standard deviation of ratio (if merging replicates)
#'
#' @details
#' - The function calculates per-position rpm ratio (IP/Total) and optionally smooths using a centered rolling sum (sliding window) with size `window_size`.
#' - When `merge_rep = TRUE`, replicates are averaged by `sample_group`, and an error ribbon ±sd is drawn for variability visualization.
#' - Supports two modes of display:
#'   * `"codon"`: collapsing 3-nt to codon units (only for CDS retained region).
#'   * `"nt"`: nucleotide units.
#'
#' - Multiple transcripts can be plotted together using ggplot2 facetting.
#' - Annotations (exons and CDS regions) are displayed using `ggside` panel extensions.
#'
#' Required packages: `dplyr`, `purrr`, `ggplot2`, `ggside`, `cowplot`, `zoo`, and optionally `ggpp` if signal range labeling is enabled.
#'
#' @section Notes:
#' - Smoothing uses center aligning for better positional accuracy (`align = "center"` in rolling sum).
#' - If the `zoo` package is missing when `smooth = TRUE`, a warning will be issued.
#' - Replicate variance visualization is only available when `merge_rep = TRUE`.
#'
#' @examples
#' \dontrun{
#' # Given a serp object
#' data(serp_example)
#'
#' # Basic plot with smoothing
#' enrichment_plot2(serp_example)
#'
#' # Plot without smoothing
#' enrichment_plot2(serp_example, smooth = FALSE)
#'
#' # Merge replicates and plot mean ± SD
#' enrichment_plot2(serp_example, merge_rep = TRUE)
#'
#' # Plot nucleotide positions (including UTRs)
#' enrichment_plot2(serp_example, mode = "nt", retain_cds = FALSE)
#'
#' # Return processed data
#' df <- enrichment_plot2(serp_example, return_data = TRUE)
#' }
#'
#' @importFrom stats sd
#'
#'
#' @name enrichment_plot2
#' @rdname enrichment_plot2
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
                if(retain_cds == TRUE){
                  struc_cds <- data.frame(x = 1,xend = tmp_feature$cds/3)

                  sly2 <- ggside::geom_xsidesegment(data = struc_cds,
                                                    aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                    inherit.aes = F,
                                                    linewidth = as.numeric(cds_line[2]),color = cds_line[1])

                  sly1 <- NULL
                }else{
                  stop('mode == "codon" is only used when retain_cds == TRUE!')
                }

              }else{
                if(retain_cds == TRUE){
                  struc_exon <- data.frame(x = 1,xend = tmp_feature$cds)

                  sly1 <- ggside::geom_xsidesegment(data = struc_exon,
                                                    aes(x = x,xend = xend,y = 0.5,yend = 0.5),
                                                    inherit.aes = F,
                                                    linewidth = as.numeric(cds_line[2]),color = cds_line[1])

                  sly2 <- NULL
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
#' @param ... Additional arguments (currently unused).
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
#' trans_plot2(serp_obj, mode = "codon")
#'
#' # Merge replicates and plot with standard deviation shading
#' trans_plot2(serp_obj, merge_rep = TRUE)
#'
#' # Return processed data without plotting
#' data <- trans_plot2(serp_obj, return_data = TRUE)
#' }
#' @export
setGeneric("trans_plot2",function(object,...) standardGeneric("trans_plot2"))



#' @rdname trans_plot2
#' @export
setMethod("trans_plot2",
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
