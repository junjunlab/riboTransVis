# ==============================================================================
# function for metagene plot
# ==============================================================================

#' Generate metagene plots for ribosome profiling data
#'
#' @description
#' Creates a metagene plot showing the average ribosome occupancy relative to start or stop codons
#' from ribosome profiling data. This function allows for various customizations including data filtering,
#' statistical methods, and visualization options.
#'
#' @param object A ribotrans object containing ribosome profiling data
#' @param merge_rep Logical; whether to merge replicates. Default is FALSE.
#' @param var_alpha Numeric; transparency level for variance ribbon (standard deviation or confidence intervals). Default is 0.5.
#' @param selected_genes Character vector; subset of genes to analyze. Default is NULL (all genes).
#' @param do_offset_correct Logical; whether to apply offset correction to read positions. Default is FALSE.
#' @param position_shift Numeric; number of nucleotides to shift read positions. Default is 0.
#' @param method Character; method to aggregate data ("median", "mean", or "sum"). Default is "median".
#' @param do_bootstrap Logical; whether to perform bootstrap resampling. Default is FALSE.
#' @param boot_n Numeric; number of bootstrap iterations. Default is 1000.
#' @param conf Numeric; confidence level for bootstrap confidence intervals. Default is 0.95.
#' @param min_cds_length Numeric; minimum CDS length to include in analysis. Default is 600.
#' @param min_counts Numeric; minimum read count threshold for genes. Default is 64.
#' @param exclude_length Numeric vector; regions to exclude near start/stop codons, c(start_region, stop_region). Default is c(90,90).
#' @param type Character; reference point for alignment ("rel2start" or "rel2stop"). Default is "rel2start".
#' @param return_data Logical; if TRUE, returns the data frame instead of the plot. Default is FALSE.
#' @param mode Character; count mode ("nt" for nucleotide-level or "codon" for codon-level aggregation). Default is "nt".
#' @param read_length Numeric vector; range of read lengths to include c(min_length, max_length). Default is c(25,35).
#' @param rel2st_dist Numeric vector; region to plot around start codon c(upstream, downstream). Default is c(-50,100).
#' @param rel2sp_dist Numeric vector; region to plot around stop codon c(upstream, downstream). Default is c(-100,50).
#' @param facet_wrap A facet_wrap specification for ggplot2. Default is ggplot2::facet_wrap(~sample).
#' @param ... Additional arguments (currently unused).
#'
#'
#' @return If return_data is FALSE (default), returns a ggplot object. If return_data is TRUE, returns a data frame with the processed data.
#'
#' @details
#' This function analyzes ribosome density patterns around start or stop codons by normalizing
#' read counts and plotting the average profile. It supports various filtering criteria to focus on
#' specific subsets of the data.
#'
#' When merge_rep=TRUE, replicates are aggregated by their sample_group and displayed with standard deviation.
#' When do_bootstrap=TRUE, bootstrap confidence intervals are calculated and displayed.
#'
#' The function normalizes raw counts by the average count per nucleotide for each transcript before aggregation.
#' The final plot shows relative ribosome occupancy across the specified region.
#'
#' @note
#' - merge_rep=TRUE and do_bootstrap=TRUE cannot be used together.
#' - For codon-level analysis (mode="codon"), positions are converted from nucleotides to codons.
#'
#' @examples
#' \dontrun{
#' # Basic metagene plot around start codons
#' metagene_plot(ribo_obj)
#'
#' # Customize plot with bootstrap confidence intervals
#' metagene_plot(ribo_obj, do_bootstrap = TRUE,
#'               boot_n = 500, conf = 0.9)
#'
#' # Plot around stop codons with merged replicates
#' metagene_plot(ribo_obj, type = "rel2stop",
#'               merge_rep = TRUE, var_alpha = 0.3)
#'
#' # Analyze specific genes at codon resolution
#' metagene_plot(ribo_obj, mode = "codon",
#'               selected_genes = c("YAL001C", "YAL002W"))
#'
#' # Return data instead of plot for custom visualization
#' plot_data <- metagene_plot(ribo_obj, return_data = TRUE)
#' }
#'
#' @importFrom ggplot2 ggplot geom_path geom_ribbon aes theme_bw theme element_blank element_text xlab ylab
#' @importFrom dplyr mutate filter inner_join summarise group_by rename
#' @importFrom stats median sd
#' @importFrom fastplyr f_filter f_group_by f_summarise
#' @importFrom data.table `:=`
#'
#' @export
setGeneric("metagene_plot",function(object,...) standardGeneric("metagene_plot"))






#' @rdname metagene_plot
#' @export
setMethod("metagene_plot",
          signature(object = "ribotrans"),
          function(object,
                   merge_rep = FALSE,
                   var_alpha = 0.5,
                   selected_genes = NULL,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   method = c("median", "mean","sum"),
                   do_bootstrap = FALSE,
                   boot_n = 1000,
                   conf = 0.95,
                   min_cds_length = 600,
                   min_counts = 64,
                   exclude_length = c(90,90),
                   type = c("rel2start","rel2stop"),
                   return_data = FALSE,
                   mode = c("nt", "codon"),
                   read_length = c(25,35),
                   rel2st_dist = c(-50,100),
                   rel2sp_dist = c(-100,50),
                   facet_wrap = ggplot2::facet_wrap(~sample)){
            method <- match.arg(method, choices = c("median", "mean","sum"))
            type <- match.arg(type, choices = c("rel2start","rel2stop"))
            mode <- match.arg(mode, choices = c("nt", "codon"))


            if(merge_rep == TRUE & do_bootstrap == TRUE){
              stop("merge_rep is useless when do_bootstrap=TRUE!")
            }
            # ==================================================================
            # check type
            if(type == "rel2start"){
              dist <- rel2st_dist
              var <- "mstart"
            }else{
              dist <- rel2sp_dist
              var <- "mstop"
            }

            features <- object@features

            # whether do reads offset correction
            if(do_offset_correct == TRUE){
              sry <- do_offset_correction(object = object,shift = position_shift)
            }else{
              sry <- object@summary_info
            }

            sry <- data.table::as.data.table(sry)


            sry <- sry[(mstart != 0 | mstop != 0) &
                         qwidth >= read_length[1] & qwidth <= read_length[2]]

            # filter genes
            if(!is.null(selected_genes)){
              ft <- subset(object@features, gene %in% selected_genes)

              sry <- subset(sry, rname %in% ft$idnew)
            }

            # average counts per position
            avg.ct <- sry %>%
              dplyr::mutate(cdslen = mstop - mstart + 1) %>%
              dplyr::mutate(relst = pos - mstart, relsp = pos - mstop) %>%
              fastplyr::f_filter(relst > exclude_length[1] & relsp < -exclude_length[2]) %>%
              fastplyr::f_group_by(sample,rname,cdslen) %>%
              fastplyr::f_summarise(counts = sum(count)) %>%
              fastplyr::f_filter(counts > min_counts) %>%
              dplyr::mutate(avg_ct = counts/cdslen)


            # filter relative distance
            if (requireNamespace("dtplyr", quietly = TRUE)) {
              dt_sry <- dtplyr::lazy_dt(sry)
              dt_avgct <- dtplyr::lazy_dt(avg.ct)
            } else {
              warning("Package 'dtplyr' is needed for this function to work.")
            }


            pltdf.tmp <- dt_sry %>%
              dplyr::filter(mstop - mstart >= min_cds_length) %>%
              dplyr::inner_join(dt_avgct, by = c("sample","rname")) %>%
              dplyr::mutate(rel  = pos - .data[[var]],
                            norm = count / avg_ct) %>%
              dplyr::filter(rel >= dist[1], rel <= dist[2]) %>%
              data.table::as.data.table()


            # check show mode
            if(mode == "codon"){
              pltdf.tmp[, rel := rel %/% 3 + 1]
              pltdf.tmp <- pltdf.tmp[
                , .(norm = sum(norm)),
                by = .(sample, sample_group, rname, rel)
              ]

              flen <- (abs(dist[2]) + abs(dist[1]) + 1) %/% 3 + 1
            }else{
              flen <- abs(dist[2]) + abs(dist[1]) + 1
            }


            # check whether do bootstraps
            if(do_bootstrap == TRUE){
              # check method
              sp <- unique(pltdf.tmp$sample)

              # x = 1
              lapply(seq_along(sp),function(x){
                tmp <- subset(pltdf.tmp, sample == sp[x])

                if (requireNamespace(c("tidyr","tibble"), quietly = TRUE)) {
                  mat <- tmp %>%
                    fastplyr::f_group_by(rname,rel) %>%
                    fastplyr::f_summarise(norm = sum(norm)) %>%
                    tidyr::pivot_wider(names_from = rel, values_from = norm, values_fill = NA) %>%
                    tibble::column_to_rownames(var = "rname")
                } else {
                  warning("Package 'tidyr' and 'tibble' are needed for this function to work.")
                }


                mat <- mat[,as.character(sort(unique(tmp$rel)))]

                if(method == "median"){
                  stat.df <- apply(mat, 2, stats::median, na.rm = T)
                }else if(method == "mean"){
                  stat.df <- apply(mat, 2, mean, na.rm = T)
                }else{
                  stat.df <- apply(mat, 2, sum, na.rm = T)
                }

                if (requireNamespace("tibble", quietly = TRUE)) {
                  stat.df <- data.frame(stat.df) %>%
                    tibble::rownames_to_column(var = "rel") %>%
                    dplyr::mutate(rel = as.numeric(rel))
                } else {
                  warning("Package 'tibble' is needed for this function to work.")
                }

                colnames(stat.df) <- c("rel","normsm")

                boot.df <- do_boot(mat = mat, boot_n = boot_n, method = method, conf = conf)

                boot.df <- boot.df %>%
                  dplyr::left_join(y = stat.df,by = "rel")

                boot.df$sample <- tmp$sample[1]
                boot.df$sample_group <- tmp$sample_group[1]

                return(boot.df)
              }) %>% do.call("rbind", .) %>% data.frame() -> pltdf
            }else{
              pltdf.tmp <- pltdf.tmp %>%
                fastplyr::f_group_by(sample,sample_group,rel)

              # check methods
              if(method == "median"){
                pltdf <- pltdf.tmp %>%
                  fastplyr::f_summarise(normsm = median(norm, na.rm = T))
              }else if(method == "mean"){
                pltdf <- pltdf.tmp %>%
                  fastplyr::f_summarise(normsm = mean(norm, na.rm = T))
              }else{
                pltdf <- pltdf.tmp %>%
                  fastplyr::f_summarise(normsm = sum(norm, na.rm = T))
              }
            }


            # average occupancy
            sm <- pltdf %>%
              fastplyr::f_group_by(sample) %>%
              fastplyr::f_summarise(norm_avg = sum(normsm)/flen)

            pltdf <- pltdf %>%
              dplyr::inner_join(y = sm,by = "sample")

            if(do_bootstrap == TRUE){
              pltdf <- pltdf %>%
                dplyr::mutate(avg = normsm/norm_avg,
                              ci_low = ci_low/norm_avg,
                              ci_high = ci_high/norm_avg)
            }else{
              pltdf <- pltdf %>%
                dplyr::mutate(avg = normsm/norm_avg)
            }


            # whether aggregate replicates
            if(merge_rep == TRUE){
              if(do_bootstrap == FALSE){
                pltdf <- pltdf %>%
                  fastplyr::f_group_by(sample_group,rel) %>%
                  fastplyr::f_summarise(avg = mean(avg),sd = stats::sd(avg)) %>%
                  dplyr::rename(sample = sample_group)

                shadow <- geom_ribbon(aes(ymin = avg - sd,
                                          ymax = avg + sd,
                                          x = rel,y = avg,
                                          fill = sample), alpha = var_alpha,
                                      show.legend = F)
              }

            }else{
              if(do_bootstrap == TRUE){
                pltdf <- pltdf

                shadow <- geom_ribbon(aes(ymin = ci_low,
                                          ymax = ci_high,
                                          x = rel,y = avg,
                                          fill = sample), alpha = var_alpha,
                                      show.legend = F)
              }else{
                shadow <- NULL
              }

            }

            # check type
            if(type == "rel2start"){
              if(mode == "codon"){
                xlb <- "Distance to start codon (codon)"
              }else{
                xlb <- "Distance to start codon (nt)"
              }
            }else{
              if(mode == "codon"){
                xlb <- "Distance to stop codon (codon)"
              }else{
                xlb <- "Distance to stop codon (nt)"
              }
            }
            # ==================================================================
            # plot
            p <-
              ggplot(pltdf) +
              shadow +
              geom_path(aes(x = rel,y = avg, color = sample)) +
              theme_bw() +
              theme(panel.grid = element_blank(),
                    strip.text = element_text(colour = "black",face = "bold",size = rel(1)),
                    axis.text = element_text(colour = "black")) +
              # facet_wrap(~sample) +
              facet_wrap +
              xlab(xlb) +
              ylab("Average ribosome occupancy")

            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(pltdf)
            }
          }
)

