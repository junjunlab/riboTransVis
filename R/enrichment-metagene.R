
#' Create enrichment metagene plots for ribosome profiling data
#'
#' @description
#' This function generates metagene plots showing the enrichment of ribosome footprints
#' around start or stop codons. It calculates the ratio of IP to total reads for each position
#' relative to start/stop codons, and can display the average profile across genes with various
#' statistical approaches.
#'
#' @param object A \code{ribotrans} object containing ribosome profiling data.
#' @param merge_rep Logical. If \code{TRUE}, replicates will be merged by sample group. Default is \code{FALSE}.
#' @param var_alpha Numeric value between 0 and 1 indicating the transparency of the confidence interval
#'   ribbons. Default is 0.5.
#' @param selected_genes Character vector of gene names to include in the analysis. If \code{NULL} (default),
#'   all genes will be used.
#' @param do_offset_correct Logical. If \code{TRUE}, apply offset correction to read positions. Default is \code{FALSE}.
#' @param position_shift Integer indicating the number of nucleotides to shift positions. Default is 0.
#' @param method Character string specifying the method for aggregating values across genes.
#'   Must be one of "median" (default), "mean", or "sum".
#' @param do_bootstrap Logical. If \code{TRUE}, perform bootstrap analysis to estimate confidence intervals.
#'   Default is \code{FALSE}.
#' @param boot_n Integer specifying the number of bootstrap iterations. Default is 1000.
#' @param conf Numeric value between 0 and 1 specifying the confidence level for bootstrap intervals.
#'   Default is 0.95.
#' @param min_cds_length Integer specifying the minimum CDS length in nucleotides to include in analysis.
#'   Default is 600.
#' @param min_counts Integer specifying the minimum number of counts per gene to include in analysis.
#'   Default is 64.
#' @param exclude_length Numeric vector of length 2 specifying regions to exclude from start 0 and end 1
#'   of CDS. Default is c(90, 90).
#' @param type Character string specifying the reference point for alignment.
#'   Must be one of "rel2start" (default) or "rel2stop".
#' @param return_data Logical. If \code{TRUE}, return the plot data instead of the plot. Default is \code{FALSE}.
#' @param mode Character string specifying the x-axis units.
#'   Must be one of "nt" (nucleotides, default) or "codon".
#' @param read_length Numeric vector of length 2 specifying the range of read lengths to include "min, max".
#'   Default is c(25, 35).
#' @param rel2st_dist Numeric vector of length 2 specifying the distance range relative to start codon "min, max".
#'   Default is c(-50, 100).
#' @param rel2sp_dist Numeric vector of length 2 specifying the distance range relative to stop codon "min, max".
#'   Default is c(-100, 50).
#' @param facet_wrap A facet specification for \code{ggplot2} to control faceting of the plot.
#' @param ... Additional arguments (currently unused).
#'   Default is \code{ggplot2::facet_wrap(~sample)}.
#'
#' @return If \code{return_data = FALSE} (default), returns a \code{ggplot} object of the metagene plot.
#'   If \code{return_data = TRUE}, returns a data frame containing the processed data used for plotting.
#'
#' @details
#' The function calculates the enrichment of ribosome footprints around start or stop codons by
#' dividing IP signals by the corresponding total RNA signals at each position. It can handle
#' multiple experimental conditions and replicates, and provides various options to customize the
#' visualization and statistical analysis.
#'
#' The enrichment values are normalized to facilitate comparison between samples. The normalization
#' is performed by dividing each position's enrichment by the average enrichment across all positions
#' in the analyzed range.
#'
#' When \code{do_bootstrap = TRUE}, confidence intervals are estimated using bootstrap resampling
#' of genes. When \code{merge_rep = TRUE}, data from replicates are averaged and standard deviation
#' is shown.
#'
#' @note
#' The parameters \code{merge_rep} and \code{do_bootstrap} cannot both be \code{TRUE} simultaneously.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' enrichment_metagene_plot(ribo_obj)
#'
#' # Plot enrichment around stop codons with bootstrapping
#' enrichment_metagene_plot(ribo_obj,
#'                         type = "rel2stop",
#'                         do_bootstrap = TRUE,
#'                         boot_n = 500)
#'
#' # Plot enrichment for specific genes, return data instead of plot
#' gene_data <- enrichment_metagene_plot(ribo_obj,
#'                                      selected_genes = c("ACTB", "GAPDH"),
#'                                      return_data = TRUE)
#'
#' # Plot with custom parameters and merge replicates
#' enrichment_metagene_plot(ribo_obj,
#'                         merge_rep = TRUE,
#'                         method = "mean",
#'                         min_cds_length = 1000,
#'                         min_counts = 100,
#'                         mode = "codon")
#'
#' # Plot with custom faceting
#' enrichment_metagene_plot(ribo_obj,
#'                         facet_wrap = ggplot2::facet_wrap(~sample_group, ncol = 2))
#' }
#'
#' @seealso
#' \code{\link{do_offset_correction}} for offset correction methods
#'
#' @importFrom data.table as.data.table tstrsplit
#' @importFrom dplyr mutate filter inner_join left_join group_by summarise rename
#' @importFrom fastplyr f_filter f_group_by f_summarise f_inner_join
#' @importFrom ggplot2 ggplot geom_path geom_ribbon geom_hline theme_bw theme element_blank element_text facet_wrap xlab ylab aes
#' @importFrom stats median sd
#'
#' @export
setGeneric("enrichment_metagene_plot",function(object,...) standardGeneric("enrichment_metagene_plot"))




#' @rdname enrichment_metagene_plot
#' @export
setMethod("enrichment_metagene_plot",
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

            sry <- as.data.table(sry)

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
              pltdf.tmp <- pltdf.tmp[, .(norm = sum(norm)),
                by = .(sample, sample_group, rname, rel)]

              flen <- (abs(dist[2]) + abs(dist[1]) + 1) %/% 3 + 1
            }else{
              flen <- abs(dist[2]) + abs(dist[1]) + 1
            }


            # check whether do bootstraps
            if(do_bootstrap == TRUE){
              # check method
              sp <- unique(pltdf.tmp$sample)

              lib <- object@library
              tt_sp <- subset(lib, type == "total")[,c("type", "sample", "sample_group")]

              ip_sp <- subset(lib, type == "ip")[,c("type", "sample", "sample_group")]

              sp_mr <- tt_sp %>% dplyr::left_join(y = ip_sp,by = "sample")

              sp.tt <- paste(sp_mr$sample,sp_mr$type.x,sep = "-")
              sp.ip <- paste(sp_mr$sample,sp_mr$type.y,sep = "-")

              # x = 1
              lapply(seq_along(sp.tt),function(x){
                tmp.tt <- subset(pltdf.tmp, sample == sp.tt[x])
                tmp.ip <- subset(pltdf.tmp, sample == sp.ip[x])

                if (requireNamespace(c("tidyr","tibble"), quietly = TRUE)) {
                  mat.tt <- tmp.tt %>%
                    fastplyr::f_group_by(rname,rel) %>%
                    fastplyr::f_summarise(norm = sum(norm)) %>%
                    tidyr::pivot_wider(names_from = rel, values_from = norm, values_fill = NA) %>%
                    tibble::column_to_rownames(var = "rname")
                } else {
                  warning("Package 'tidyr' and 'tibble' are needed for this function to work.")
                }

                mat.tt <- mat.tt[,as.character(sort(unique(tmp.tt$rel)))]

                if (requireNamespace(c("tidyr","tibble"), quietly = TRUE)) {
                  mat.ip <- tmp.ip %>%
                    fastplyr::f_group_by(rname,rel) %>%
                    fastplyr::f_summarise(norm = sum(norm)) %>%
                    tidyr::pivot_wider(names_from = rel, values_from = norm, values_fill = NA) %>%
                    tibble::column_to_rownames(var = "rname")
                } else {
                  warning("Package 'tidyr' and 'tibble' are needed for this function to work.")
                }

                mat.ip <- mat.ip[,as.character(sort(unique(tmp.ip$rel)))]

                ids <- intersect(rownames(mat.tt),rownames(mat.ip))

                mat.tt <- mat.tt[ids,]
                mat.ip <- mat.ip[ids,]

                # do ratio
                enrich.mat <- mat.ip / mat.tt

                if(method == "median"){
                  stat.df <- apply(enrich.mat, 2, stats::median, na.rm = T)
                }else if(method == "mean"){
                  stat.df <- apply(enrich.mat, 2, mean, na.rm = T)
                }else{
                  stat.df <- apply(enrich.mat, 2, sum, na.rm = T)
                }

                if (requireNamespace("tibble", quietly = TRUE)) {
                  stat.df <- data.frame(stat.df) %>%
                    tibble::rownames_to_column(var = "rel") %>%
                    dplyr::mutate(rel = as.numeric(rel))

                } else {
                  warning("Package 'tibble' is needed for this function to work.")
                }

                colnames(stat.df) <- c("rel","normsm")

                boot.df <- do_boot(mat = enrich.mat, boot_n = boot_n, method = method, conf = conf)

                boot.df <- boot.df %>%
                  dplyr::left_join(y = stat.df,by = "rel")

                boot.df$sample <- sp_mr$sample[x]
                boot.df$sample_group <- sp_mr$sample_group.x[x]

                return(boot.df)
              }) %>% do.call("rbind", .) %>% data.frame() -> pltdf
            }else{
              pltdf.tmp.tt <- pltdf.tmp %>%
                fastplyr::f_filter(sample %in% sp.tt)

              # re-name total samples
              pltdf.tmp.tt[, sample := data.table::tstrsplit(sample, "-total", fixed = TRUE)[[1]]]
              pltdf.tmp.tt[, sample_group := data.table::tstrsplit(sample_group, "-total", fixed = TRUE)[[1]]]

              # re-name ip samples
              pltdf.tmp.ip <- pltdf.tmp %>%
                fastplyr::f_filter(sample %in% sp.ip)

              pltdf.tmp.ip[, sample := data.table::tstrsplit(sample, "-ip", fixed = TRUE)[[1]]]
              pltdf.tmp.ip[, sample_group := data.table::tstrsplit(sample_group, "-ip", fixed = TRUE)[[1]]]

              # merge total and ip
              pltdf.tmp.mer <- pltdf.tmp.tt %>%
                fastplyr::f_inner_join(y = pltdf.tmp.ip,by = c("sample","sample_group","rname","rel")) %>%
                dplyr::mutate(enrich = norm.y / norm.x) %>%
                fastplyr::f_group_by(sample,sample_group,rel)

              # check methods
              if(method == "median"){
                pltdf <- pltdf.tmp.mer %>%
                  fastplyr::f_summarise(normsm = median(enrich))
              }else if(method == "mean"){
                pltdf <- pltdf.tmp.mer %>%
                  fastplyr::f_summarise(normsm = mean(enrich))
              }else{
                pltdf <- pltdf.tmp.mer %>%
                  fastplyr::f_summarise(normsm = sum(enrich))
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
              geom_hline(yintercept = 1,lty = "dashed", color = "black") +
              theme_bw() +
              theme(panel.grid = element_blank(),
                    strip.text = element_text(colour = "black",face = "bold",size = rel(1)),
                    axis.text = element_text(colour = "black")) +
              # facet_wrap(~sample) +
              facet_wrap +
              xlab(xlb) +
              ylab("Normalized enrichment (ip/total)[AU]")

            # return
            if(return_data == FALSE){
              return(p)
            }else{
              return(pltdf)
            }
          }
)
