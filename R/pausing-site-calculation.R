
#' @title Identify Ribosome Pausing Sites from Ribosome Profiling Data
#'
#' @description
#' This function identifies ribosome pausing sites by calculating pause scores
#' based on the PausePred algorithm. It processes ribosome profiling data to
#' detect positions where ribosomes pause during translation elongation.
#'
#' @param object A RiboSeq object containing ribosome profiling data with
#'   summary_info and features slots
#' @param sample_selected Character vector. Names of samples to include in the
#'   analysis. If NULL (default), all samples will be processed
#' @param merge_rep Logical. Whether to merge technical replicates by taking
#'   the mean of normalized values within sample groups. Default is FALSE
#' @param do_offset_correct Logical. Whether to perform read offset correction
#'   to account for ribosome A-site positioning. Default is FALSE
#' @param position_shift Numeric. The number of nucleotides to shift reads
#'   when offset correction is applied. Only used when do_offset_correct = TRUE.
#'   Default is 0
#' @param exclude_length Numeric vector of length 2. Number of nucleotides to
#'   exclude from the start and end of coding sequences.
#'   Default is c(0,0)
#' @param min_counts Numeric. Minimum total read counts required for a gene
#'   to be included in the analysis. Genes with fewer counts are filtered out.
#'   Default is 64
#' @param window Numeric. Window size (in nucleotides) used for calculating
#'   pause scores. The algorithm uses overlapping windows to assess local
#'   ribosome density. Default is 200
#' @param ... Additional arguments (currently unused).
#'
#' @return A data.table containing pausing site information with the following columns:
#' \describe{
#'   \item{sample}{Sample identifier}
#'   \item{rname}{Gene/transcript identifier}
#'   \item{relst}{Relative position from CDS start}
#'   \item{normsm}{Normalized ribosome density (read count / average count)}
#'   \item{sum_1}{Sum of normalized densities in window pos, pos + window}
#'   \item{sum_2}{Sum of normalized densities in window pos + window/2, pos + 1.5*window}
#'   \item{pause_score}{Calculated pause score using PausePred algorithm}
#'   \item{total_sum}{Total sum of sum_1 and sum_2}
#'   \item{min_sum}{Minimum of sum_1 and sum_2}
#'   \item{max_sum}{Maximum of sum_1 and sum_2}
#'   \item{balance_score}{Balance score (min_sum / total_sum) for quality control}
#'   \item{imbalance_ratio}{Ratio of max_sum to min_sum for quality assessment}
#'   \item{z_score}{Z-score of pause score relative to the gene}
#' }
#'
#' @details
#' The function implements the PausePred algorithm for identifying ribosome pausing
#' sites with the following steps:
#'
#' \enumerate{
#'   \item \strong{Data preprocessing}: Applies optional offset correction and filters
#'         samples and genomic regions based on user specifications
#'   \item \strong{Normalization}: Calculates normalized ribosome density by dividing
#'         read counts by the average count per position for each gene
#'   \item \strong{Quality filtering}: Excludes genes with low total read counts
#'         (< min_counts) to ensure reliable statistical inference
#'   \item \strong{Window-based scoring}: For each position, calculates sums in two
#'         overlapping windows and computes pause scores using the formula:
#'         \deqn{pause\_score = \frac{window}{2} \times normsm \times \frac{sum_1 + sum_2}{sum_1 \times sum_2}}
#'   \item \strong{Quality metrics}: Adds balance and imbalance ratio metrics for
#'         downstream quality control and filtering
#' }
#'
#' The pause score reflects the likelihood that a given position represents a
#' genuine pausing site, with higher scores indicating stronger evidence for
#' ribosome pausing. The algorithm is designed to identify positions where
#' ribosome density is elevated relative to neighboring regions.
#'
#' \strong{Window Strategy:}
#' \itemize{
#'   \item Window 1: position, position + window - captures downstream density
#'   \item Window 2: position + window/2, position + 1.5*window - overlapping window
#'   \item The overlapping design helps distinguish genuine pauses from noise
#' }
#'
#' \strong{Quality Control Metrics:}
#' \itemize{
#'   \item \code{balance_score}: Ranges from 0 to 0.5, with higher values indicating
#'         more balanced window sums (recommended threshold: > 0.1-0.2)
#'   \item \code{imbalance_ratio}: Ratio of larger to smaller window sum
#'         (recommended threshold: < 5-10)
#' }
#'
#' @section Filtering Recommendations:
#' For downstream analysis, consider filtering pause sites based on:
#' \itemize{
#'   \item Pause score threshold (typically 10-20 for initial screening)
#'   \item Balance score > 0.15 to exclude sites with extremely unbalanced windows
#'   \item Imbalance ratio < 5 to remove artifacts from coverage imbalances
#'   \item Minimum total coverage to ensure statistical reliability
#' }
#'
#' @seealso
#' \code{\link{do_offset_correction}} for read offset correction
#'
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' pause_sites <- get_pausing_sites(ribo_object)
#'
#' # Analyze specific samples with offset correction
#' pause_sites <- get_pausing_sites(
#'   object = ribo_object,
#'   sample_selected = c("sample1", "sample2"),
#'   do_offset_correct = TRUE,
#'   position_shift = 12
#' )
#'
#' # Merge replicates and use stricter filtering
#' pause_sites <- get_pausing_sites(
#'   object = ribo_object,
#'   merge_rep = TRUE,
#'   exclude_length = c(30, 30),  # Exclude first/last 30 nt
#'   min_counts = 100,            # Higher count threshold
#'   window = 150                 # Smaller window size
#' )
#'
#' # Filter high-quality pause sites
#' high_quality_pauses <- pause_sites[
#'   pause_score >= 20 &
#'   balance_score >= 0.15 &
#'   imbalance_ratio <= 5
#' ]
#'
#' # Summary statistics
#' summary(pause_sites$pause_score)
#' table(pause_sites$sample)
#' }
#'
#'
#'
#' @export
#' @importFrom dplyr mutate inner_join select rename filter full_join
#' @importFrom fastplyr f_filter f_select f_group_by f_summarise
#' @importFrom data.table as.data.table fifelse
setGeneric("get_pausing_sites",function(object,...) standardGeneric("get_pausing_sites"))




#' @rdname get_pausing_sites
#' @export
setMethod("get_pausing_sites",
          signature(object = "ribotrans"),
          function(object,
                   sample_selected = NULL,
                   merge_rep = FALSE,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   exclude_length = c(0,0),
                   min_counts = 64,
                   window = 200){
            # ==================================================================
            # whether do reads offset correction
            if(do_offset_correct == TRUE){
              sry <- do_offset_correction(object = object,shift = position_shift)
            }else{
              sry <- object@summary_info
            }

            # filter samples
            if(!is.null(sample_selected)){
              sry <- subset(sry, sample %in% sample_selected)
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
              fastplyr::f_summarise(normsm = sum(norm),total_counts = sum(counts))


            # whether aggregate replicates
            if(merge_rep == TRUE){
              density.tt <- density.tt %>%
                fastplyr::f_group_by(sample_group,rname,relst) %>%
                fastplyr::f_summarise(normsm = mean(normsm),total_counts = mean(total_counts)) %>%
                dplyr::rename(sample = sample_group)
            }else{
              density.tt <- density.tt %>% dplyr::select(-sample_group)
            }

            # get all pos
            rn <- unique(density.tt$rname)
            ft <- object@features %>% dplyr::filter(idnew %in% rn)

            # x= 1
            ft.all <- lapply(seq_along(rn),function(x){
              tmp <- subset(ft , idnew == rn[x])
              df <- data.frame(rname = rn[x], relst = 1:tmp$cds)
            }) %>% do.call("rbind",.) %>% data.frame()

            # add samples
            sps <- unique(density.tt$sample)

            ft.all.sp <- lapply(seq_along(sps),function(x){
              tmp <- ft.all
              tmp$sample <- sps[x]
              return(tmp)
            }) %>% do.call("rbind",.) %>% data.frame()

            # merge density
            ft.all.dt <- ft.all.sp %>%
              dplyr::full_join(y = density.tt,by = c("sample","rname","relst"))
            ft.all.dt[is.na(ft.all.dt)] <- 0

            # ==================================================================
            # calculate pause score
            dt <- data.table::as.data.table(ft.all.dt)

            # get window sums
            dt[, `:=`(
              sum_1 = rolling_window_sum(relst, normsm, window, "sum1"),
              sum_2 = rolling_window_sum(relst, normsm, window, "sum2")
            ), by = .(rname, sample)]

            dt[, pause_score := data.table::fifelse(
              sum_1 * sum_2 == 0, 0,
              (window/2) * normsm * ((sum_1 + sum_2) / (sum_1 * sum_2))
            )]

            dt <- dt %>%
              dplyr::mutate(total_sum = sum_1 + sum_2,
                            min_sum = pmin(sum_1, sum_2),
                            max_sum = pmax(sum_1, sum_2),
                            balance_score = ifelse(total_sum == 0, 0, min_sum / total_sum),
                            imbalance_ratio = ifelse(min_sum == 0, Inf, max_sum / min_sum))

            # add zscore
            dt <- dt %>%
              dplyr::group_by(rname, sample) %>%
              dplyr::mutate(pause_mean = mean(pause_score, na.rm = TRUE),
                            pause_sd = sd(pause_score, na.rm = TRUE),
                            z_score = ifelse(pause_sd == 0, 0, (pause_score - pause_mean) / pause_sd))

            return(dt)

          }
)


# ==============================================================================


# Rcpp::cppFunction('
# NumericVector rolling_window_sum(NumericVector positions, NumericVector expressions,
#                                 double window, String type) {
#   int n = positions.size();
#   NumericVector result(n);
#
#   for(int i = 0; i < n; i++) {
#     double pos_i = positions[i];
#     double sum = 0.0;
#
#     if(type == "sum1") {
#       // window 1: [pos_i, pos_i + window]
#       for(int j = 0; j < n; j++) {
#         if(positions[j] >= pos_i && positions[j] <= pos_i + window) {
#           sum += expressions[j];
#         }
#       }
#     } else {
#       // window 2: [pos_i + window/2, pos_i + 1.5*window]
#       for(int j = 0; j < n; j++) {
#         if(positions[j] >= pos_i + window/2 && positions[j] <= pos_i + 1.5*window) {
#           sum += expressions[j];
#         }
#       }
#     }
#
#     result[i] = sum;
#   }
#
#   return result;
# }
# ')
