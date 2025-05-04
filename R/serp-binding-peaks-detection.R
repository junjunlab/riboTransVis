# ==============================================================================
# find strong peaks
# ==============================================================================

#' Find Strong Binding Peaks from SeRP Data
#'
#' Identifies strong signal peaks in selective ribosome profiling (SeRP) data by comparing the IP (immunoprecipitated) and total samples. This function can optionally apply P-site offset correction, smooth read signal, calculate enrichment ratio, filter low-count features, and detect enriched regions (peaks) within coding sequences (CDSs).
#'
#' @param object A \code{serp} object containing summarized read alignment data, sample metadata, and transcript annotations.
#' @param do_offset_correct Logical. Whether to apply P-site offset correction. Default is \code{FALSE}.
#' @param position_shift Integer. Fixed offset to apply when \code{do_offset_correct = TRUE}. Default is \code{0}.
#' @param min_counts Integer. Minimum total count required for a gene in a sample to be included in the analysis. Default is \code{32}.
#' @param read_length Integer vector of length 2. Specifies the range of read lengths to include in the analysis. Default is \code{c(25, 35)}.
#' @param smooth Logical. Whether to apply signal smoothing using a rolling sum of counts. Default is \code{TRUE}.
#' @param window_size Integer. Size of the rolling window used for smoothing. Only used if \code{smooth = TRUE}. Default is \code{45}.
#' @param enrichment_threshold Numeric. Minimum enrichment ratio (IP / total) required to consider a region as a strong peak. Default is \code{1.5}.
#' @param binding_width Integer. Minimum width (in nucleotides or codons) of a contiguous enriched region to be considered a peak. Default is \code{15}.
#' @param mode Character. Output resolution of the coordinates. Either \code{"nt"} (nucleotide-level) or \code{"codon"} (codon-level). Default is \code{"nt"}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @return A \code{data.frame} containing the detected strong binding peaks with the following columns:
#' \item{sample}{Sample identifier}
#' \item{sample_group}{Group corresponding to the sample}
#' \item{rname}{Transcript or gene name}
#' \item{start}{Start position of the peak (in nt or codon)}
#' \item{end}{End position of the peak (in nt or codon)}
#' \item{mean_ratio}{Average enrichment ratio in the peak window}
#' \item{binding_width}{Width of the binding peak (number of positions)}
#'
#' @details
#' The function compares RPM-normalized read counts from total and IP libraries. It identifies regions where the enrichment ratio (IP/total) exceeds a given threshold for at least a minimum number of nucleotides (or codons). Optional signal smoothing via RcppRoll can help reduce local noise.
#'
#' When \code{mode = "codon"}, the data is downsampled to codon resolution after peak calling, with positions divided by 3 (frame-shifted if needed).
#'
#' Read filtering is applied based on read length and minimum count thresholds.
#'
#' Offset correction (if enabled) aligns reads based on inferred P-site positions.
#'
#' @importFrom dplyr %>% filter select mutate group_by summarise left_join ungroup
#' @importFrom data.table as.data.table rleid tstrsplit
#' @importFrom fastplyr f_inner_join f_left_join f_full_join
#' @importFrom purrr map
#'
#'
#' @seealso \code{\link{do_offset_correction}}, \code{\link{serp-class}}
#'
#' @examples
#' \dontrun{
#' peaks <- find_strong_peaks(object = my_serp,
#' do_offset_correct = TRUE,
#' position_shift = 12,
#' min_counts = 32,
#' read_length = c(25, 35),
#' smooth = TRUE,
#' window_size = 45,
#' enrichment_threshold = 1.5,
#' binding_width = 15,
#' mode = "codon")
#' head(peaks)
#' }
#'
#' @rdname find_strong_peaks
#' @export
setGeneric("find_strong_peaks",function(object,...) standardGeneric("find_strong_peaks"))





#' @rdname find_strong_peaks
#' @export
setMethod("find_strong_peaks",
          signature(object = "serp"),
          function(object,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   min_counts = 32,
                   read_length = c(25,35),
                   smooth = TRUE,
                   window_size = 45,
                   enrichment_threshold = 1.5,
                   binding_width = 15,
                   mode = c("nt", "codon")){
            mode <- match.arg(mode, choices = c("nt", "codon"))

            # ==================================================================
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


            # retain CDS region
            sry <- sry[pos >= mstart & pos <= mstop][, pos := pos - mstart + 1][, .(count = sum(count)), by = .(sample_group, sample, rname, pos)]

            # filter low counts for each sample
            ct <- sry[, .(counts = sum(count)), by = .(sample, rname)][counts >= min_counts][,`:=`(counts = NULL)]
            sry <- sry %>% fastplyr::f_inner_join(y = ct,by = c("sample", "rname"))

            # rpm normalization
            lib <- object@library %>%
              dplyr::mutate(sample = paste(sample,type,sep = "-")) %>%
              dplyr::select(sample,mappped_reads)

            sry.tmp <- sry[lib, on = .(sample), nomatch = 0][, rpm := (count / mappped_reads) * 10^6][, `:=`(mappped_reads = NULL, count = NULL)]

            # ===================================================================
            # run finding peaks
            lib <- object@library
            tt_sp <- subset(lib, type == "total")[,c("type", "sample", "sample_group")]

            ip_sp <- subset(lib, type == "ip")[,c("type", "sample", "sample_group")]

            sp_mr <- tt_sp %>% dplyr::left_join(y = ip_sp,by = "sample")

            sp.tt <- paste(sp_mr$sample,sp_mr$type.x,sep = "-")
            sp.ip <- paste(sp_mr$sample,sp_mr$type.y,sep = "-")

            tmp.tt <- subset(sry.tmp, sample == sp.tt)
            tmp.tt[, sample := data.table::tstrsplit(sample, "-total", fixed = TRUE)[[1]]]
            tmp.tt[, sample_group := data.table::tstrsplit(sample_group, "-total", fixed = TRUE)[[1]]]

            tmp.ip <- subset(sry.tmp, sample == sp.ip)
            tmp.ip[, sample := data.table::tstrsplit(sample, "-ip", fixed = TRUE)[[1]]]
            tmp.ip[, sample_group := data.table::tstrsplit(sample_group, "-ip", fixed = TRUE)[[1]]]

            # merge
            serp.mer <- tmp.tt %>%
              fastplyr::f_full_join(y = tmp.ip,by = c("sample","sample_group","rname","pos"))

            ft <- object@features %>%
              dplyr::ungroup() %>%
              dplyr::select(idnew,cds) %>%
              dplyr::rename(rname = idnew)

            sp <- unique(serp.mer$sample)

            # loop for each sample to find peaks
            # x = 1
            lapply(seq_along(sp),function(x){
              tmp <- subset(serp.mer,sample == sp[x])

              ids <- unique(tmp$rname)

              gpos <- ft %>%
                dplyr::filter(rname %in% ids) %>%
                dplyr::mutate(pos = purrr::map(cds, ~1:.x)) %>%
                tidyr::unnest(pos) %>%
                dplyr::select(-cds)

              serp.mer.pos <- gpos %>%
                fastplyr::f_left_join(y = tmp,by = c("rname","pos"))

              serp.mer.pos$sample_group <- tmp$sample_group[1]
              serp.mer.pos$sample <- tmp$sample[1]

              serp.mer.pos[is.na(serp.mer.pos)] <- 0

              # soomth data for each position
              if(smooth == TRUE){
                if (requireNamespace("RcppRoll", quietly = TRUE)) {
                  serp.mer.pos <- serp.mer.pos %>%
                    dplyr::group_by(sample,sample_group,rname) %>%
                    dplyr::mutate(tt = RcppRoll::roll_sum(rpm.x, n = window_size, fill = 0, align = "center"),
                                  ip = RcppRoll::roll_sum(rpm.y, n = window_size, fill = 0, align = "center"),
                                  ratio = ip / tt) %>%
                    dplyr::mutate(ratio = ifelse(is.na(ratio) | is.infinite(ratio),0,ratio))
                } else {
                  warning("Package 'RcppRoll' is needed for this function to work.")
                }

              }else{
                serp.mer.pos <- serp.mer.pos %>%
                  dplyr::mutate(tt = rpm.x,ip = rpm.y,ratio = ip / tt) %>%
                  dplyr::mutate(ratio = ifelse(is.na(ratio) | is.infinite(ratio),0,ratio))
              }


              # ============================================================
              # find strong bind peaks
              strong.peaks <- serp.mer.pos %>%
                dplyr::mutate(above = ratio >= enrichment_threshold) %>%
                # add group id
                dplyr::mutate(group = data.table::rleid(above)) %>%
                dplyr::group_by(sample,sample_group,rname,group, above) %>%
                dplyr::summarise(start = dplyr::first(pos),
                                 end = dplyr::last(pos),
                                 mean_ratio = mean(ratio),
                                 .groups = "drop") %>%
                # keep density > threshold and min_len > threshould
                dplyr::filter(above, end - start + 1 >= binding_width) %>%
                dplyr::mutate(binding_width = end - start + 1) %>%
                dplyr::select(-group)


              return(strong.peaks)
            }) %>% do.call("rbind",.) %>% data.frame() -> peaks.df


            # check show mode
            if(mode == "codon"){
              peaks.df <- peaks.df %>%
                dplyr::mutate(start = (start - 1) %/% 3 + 1,
                              end = (end - 1) %/% 3 + 1,
                              binding_width = end - start + 1)
            }

            return(peaks.df)
          }
)





# ==============================================================================
# find weak peaks
# ==============================================================================

#' Find Weak Binding Peaks in Selective Ribo-seq (SeRP) Data
#'
#' This function identifies weak binding peaks in IP (selective) ribosome profiling (SeRP) data. Peaks are regions in a transcript where the normalized read coverage (RPM) exceeds a specified fold-change relative to the transcript's background region (typically the first few codons/nucleotides) for at least a minimum width. This analysis helps discover candidate RNA-binding protein (RBP) footprints that are weaker but reproducible.
#'
#' @param object A \code{serp} object containing processed SeRP data, including normalized read coverage (RPM), transcript annotations, and sample metadata.
#' @param do_offset_correct Logical. Whether to apply P-site offset correction. Default is \code{FALSE}.
#' @param position_shift Integer. Position shift applied during offset correction. Only used when \code{do_offset_correct = TRUE}. Default is \code{0}.
#' @param min_counts Integer. Minimum total read counts required per gene per sample to be included in the analysis. Default is \code{32}.
#' @param read_length Integer vector of length 2 representing the minimum and maximum read lengths to retain for analysis. Default is \code{c(25, 35)}.
#' @param smooth Logical. Whether to smooth read signal using a rolling window. Default is \code{TRUE}.
#' @param window_size Integer. Size of the rolling window width for smoothing. Only used if \code{smooth = TRUE}. Default is \code{45}.
#' @param background_width Integer. Width of the initial region in each transcript used to estimate the background RPM signal. Default is \code{90}.
#' @param enrichment_threshold Numeric. Minimum enrichment ratio (rpm/background_mean) required for a genome region to be classified as a peak. Typical value is \code{3}. Default is \code{3}.
#' @param binding_width Integer. Minimum width (in nucleotides or codons) of a contiguous high-ratio region for it to qualify as a peak. Default is \code{15}.
#' @param mode Character. Indicates whether the output coordinates refer to "nt" (nucleotides) or "codon" (triplet resolution). Default is \code{"nt"}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @return A \code{data.frame} with identified weak binding peaks, containing the following columns:
#' \item{sample}{Sample identifier}
#' \item{sample_group}{Group to which the sample belongs}
#' \item{rname}{Transcript or gene name}
#' \item{start}{Start position of the peak (in nucleotides or codons)}
#' \item{end}{End position of the peak (in nucleotides or codons)}
#' \item{mean_ratio}{Average enrichment ratio within the peak}
#' \item{binding_width}{Peak width (in nucleotides or codons)}
#'
#' @details
#' This function only analyzes "-ip" type SeRP samples to find weak binding events. The RPM signal is compared to a background signal estimated from the first \code{background_width} nucleotides of the coding region. Smoothing can enhance reproducibility and reduce noise.
#'
#' Peaks are filtered by a minimum enrichment threshold and length. Optionally, codon-level resolution can be returned.
#'
#' @section Signal Smoothing:
#' If \code{smooth = TRUE}, the function uses a centered rolling sum with window size \code{window_size} via the \pkg{RcppRoll} package.
#'
#' @seealso \code{\link{find_strong_peaks}}, \code{\link{do_offset_correction}}, \code{\link{serp-class}}
#'
#' @importFrom dplyr %>% filter mutate select group_by ungroup summarise left_join
#' @importFrom data.table as.data.table rleid tstrsplit
#' @importFrom purrr map
#' @importFrom fastplyr f_inner_join f_left_join f_group_by
#'
#'
#' @examples
#' \dontrun{
#' peaks <- find_weak_peaks(object = my_serp,
#' do_offset_correct = TRUE,
#' position_shift = 12,
#' min_counts = 20,
#' read_length = c(25, 34),
#' smooth = TRUE,
#' window_size = 45,
#' background_width = 90,
#' enrichment_threshold = 3,
#' binding_width = 15,
#' mode = "codon")
#' head(peaks)
#' }
#'
#' @rdname find_weak_peaks
#' @export
setGeneric("find_weak_peaks",function(object,...) standardGeneric("find_weak_peaks"))






#' @rdname find_weak_peaks
#' @export
setMethod("find_weak_peaks",
          signature(object = "serp"),
          function(object,
                   do_offset_correct = FALSE,
                   position_shift = 0,
                   min_counts = 32,
                   read_length = c(25,35),
                   smooth = TRUE,
                   window_size = 45,
                   background_width = 90,
                   enrichment_threshold = 3,
                   binding_width = 15,
                   mode = c("nt", "codon")){
            mode <- match.arg(mode, choices = c("nt", "codon"))

            # ==================================================================
            features <- object@features

            # whether do reads offset correction
            if(do_offset_correct == TRUE){
              sry <- do_offset_correction(object = object,shift = position_shift)
            }else{
              sry <- object@summary_info
            }

            sry <- data.table::as.data.table(sry)

            # filter data
            sry <- sry[(mstart != 0 | mstop != 0)
                       & between(qwidth, read_length[1], read_length[2])
                       & endsWith(sample, "-ip")]

            # retain CDS region
            sry <- sry[pos >= mstart & pos <= mstop][, pos := pos - mstart + 1][, .(count = sum(count)), by = .(sample_group, sample, rname, pos)]

            # filter low counts for each sample
            ct <- sry[, .(counts = sum(count)), by = .(sample, rname)][counts >= min_counts][,`:=`(counts = NULL)]
            sry <- sry %>% fastplyr::f_inner_join(y = ct,by = c("sample", "rname"))

            # rpm normalization
            lib <- object@library %>%
              dplyr::mutate(sample = paste(sample,type,sep = "-")) %>%
              dplyr::select(sample,mappped_reads)

            sry.tmp <- sry[lib, on = .(sample), nomatch = 0][, rpm := (count / mappped_reads) * 10^6][, `:=`(mappped_reads = NULL, count = NULL)]

            # ===================================================================
            # run finding peaks

            tmp.ip <- sry.tmp
            tmp.ip[, sample := data.table::tstrsplit(sample, "-ip", fixed = TRUE)[[1]]]
            tmp.ip[, sample_group := data.table::tstrsplit(sample_group, "-ip", fixed = TRUE)[[1]]]


            ft <- object@features %>%
              dplyr::ungroup() %>%
              dplyr::select(idnew,cds) %>%
              dplyr::rename(rname = idnew)

            sp <- unique(tmp.ip$sample)

            # loop for each sample to find peaks
            # x = 1
            lapply(seq_along(sp),function(x){
              tmp <- subset(tmp.ip,sample == sp[x])

              ids <- unique(tmp$rname)

              gpos <- ft %>%
                dplyr::filter(rname %in% ids) %>%
                dplyr::mutate(pos = purrr::map(cds, ~1:.x)) %>%
                tidyr::unnest(pos) %>%
                dplyr::select(-cds)

              serp.mer.pos <- gpos %>%
                fastplyr::f_left_join(y = tmp,by = c("rname","pos"))

              serp.mer.pos$sample_group <- tmp$sample_group[1]
              serp.mer.pos$sample <- tmp$sample[1]

              serp.mer.pos[is.na(serp.mer.pos)] <- 0

              # soomth data for each position
              if(smooth == TRUE){
                if (requireNamespace("RcppRoll", quietly = TRUE)) {
                  serp.mer.pos <- serp.mer.pos %>%
                    dplyr::group_by(sample,sample_group,rname) %>%
                    dplyr::mutate(ip = RcppRoll::roll_sum(rpm, n = window_size, fill = 0, align = "center"))
                } else {
                  warning("Package 'RcppRoll' is needed for this function to work.")
                }

              }else{
                serp.mer.pos <- serp.mer.pos %>%
                  dplyr::mutate(ip = rpm)
              }


              # ============================================================
              # find weak bind peaks
              weak.peaks.tmp <- serp.mer.pos %>%
                fastplyr::f_group_by(sample,sample_group,rname) %>%
                dplyr::mutate(background_mean = ifelse(max(pos) <= background_width,
                                                       mean(ip, na.rm = TRUE),
                                                       mean(ip[pos <= background_width & pos >= 1], na.rm = TRUE))) %>%
                dplyr::mutate(ratio = ip / background_mean) %>%
                dplyr::mutate(ratio = ifelse(is.na(ratio) | is.infinite(ratio),0,ratio)) %>%
                dplyr::ungroup()

              # find peaks
              weak.peaks <- weak.peaks.tmp %>%
                dplyr::mutate(above = ratio >= enrichment_threshold) %>%
                # add group id
                dplyr::mutate(group = data.table::rleid(above)) %>%
                dplyr::group_by(sample,sample_group,rname,group, above) %>%
                dplyr::summarise(start = dplyr::first(pos),end = dplyr::last(pos),
                                 mean_ratio = mean(ratio),
                                 .groups = "drop") %>%
                # keep density > threshold and min_len > threshould
                dplyr::filter(above, end - start + 1 >= binding_width) %>%
                dplyr::mutate(binding_width = end - start + 1) %>%
                dplyr::select(-group)

              return(weak.peaks)
            }) %>% do.call("rbind",.) %>% data.frame() -> peaks.df


            # check show mode
            if(mode == "codon"){
              peaks.df <- peaks.df %>%
                dplyr::mutate(start = (start - 1) %/% 3 + 1,
                              end = (end - 1) %/% 3 + 1,
                              binding_width = end - start + 1)
            }

            return(peaks.df)
          }
)
