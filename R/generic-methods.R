# ==============================================================================
# function to get detailed information for reads
# ==============================================================================

#' @title Generate Summary Information from BAM Files
#'
#' @description
#' This method processes ribosome profiling BAM files and generates positional read count data (footprints)
#' across annotated features. It supports multiple experiment types (such as ribosome footprints, RNA-seq, total RNA, and IP).
#' Mapping positions are adjusted based on user-defined assignment mode (5'-end or 3'-end), and if genome mapping is used,
#' reads are converted to transcriptomic coordinates.
#'
#' @param object A \code{ribotrans} S4 object that contains experiment library information, features, and settings such as mapping type.
#' @param exp_type Character vector specifying experiment types to include. Options are \code{"ribo"}, \code{"rna"}, \code{"total"}, \code{"ip"}.
#' @param load_local Logical; if \code{TRUE}, tries to load precomputed summary data \code{tinfo.anno.rda} from local directory.
#' @param nThreads Integer; number of threads to use when reading BAM files. Passed to \code{Rsamtools::scanBam}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @details
#' If \code{exp_type} includes both \code{"total"} and \code{"ip"}, the \code{sample} and \code{sample_group} columns are suffixed
#' with the experiment type to differentiate datasets.
#'
#' For genome-mapped data, genomic coordinates are translated to transcriptomic positions using a transcript feature annotation table
#' provided in \code{object@genome_trans_features}. If not genome-mapped, data are grouped and summarized directly by transcript ID and position.
#'
#' If \code{load_local = FALSE}, coverage information is read from BAM files and annotated against the transcript features.
#' Data is summarized and merged with transcript metadata (like CDS start and stop).
#' The result is saved as \code{tinfo.anno.rda} in the working directory for future reuse.
#'
#' @return A modified \code{ribotrans} object with the slot \code{summary_info} populated with annotated read counts.
#'
#' @importFrom Rsamtools scanBam ScanBamParam scanBamFlag
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom fastplyr f_group_by f_summarise f_left_join f_select
#' @importFrom dplyr mutate case_when rename n
#'
#'
#' @examples
#' \dontrun{
#' object <- generate_summary(object, exp_type = c("ribo", "rna"), load_local = FALSE)
#' }
#'
#' @export
setGeneric("generate_summary",function(object,...) standardGeneric("generate_summary"))



#' @rdname generate_summary
#' @export
setMethod("generate_summary",
          signature(object = "ribotrans"),
          function(object,
                   exp_type = c("ribo", "rna", "total", "ip"),
                   load_local = FALSE,
                   nThreads = 1){
            options(fastplyr.inform = FALSE)
            # exp_type <- match.arg(exp_type, c("ribo", "rna", "total", "ip"))

            bams <- subset(object@library, type %in% exp_type)
            bfn <- bams$bam
            sp <- bams$sample
            gp <- bams$sample_group

            if(length(exp_type) == 2 & all(exp_type %in% c("total", "ip"))){
              sp <- paste(sp,bams$type, sep = "-")
              gp <- paste(gp,bams$type, sep = "-")
            }


            features <- object@features

            # whether load saved data
            if(load_local == TRUE){
              load(file = "tinfo.anno.rda")
            }else{
              # x = 1
              purrr::map_df(seq_along(bfn),function(x){
                # get position
                bam_data <- Rsamtools::scanBam(file = bfn[x],
                                               nThreads = nThreads,
                                               param = Rsamtools::ScanBamParam(
                                                 what = c("rname", "pos", "strand", "qwidth"),
                                                 flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
                                                                               isSecondaryAlignment = FALSE,
                                                                               isDuplicate = FALSE))
                )

                tinfo <- data.frame(bam_data[[1]])

                # check mapping_type
                if(object@mapping_type == "genome"){
                  if(object@assignment_mode == "end5"){
                    # check strand to assign 5'end position
                    tinfo <- tinfo %>%
                      dplyr::mutate(pos = ifelse(strand == "-", pos + qwidth - 1, pos))
                  }else{
                    # check strand to assign 3'end position
                    tinfo <- tinfo %>%
                      dplyr::mutate(pos = ifelse(strand == "+", pos + qwidth - 1, pos))
                  }

                  tinfo.gr <- tinfo %>%
                    fastplyr::f_group_by(rname,pos,qwidth,strand) %>%
                    fastplyr::f_summarise(count = dplyr::n()) %>%
                    dplyr::rename(seqnames = rname, start = pos) %>%
                    dplyr::mutate(end = start,.after = "start")

                  # to GRanges
                  tinfo.gr <- GenomicRanges::GRanges(tinfo.gr)

                  # transcript annotation regions
                  trans_rg <- GenomicRanges::GRanges(object@genome_trans_features)

                  # genomic position to transcriptome position
                  ov <- IRanges::findOverlaps(query = tinfo.gr,subject = trans_rg, ignore.strand = TRUE)

                  # get overlap data
                  if (requireNamespace("S4Vectors", quietly = TRUE)) {
                    lo <- cbind(as.data.frame(tinfo.gr[S4Vectors::queryHits(ov)]),
                                as.data.frame(trans_rg[S4Vectors::subjectHits(ov)]))
                  } else {
                    warning("Package 'S4Vectors' is needed for this function to work.")
                  }


                  # make unique names
                  names(lo) <- make.names(names(lo),unique = T)

                  # calculate transcript position
                  tinfo <- lo %>%
                    dplyr::mutate(pos = dplyr::case_when(strand == "+" & strand.1 == "+" ~ tx_len - abs(end.1 - start),
                                                         strand == "-" & strand.1 == "+" ~ tx_len - abs(end.1 - (start - qwidth + 1)),
                                                         strand == "-" & strand.1 == "-" ~ tx_len - abs(start - start.1),
                                                         strand == "+" & strand.1 == "-" ~ tx_len - abs(start - (start + qwidth - 1))),
                                  rname = paste(transcript_id,gene_name,sep = "|")) %>%
                    fastplyr::f_select(rname,pos,qwidth,count)

                }else{
                  if(object@assignment_mode == "end3"){
                    # check strand to assign 5'end position
                    tinfo <- tinfo %>%
                      dplyr::mutate(pos = ifelse(strand == "-", pos + qwidth - 1, pos))
                  }

                  tinfo <- tinfo %>%
                    fastplyr::f_group_by(rname,pos,qwidth) %>%
                    fastplyr::f_summarise(count = dplyr::n())

                }

                # add sample and group info
                tinfo$sample <- sp[x]
                tinfo$sample_group <- gp[x]

                return(tinfo)
              }) -> all.tinfo

              # merge with transcript annotation
              tinfo.anno <- all.tinfo %>%
                fastplyr::f_left_join(y = features[,c("idnew","mstart","mstop","translen")],
                                      by = c("rname" = "idnew"))

              # save in local
              save(tinfo.anno, file = "tinfo.anno.rda", compress = TRUE)
            }

            object@summary_info <- tinfo.anno
            return(object)
          }
)





# ==============================================================================
# function to calculate coverage
# ==============================================================================

#' Extract RNA-seq Coverage for a Specific Gene
#'
#' This function extracts transcript- or genome-level RNA-seq coverage (read count and RPM) for a given gene across samples from a \code{ribotrans} object.
#' It supports optional replicates merging, coordinate system conversion, and smoothing using a sliding window.
#'
#' @param object A \code{ribotrans} S4 object containing RNA-seq BAM file information, sample annotations, and transcript models.
#' @param gene_name Character. A gene symbol present in \code{object@features$gene}. Coverage will be extracted for this gene.
#' @param exp Character, either "rna" (default) or "ribo", indicating
#' which experiment type to use.
#' @param smooth Logical. Whether to smooth coverage values using a sliding window (especially useful for visualization). Default is \code{FALSE}.
#' @param coordinate_to_trans Logical. If \code{TRUE} and genome mapping is used, convert genomic coordinates to transcript coordinates. Default is \code{FALSE}.
#' @param merge_rep Logical. If \code{TRUE}, biological replicates in the same \code{sample_group} are merged (mean coverage). Default is \code{FALSE}.
#' @param slide_window Integer. Window size in nucleotides used for smoothing. Only applied when \code{smooth = TRUE}. Default is \code{30}.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This function is used to extract and optionally smooth RNA-seq coverage information at single-nucleotide resolution for a user-specified gene.
#'
#' For each RNA-seq BAM file linked to the \code{ribotrans} object, read coverage is extracted using either transcript or genome coordinates depending on \code{mapping_type}.
#' If replicate merging is specified, then values are averaged per position and per sample group.
#'
#' The final RNA coverage information is stored in the \code{RNA_coverage} slot of the object and optionally smoothed.
#'
#' @return The original \code{ribotrans} object with updated slots:
#' \describe{
#'   \item{\code{RNA_coverage}}{A \code{data.frame} with nucleotide-level RNA coverage (RPM and raw counts).}
#'   \item{\code{RNA.smoothed}}{A character string indicating whether signal smoothing was applied ("TRUE"/"FALSE").}
#'   \item{\code{gene_name}}{The target gene name used in extraction.}
#' }
#'
#' @importFrom dplyr %>% filter mutate select rename group_by summarise left_join
#' @importFrom fastplyr f_group_by f_summarise
#' @importFrom purrr map_df
#'
#' @seealso \code{\link{getCoverage}}, \code{\link{getCoverageGenome}}, \code{\link{smoothEachPosition}}, \code{\link{get_occupancy}}
#'
#' @examples
#' \dontrun{
#'   # Extract raw and smoothed RNA-seq coverage for gene "CDC28"
#'   ribo_obj <- get_coverage(object = ribo_obj,
#'                            gene_name = "CDC28",
#'                            smooth = TRUE,
#'                            merge_rep = TRUE)
#'
#'   # View smoothed signal
#'   head(ribo_obj@RNA_coverage)
#'
#' }
#'
#' @export
setGeneric("get_coverage",function(object,...) standardGeneric("get_coverage"))



#' @rdname get_coverage
#' @export
setMethod("get_coverage",
          signature(object = "ribotrans"),
          function(object, gene_name = NULL,
                   exp = c("rna","ribo"),
                   smooth = FALSE,
                   coordinate_to_trans = FALSE,
                   merge_rep = FALSE,
                   slide_window = 30){
            exp <- match.arg(exp,choices = c("rna","ribo"))

            # check gene name
            if(!(gene_name %in% object@features$gene)){
              stop("Can't find this gene symbol, please have a check!")
            }

            # smooth <- match.arg(smooth, c(FALSE, TRUE))

            if(exp == "ribo"){
              bf <- subset(object@bam_file, type == "ribo")
              bams <- bf$bam
              sp <- bf$sample
              gp <- bf$sample_group

              lib <- subset(object@library, type == "ribo")
              lib <- lib[match(bams, lib$bam),]
              totalreads <- lib$mappped_reads
            }else{
              bf <- subset(object@bam_file, type == "rna")
              bams <- bf$bam
              sp <- bf$sample
              gp <- bf$sample_group

              lib <- subset(object@library, type == "rna")
              lib <- lib[match(bams, lib$bam),]
              totalreads <- lib$mappped_reads
            }


            # loop for each sample
            purrr::map_df(seq_along(bams),function(x){
              # check mapping type
              if(object@mapping_type == "genome"){
                tmp <- getCoverageGenome(bam_file = bams[x],
                                         gene_name = gene_name,
                                         features = object@genome_trans_features,
                                         total_reads = totalreads[x],
                                         coordinate_to_trans = coordinate_to_trans)


                # check tmp
                if(nrow(tmp) == 0){
                  if(coordinate_to_trans == TRUE){
                    gn <- gene_name
                    gf <- subset(object@features, gene == gn)

                    tmp <- data.frame(rname = gf$idnew[1], pos = 1,
                                      count = 0, rpm = 0)
                  }else{
                    gn <- gene_name
                    gf <- subset(object@genome_trans_features, gene_name == gn)

                    tmp <- data.frame(rname = gf$seqnames[1], pos = gf$start,
                                      count = 0, rpm = 0)
                  }

                }

              }else{
                tmp <- getCoverage(bam_file = bams[x],
                                   gene_name = gene_name,
                                   features = object@features,
                                   total_reads = totalreads[x])

                # check tmp
                if(nrow(tmp) == 0){
                  gn <- gene_name
                  gf <- subset(object@features, gene == gn)

                  tmp <- data.frame(rname = gf$idnew[1], pos = 1,
                                    count = 0, rpm = 0)

                }
              }



              # add sample name
              tmp$sample <- sp[x]
              tmp$sample_group <- gp[x]

              return(tmp)
            }) -> rna.df

            # whether aggregate replicates
            if(merge_rep == TRUE){
              rna.df <- rna.df %>%
                fastplyr::f_group_by(sample_group,rname,pos) %>%
                fastplyr::f_summarise(count = mean(count),rpm = mean(rpm)) %>%
                dplyr::rename(sample = sample_group)
            }


            # smooth for each position
            cond <- object@mapping_type == "genome" & coordinate_to_trans == TRUE | object@mapping_type == "transcriptome"
            if(cond){
              if(smooth == TRUE){
                sm.df <- smoothEachPosition(features = object@features,
                                            posdf = rna.df,
                                            slide_window = slide_window)
              }else{
                sm.df <- rna.df
              }
            }else{
              sm.df <- rna.df
            }

            if(exp == "ribo"){
              if(smooth == TRUE){
                object@ribo_occupancy <- sm.df
                object@ribo.smoothed <- "TRUE"
              }else{
                sm.df$smooth <- sm.df$rpm
                object@ribo_occupancy <- sm.df
                object@ribo.smoothed <- "FALSE"
              }
            }else{
              if(smooth == TRUE){
                object@RNA_coverage <- sm.df
                object@RNA.smoothed <- "TRUE"
              }else{
                sm.df$smooth <- sm.df$rpm
                object@RNA_coverage <- sm.df
                object@RNA.smoothed <- "FALSE"
              }
            }

            object@gene_name <- gene_name

            return(object)
          }
)




# ==============================================================================
# function to calculate ribosome occupancy
# ==============================================================================

#' Extract Ribosome Occupancy for a Specific Gene
#'
#' Retrieves position-wise ribosome footprint occupancy and RPM values from a ribotrans object for a specified gene. Supports optional offset correction, replicate merging, coordinate transformations, and sliding window smoothing.
#'
#' @param object A \code{ribotrans} S4 object containing processed ribo-seq alignment data and metadata.
#' @param gene_name Character. The gene symbol (as in \code{features$gene}) for which occupancy should be extracted.
#' @param smooth Logical. Whether to smooth ribosome density (RPM) using a sliding window across positions. Default: \code{FALSE}.
#' @param coordinate_to_trans Logical. If \code{TRUE} and mapping is from genome, will convert coordinates to the transcript coordinate space. Default: \code{FALSE}.
#' @param do_reads_offset Logical. Whether to apply P-site offset correction based on pre-learned offsets stored in \code{object@reads_offset_info}. Default: \code{FALSE}.
#' @param merge_rep Logical. If \code{TRUE}, average replicate samples within the same sample_group. Default: \code{FALSE}.
#' @param slide_window Integer. The window size to use in smoothing. Only relevant if \code{smooth = TRUE}. Default: \code{30}.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This function is intended to extract transcript- or genome-level ribosome occupancy profiles for a single gene across multiple samples.
#' It supports:
#' \itemize{
#'   \item Extracting counts/RPM mapped to either genome or transcript coordinates.
#'   \item P-site offset correction per sample/read length.
#'   \item Merging replicates based on sample group.
#'   \item Optional smoothing of RPM signal using a sliding window.
#' }
#'
#' If offset correction is enabled (\code{do_reads_offset = TRUE}), the method uses the offset table in \code{object@reads_offset_info} to adjust footprint positions.
#'
#' @return The original \code{ribotrans} object with added/updated slots:
#' \itemize{
#'   \item \code{ribo_occupancy}: A \code{data.frame} containing gene-specific ribosome profiles (position, count, rpm).
#'   \item \code{ribo.smoothed}: A character indicating whether smoothing was applied ("TRUE"/"FALSE").
#'   \item \code{gene_name}: The gene symbol analyzed.
#' }
#'
#' @importFrom dplyr %>% filter mutate select rename group_by summarise left_join
#' @importFrom fastplyr f_group_by f_summarise
#' @importFrom purrr map_df
#' @importFrom stats na.omit
#'
#' @seealso \code{\link{getOccupancy}}, \code{\link{getOccupancyGenome}}, \code{\link{smoothEachPosition}}
#'
#' @examples
#' \dontrun{
#'   # Extract ribosome density for gene "CDC4" across samples
#'   ribo_obj <- get_occupancy(object = ribo_obj,
#'                             gene_name = "CDC4",
#'                             smooth = TRUE,
#'                             do_reads_offset = TRUE,
#'                             merge_rep = FALSE)
#'
#'   head(ribo_obj@ribo_occupancy)
#'
#' }
#'
#' @export
setGeneric("get_occupancy",function(object,...) standardGeneric("get_occupancy"))



#' @rdname get_occupancy
#' @export
setMethod("get_occupancy",
          signature(object = "ribotrans"),
          function(object, gene_name = NULL,
                   smooth = FALSE,
                   coordinate_to_trans = FALSE,
                   do_reads_offset = FALSE,
                   merge_rep = FALSE,
                   slide_window = 30){
            # check gene name
            if(!(gene_name %in% object@features$gene)){
              stop("Can't find this gene symbol, please have a check!")
            }

            # smooth <- match.arg(smooth, c(FALSE, TRUE))
            assignment_mode <- object@assignment_mode

            bf <- subset(object@bam_file, type == "ribo")
            ribobams <- bf$bam
            sp <- bf$sample
            gp <- bf$sample_group

            lib <- subset(object@library, type == "ribo")
            lib <- lib[match(ribobams, lib$bam),]
            totalreads <- lib$mappped_reads

            # loop for each sample
            purrr::map_df(seq_along(ribobams),function(x){
              # check mapping type
              if(object@mapping_type == "genome"){
                tmp <- getOccupancyGenome(bam_file = ribobams[x],
                                          gene_name = gene_name,
                                          features = object@genome_trans_features,
                                          total_reads = totalreads[x],
                                          assignment_mode = assignment_mode,
                                          coordinate_to_trans = coordinate_to_trans)

                # check tmp
                if(nrow(tmp) == 0){
                  if(coordinate_to_trans == TRUE){
                    gn <- gene_name
                    gf <- subset(object@features, gene == gn)

                    tmp <- data.frame(rname = gf$idnew[1], pos = 1,
                                      strand = "+", qwidth = 29,
                                      count = 0, rpm = 0)
                  }else{
                    gn <- gene_name
                    gf <- subset(object@genome_trans_features, gene_name == gn)

                    tmp <- data.frame(rname = gf$seqnames[1], pos = gf$start,
                                      strand = gf$strand, qwidth = 29,
                                      count = 0, rpm = 0)
                  }
                }
              }else{
                tmp <- getOccupancy(bam_file = ribobams[x],
                                    gene_name = gene_name,
                                    features = object@features,
                                    total_reads = totalreads[x])

                # check tmp
                if(nrow(tmp) == 0){
                  gn <- gene_name
                  gf <- subset(object@features, gene == gn)

                  tmp <- data.frame(rname = gf$idnew[1], pos = 1,
                                    strand = "+", qwidth = 29,
                                    count = 0, rpm = 0)

                }
              }

              # add sample name
              tmp$sample <- sp[x]
              tmp$sample_group <- gp[x]

              # ================================================================
              # whether do read offset
              if(do_reads_offset == TRUE){
                # check offset data
                if(nrow(object@reads_offset_info) == 0){
                  stop("No reads offset information in object!")
                }

                offset <- object@reads_offset_info

                if(object@mapping_type == "genome" & coordinate_to_trans == FALSE){
                  adjusted <- tmp %>%
                    dplyr::left_join(y = offset[,c("sample","rel_pos","qwidth")],by = c("sample","qwidth")) %>%
                    stats::na.omit() %>%
                    dplyr::mutate(pos = dplyr::if_else(strand == "+",
                                                       pos - rel_pos,pos + rel_pos)) %>%
                    dplyr::select(-strand,-qwidth,-rel_pos)
                }else{
                  adjusted <- tmp %>%
                    dplyr::left_join(y = offset[,c("sample","rel_pos","qwidth")],by = c("sample","qwidth")) %>%
                    stats::na.omit() %>%
                    dplyr::mutate(pos = pos - rel_pos) %>%
                    dplyr::select(-qwidth,-rel_pos)
                }

              }else{
                adjusted <- tmp
              }

              adjusted <- adjusted %>%
                dplyr::group_by(sample_group,sample,rname,pos) %>%
                dplyr::summarise(count = sum(count),rpm = sum(rpm), .groups = "drop")

              return(adjusted)
            }) -> ribo.df

            # whether aggregate replicates
            if(merge_rep == TRUE){
              ribo.df <- ribo.df %>%
                fastplyr::f_group_by(sample_group,rname,pos) %>%
                fastplyr::f_summarise(count = mean(count),rpm = mean(rpm)) %>%
                dplyr::rename(sample = sample_group)
            }

            # smooth for each position
            cond <- object@mapping_type == "genome" & coordinate_to_trans == TRUE | object@mapping_type == "transcriptome"
            if(cond){
              if(smooth == TRUE){
                sm.df <- smoothEachPosition(features = object@features,
                                            posdf = ribo.df,
                                            slide_window = slide_window)
              }else{
                sm.df <- ribo.df
              }
            }else{
              sm.df <- ribo.df
            }

            if(smooth == TRUE){
              object@ribo_occupancy <- sm.df
              object@ribo.smoothed <- "TRUE"
            }else{
              sm.df$smooth <- sm.df$rpm
              object@ribo_occupancy <- sm.df
              object@ribo.smoothed <- "FALSE"
            }

            object@gene_name <- gene_name

            return(object)
          }
)



# ==============================================================================
# scale to mRNA abundance
# ==============================================================================

#' @title Compute Scaled Ribosome Occupancy
#' @description This generic function calculates scaled ribosome occupancy by normalizing ribosome occupancy
#'     with RNA-seq coverage.
#' @param object A `ribotrans` object containing RNA-seq coverage and ribosome occupancy.
#' @param ... Additional arguments passed to specific methods.
#' @return A modified `ribotrans` object with scaled occupancy stored in `object@scaled_occupancy`.
#' @seealso \code{\link{get_scaled_occupancy,ribotrans-method}} for the method implementation.
#' @export
#' @rdname get_scaled_occupancy
setGeneric("get_scaled_occupancy",function(object,...) standardGeneric("get_scaled_occupancy"))


#' @title Compute Scaled Ribosome Occupancy for a `ribotrans` Object
#' @description This method calculates the scaled ribosome occupancy by dividing ribosome density
#'     (`ribo_occupancy`) by RNA-seq coverage (`RNA_coverage`), optionally applying a smoothing function.
#'
#' @param object A `ribotrans` object containing both `ribo_occupancy` and `RNA_coverage` data.
#' @param smooth Character. Should the data be smoothed? Options: `TRUE` or `FALSE`. Default is `TRUE`.
#' @param slide_window Integer. The window size for smoothing (used only if `smooth = "TRUE"`). Default: `30`.
#'
#' @details
#' **Functionality:**
#' This method normalizes ribosome occupancy by RNA-seq coverage per position, producing a **scaled occupancy profile**.
#'
#' **Workflow:**
#' 1. Ensure both `RNA_coverage` (`object@RNA_coverage`) and `ribo_occupancy` (`object@ribo_occupancy`) exist.
#' 2. Compute: `scaled occupancy = (ribo occupancy) / (RNA coverage)`.
#' 3. Replace `NA` or `Inf` values with `0` for stability.
#' 4. **Optional:** Apply a rolling mean (`zoo::rollmean()`) for smoothing.
#'
#' **Error Handling:**
#' - If `RNA_coverage` or `ribo_occupancy` is missing, the function stops with an error.
#'
#' @return A modified `ribotrans` object containing `object@scaled_occupancy`, a `data.frame` with the scaled enrichment values.
#' @seealso
#' \code{\link{get_occupancy}}: Extract raw ribosome occupancy.
#' \code{\link{get_coverage}}: Extract RNA coverage.
#' \code{\link{ribotrans}}: The S4 class storing input data.
#'
#' @examples
#' \dontrun{
#' # Assume `ribo_obj` is a ribotrans object with RNA and ribosome data
#' result <- get_scaled_occupancy(ribo_obj, smooth = "TRUE", slide_window = 50)
#' print(result@scaled_occupancy)
#' }
#'
#' @importFrom dplyr select left_join mutate filter
#' @importFrom purrr map_df
#' @export
#' @rdname get_scaled_occupancy
setMethod("get_scaled_occupancy",
          signature(object = "ribotrans"),
          function(object,
                   smooth = TRUE,
                   slide_window = 30){
            # smooth <- match.arg(smooth, c(TRUE, FALSE))

            # ribo rpm/rna rpm for each position

            # check rna coverage
            if(is.null(object@RNA_coverage)){
              stop("Please supply RNA coverage data!")
            }

            if(is.null(object@ribo_occupancy)){
              stop("Please supply Ribo occupancy data!")
            }

            # run
            mb <- object@ribo_occupancy %>%
              dplyr::select(sample,rname,pos,rpm) %>%
              dplyr::left_join(y = object@RNA_coverage[,c("sample","rname","pos","rpm")],
                               by = c("sample","rname","pos")) %>%
              dplyr::mutate(enrich = (rpm.x/rpm.y)) %>%
              dplyr::mutate(enrich = ifelse(is.na(enrich) | is.infinite(enrich),0,enrich))

            # loop for each sample and transcript to smooth
            tanno <- subset(object@features, gene == object@gene_name)

            # check whether smooth data
            purrr::map_df(1:nrow(tanno),function(x){
              tmp <- tanno[x,]

              tmp2 <- mb %>%
                dplyr::filter(rname == tmp$idnew)

              sp <- unique(tmp2$sample)

              # loop for each sample
              purrr::map_df(seq_along(sp),function(x){
                tmp3 <- tmp2 %>%
                  dplyr::filter(sample == sp[x])

                pos.df <- data.frame(rname = tmp$idnew, pos = 1:tmp$translen)

                # merge with value
                pos.df <- pos.df %>%
                  dplyr::left_join(y = tmp3,by = c("rname", "pos"))
                pos.df$sample <- tmp3$sample[1]
                pos.df[is.na(pos.df)] <- 0

                if(smooth == TRUE){
                  # smooth data
                  if (requireNamespace("zoo", quietly = TRUE)) {
                    pos.df$smooth <- zoo::rollmean(pos.df$enrich, k = slide_window, fill = 0)
                  } else {
                    warning("Package 'zoo' is needed for this function to work.")
                  }
                }else{
                  pos.df$smooth <- pos.df$enrich
                }

                return(pos.df)
              }) -> smooth.df

              return(smooth.df)
            }) -> mb.smoothed

            object@scaled_occupancy <- mb.smoothed
            return(object)
          }
)




