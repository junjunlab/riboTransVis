# ==============================================================================
# function to get detailed information for reads
# ==============================================================================

#' Generate Summary Information for Ribosome Profiling Data
#'
#' @description Generic function to generate summary information from ribosome profiling
#' or RNA-seq data.
#'
#' @param object An object containing experimental data
#' @param ... Additional arguments passed to methods
#'
#' @return An updated object containing summary information
#'
#' @export
setGeneric("generate_summary",function(object,...) standardGeneric("generate_summary"))



#' generate_summary Generate summary for ribotrans objects
#'
#' @description This method extracts read mapping information from BAM files and links it
#' to transcript features in a `ribotrans` object.
#'
#' @param object A `ribotrans` object containing library and feature information.
#' @param exp_type Character string indicating experiment type. One of `"ribo"` (ribosome profiling)
#' or `"rna"` (RNA-seq).
#' @param nThreads Integer specifying the number of threads to use for BAM file processing.
#'
#' @return The input `ribotrans` object with an updated `summary_info` slot containing read count summaries.
#'
#' @details This method processes BAM files corresponding to the specified experiment type (`ribo` or `rna`) and extracts:
#' \itemize{
#'   \item Read positions (`pos`)
#'   \item Read lengths (`qwidth`)
#'   \item Read counts per transcript position
#' }
#' The extracted mapping data is then merged with transcript annotation data from the `features` slot.
#'
#' @importFrom Rsamtools scanBam ScanBamParam scanBamFlag
#' @importFrom purrr map_df
#' @importFrom fastplyr f_group_by f_summarise f_left_join
#'
#' @export
setMethod("generate_summary",
          signature(object = "ribotrans"),
          function(object, exp_type = c("ribo", "rna"),
                   nThreads = 1){
            exp_type <- match.arg(exp_type, c("ribo", "rna"))

            bams <- subset(object@library, type == exp_type)
            bfn <- bams$bam
            gp <- bams$sample

            features <- object@features

            # x = 1
            purrr::map_df(seq_along(bfn),function(x){
              # get position
              bam_data <- Rsamtools::scanBam(file = bfn[x],
                                             nThreads = nThreads,
                                             param = Rsamtools::ScanBamParam(what = c("rname", "pos", "qwidth"),
                                                                             flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
              )

              tinfo <- data.frame(bam_data[[1]]) %>%
                fastplyr::f_group_by(rname,pos,qwidth) %>%
                fastplyr::f_summarise(count = n())

              tinfo$sample <- gp[x]

              return(tinfo)
            }) -> all.tinfo

            # merge with transcript annotation
            tinfo.anno <- all.tinfo %>%
              fastplyr::f_left_join(y = features[,c("idnew","mstart","mstop","translen")],
                                    by = c("rname" = "idnew"))


            object@summary_info <- tinfo.anno
            return(object)
          }
)





# ==============================================================================
# function to calculate coverage
# ==============================================================================

#' Generic Function for Extracting Read Coverage
#'
#' A S4 generic function to extract RNA-seq read coverage for a given gene.
#' This function dispatches S4 methods based on the class of `object`.
#'
#' @param object An object containing RNA-seq or ribosome profiling data.
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' This generic function supports **different implementations based on object type**.
#' Specific methods are defined via `setMethod()`.
#'
#' - **`get_coverage,ribotrans-method`**: Extracts per-position RNA-seq read coverage from a `ribotrans` object.
#'
#' @return The return value depends on the method. Generally, a `ribotrans` object is returned with computed coverage.
#'
#' @seealso \code{\link{getCoverage}} for coverage computation from BAM files.
#' @export
setGeneric("get_coverage",function(object,...) standardGeneric("get_coverage"))



#' Extract RNA-Seq Read Coverage for a Gene
#'
#' Retrieves per-position RNA-seq read coverage from a `ribotrans` object.
#' Normalization and optional smoothing are applied.
#'
#' @param object A `ribotrans` object containing transcript annotations, BAM files, and library size.
#' @param gene_name A `character` string specifying the gene of interest.
#' @param smooth A `character` string (`"TRUE"` or `"FALSE"`, default: `"TRUE"`)
#'     indicating whether to apply smoothing.
#' @param slide_window An `integer` specifying the window size used for smoothing (default: `30`).
#'
#' @details
#' The function **retrieves RNA-seq coverage per nucleotide position** using `getCoverage()`.
#' The process involves:
#'
#' 1. Extracting **RNA-seq BAM files** from the `ribotrans` object.
#' 2. Matching **library size (total mapped reads)** for each sample.
#' 3. Iterating over samples to compute per-position coverage.
#' 4. **Optional smoothing** using `smoothEachPosition()`, applying a sliding window.
#'
#' **Smoothing is done with a sliding mean filter over transcript positions.**
#'
#' @return A `ribotrans` object with updated slot:
#' \describe{
#'   \item{RNA_coverage}{A `data.frame` with per-position read coverage.}
#'   \item{RNA.smoothed}{A `character` (`"TRUE"`/`"FALSE"`) denoting if smoothing was applied.}
#'   \item{gene_name}{Updated with `gene_name` argument.}
#' }
#'
#' @importFrom purrr map_df
#' @export
setMethod("get_coverage",
          signature(object = "ribotrans"),
          function(object, gene_name = NULL,
                   smooth = c("TRUE", "FALSE"),
                   slide_window = 30){
            smooth <- match.arg(smooth, c("TRUE", "FALSE"))

            bf <- subset(object@bam_file, type == "rna")
            rnabams <- bf$bam
            gp <- bf$sample

            lib <- subset(object@library, type == "rna")
            lib <- lib[match(rnabams, lib$bam),]
            totalreads <- lib$mappped_reads

            # loop for each sample
            purrr::map_df(seq_along(rnabams),function(x){
              tmp <- getCoverage(bam_file = rnabams[x],
                                 gene_name = gene_name,
                                 features = object@features,
                                 total_reads = totalreads[x])

              # add sample name
              tmp$sample <- gp[x]

              return(tmp)
            }) -> rna.df

            # smooth for each position
            if(smooth == TRUE){
              sm.df <- smoothEachPosition(features = object@features,
                                          posdf = rna.df,
                                          slide_window = slide_window)

              object@RNA_coverage <- sm.df
              object@RNA.smoothed <- "TRUE"
            }else{
              rna.df$smooth <- rna.df$rpm
              object@RNA_coverage <- rna.df
              object@RNA.smoothed <- "FALSE"
            }

            object@gene_name <- gene_name

            return(object)
          }
)



# ==============================================================================
# function to calculate ribosome occupancy
# ==============================================================================

#' Generic Function for Extracting Ribosome Occupancy
#'
#' A S4 generic function for extracting ribosome footprinting occupancy from different object types.
#' It dispatches S4 methods based on the class of `object`.
#'
#' @param object An object containing ribosome profiling (`ribo-seq`) or transcriptomic data.
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' This generic function enables **ribosome occupancy extraction** from different S4 objects.
#' Specific implementations are handled via `setMethod()`.
#'
#' - **`get_occupancy,ribotrans-method`**: Extracts ribosome footprinting per-position occupancy from a `ribotrans` object.
#'
#' @return The output format depends on the method. Generally, it returns a `ribotrans` object with updated occupancy data.
#'
#' @seealso \code{\link{getOccupancy}} for computing occupancy directly from BAM files.
#' @export
setGeneric("get_occupancy",function(object,...) standardGeneric("get_occupancy"))


#' Extract Ribosomal Footprint Occupancy for a Gene
#'
#' Computes per-position ribosome occupancy from a `ribotrans` object.
#' Optionally, a **sliding window smoothing** may be applied.
#'
#' @param object A `ribotrans` object containing ribosome profiling data, BAM files, and total read counts.
#' @param gene_name A `character` string specifying the gene of interest.
#' @param smooth A `character` string (`"TRUE"` or `"FALSE"`, default: `"TRUE"`)
#'     indicating whether to apply smoothing to occupancy values.
#' @param slide_window An `integer` specifying the window size used for smoothing (default: `30`).
#'
#' @details
#' This method extracts **ribosome footprint occupancy per transcript position** using `getOccupancy()`.
#' **Workflow:**
#'
#' 1. Filters **ribo-seq BAM files** associated with the `ribotrans` object.
#' 2. Matches **total mapped reads** for normalization (RPM).
#' 3. Iterates over **all ribo-seq samples**, computing per-position occupancy.
#' 4. **Optional smoothing** using `smoothEachPosition()`, applying a sliding window.
#'
#' **Reads Per Million (RPM) Normalization:**
#' \deqn{RPM = \frac{\text{Read Count at Position}}{\text{Total Mapped Reads}} \times 10^6}
#'
#' @return A `ribotrans` object with updated slots:
#' \describe{
#'   \item{ribo_occupancy}{A `data.frame` with per-position ribosome occupancy.}
#'   \item{ribo.smoothed}{A `character` (`"TRUE"`/`"FALSE"`) marking if smoothing was applied.}
#'   \item{gene_name}{Updated with the queried `gene_name`.}
#' }
#'
#' @importFrom purrr map_df
#' @export
setMethod("get_occupancy",
          signature(object = "ribotrans"),
          function(object, gene_name = NULL,
                   smooth = c("TRUE", "FALSE"),
                   slide_window = 30){
            smooth <- match.arg(smooth, c("TRUE", "FALSE"))

            bf <- subset(object@bam_file, type == "ribo")
            ribobams <- bf$bam
            gp <- bf$sample

            lib <- subset(object@library, type == "ribo")
            lib <- lib[match(ribobams, lib$bam),]
            totalreads <- lib$mappped_reads

            # loop for each sample
            purrr::map_df(seq_along(ribobams),function(x){
              tmp <- getOccupancy(bam_file = ribobams[x],
                                  gene_name = gene_name,
                                  features = object@features,
                                  total_reads = totalreads[x])

              # add sample name
              tmp$sample <- gp[x]

              return(tmp)
            }) -> ribo.df

            # smooth for each position
            if(smooth == TRUE){
              sm.df <- smoothEachPosition(features = object@features,
                                          posdf = ribo.df,
                                          slide_window = slide_window)

              object@ribo_occupancy <- sm.df
              object@ribo.smoothed <- "TRUE"
            }else{
              ribo.df$smooth <- ribo.df$rpm
              object@ribo_occupancy <- ribo.df
              object@ribo.smoothed <- "FALSE"
            }

            object@gene_name <- gene_name

            return(object)
          }
)



# ==============================================================================
# scale to mRNA abundance
# ==============================================================================

#' Generic Function for Computing Scaled Ribosome Occupancy
#'
#' A S4 generic function to compute **scaled ribosome occupancy**
#' (ribo-seq signal normalized to RNA-seq coverage).
#' This function dispatches S4 methods based on the class of `object`.
#'
#' @param object An object containing ribosome profiling (`ribo-seq`) and RNA-seq coverage data.
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' This generic function allows different implementations of **scaled ribosome occupancy calculations**.
#' Specific implementations should be defined using `setMethod()`.
#'
#' - **`get_scaled_occupancy,ribotrans-method`**: Computes **Ribo/RNA enrichment per position** from a `ribotrans` object.
#'
#' @return The return structure depends on the method. Generally, it returns a `ribotrans` object with updated scaled occupancy data.
#'
#' @seealso \code{\link{get_occupancy}}, \code{\link{get_coverage}}
#' @export
setGeneric("get_scaled_occupancy",function(object,...) standardGeneric("get_scaled_occupancy"))


#' Compute Scaled Ribosome Occupancy for a Gene
#'
#' Computes **position-wise ribosome footprint enrichment**,
#' defined as **Ribo-seq RPM / RNA-seq RPM**, and optionally smooths the data.
#'
#' @param object A `ribotrans` object containing ribosome profiling and RNA-seq coverage.
#' @param slide_window An `integer` specifying the smoothing window size (default: `30`).
#'
#' @details
#' This function **computes scaled occupancy by normalizing ribosome occupancy with RNA-seq coverage**.
#' The computation follows the formula:
#'
#' \deqn{Scaled\ Occupancy = \frac{\text{Ribo-seq RPM at Position}}{\text{RNA-seq RPM at Position}}}
#'
#' **Processing steps:**
#' 1. **Ensures RNA-seq (`RNA_coverage`) and ribo-seq (`ribo_occupancy`) data exist**.
#' 2. **Merges RNA-seq and Ribo-seq coverage by transcript position**.
#' 3. **Computes ribo/RNA enrichment per position**.
#' 4. **Performs sliding window smoothing**.
#'
#' **Missing values (`NA` or `Inf`) are replaced with `0`** to handle unmapped positions.
#'
#' @return A `ribotrans` object with an updated slot:
#' \describe{
#'   \item{scaled_occupancy}{A `data.frame` containing position-wise **Ribo-seq / RNA-seq** enrichment.}
#' }
#'
#' @importFrom dplyr select left_join mutate filter
#' @importFrom purrr map_df
#' @importFrom zoo rollmean
#' @export
setMethod("get_scaled_occupancy",
          signature(object = "ribotrans"),
          function(object,slide_window = 30){
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
              dplyr::select(sample,rname,pos,smooth) %>%
              dplyr::left_join(y = object@RNA_coverage[,c("sample","rname","pos","smooth")],
                               by = c("sample","rname","pos")) %>%
              dplyr::mutate(enrich = (smooth.x/smooth.y)) %>%
              dplyr::mutate(enrich = ifelse(is.na(enrich) | is.infinite(enrich),0,enrich))

            # loop for each sample and transcript to smooth
            tanno <- subset(object@features, gene == object@gene_name)

            purrr::map_df(1:nrow(tanno),function(x){
              tmp <- tanno[x,]

              tmp2 <- mb %>%
                dplyr::filter(rname == tmp$idnew)

              sp <- unique(tmp2$sample)

              # loop for each sample
              purrr::map_df(seq_along(sp),function(x){
                tmp3 <- tmp2 %>%
                  dplyr::filter(sample == sp[x])

                # smooth data
                tmp3$smooth <- zoo::rollmean(tmp3$enrich, k = slide_window, fill = 0)

                return(tmp3)
              }) -> smooth.df

              return(smooth.df)
            }) -> mb.smoothed

            object@scaled_occupancy <- mb.smoothed
            return(object)
          }
)




