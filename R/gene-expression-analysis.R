# ==============================================================================
# counts reads for rna and rpf
# ==============================================================================

#' @title Get Read Counts for RNA-seq and Ribo-seq
#'
#' @description
#' Uses Rsubread::featureCounts (for genome-mapped data) or Rsamtools::idxstatsBam / Rsamtools::scanBam
#' (for transcriptome-mapped data) to quantify gene-level counts of RNA-seq and Ribo-seq reads.
#'
#' @details
#' 1. If \code{object@mapping_type} is "genome", the function calls
#'    \code{Rsubread::featureCounts()} to summarize alignments for:
#'      \itemize{
#'        \item RNA-seq: exons by default (\code{rna_feature = "exon"}).
#'        \item Ribo-seq: coding sequences (\code{ribo_feature = "CDS"}).
#'      }
#'    The GTF file path is taken from \code{object@gtf_path}, and extra gene attributes
#'    such as \code{gene_name} and \code{gene_biotype} are captured if present.
#'
#' 2. If \code{object@mapping_type} is "transcriptome", the function performs:
#'      \itemize{
#'        \item \code{idxstatsBam()} to get rapid counts for RNA-seq aligned to transcript-level BAMs,
#'          then uses string parsing and grouping to merge counts by gene_name.
#'        \item \code{scanBam()} to extract read positions for Ribo-seq alignments,
#'          then intersects them with known mRNA regions (\code{object@features}),
#'          aggregating by \code{mstart} and \code{mstop}.
#'      }
#'
#' After counting, the function stores two separate list entries in \code{object@counts}:
#' \code{$rpf} for Ribo-seq data and \code{$rna} for RNA-seq data. Each contains
#' \code{$counts}, \code{$annotation}, \code{$targets}, and \code{$stat}, when available.
#'
#' @param object A \code{\link{ribotrans}} object, which contains metadata and paths to BAM files.
#' @param GTF_attrType.extra A character vector specifying extra gene attributes (beyond the default \code{gene_id})
#'   to be included from the GTF file. Typically \code{c("gene_name","gene_biotype")}.
#' @param rna_feature The annotation feature type to use for counting RNA-seq reads
#'   (e.g. \code{"exon"}). Only used if \code{object@mapping_type} is \code{"genome"}.
#' @param rna_PairedEnd A logical indicating if the RNA reads are paired-end (default is FALSE).
#' @param ribo_feature The annotation feature type to use for counting Ribo-seq reads
#'   (e.g. \code{"CDS"}). Only used if \code{object@mapping_type} is \code{"genome"}.
#' @param ribo_PairedEnd A logical indicating if the ribo reads are paired-end (default is FALSE).
#' @param nThreads Number of threads for parallel processing (default is 1).
#' @param ... Additional arguments (currently unused).
#'
#' @return The input \code{object} (of class \code{\link{ribotrans}}) is returned,
#'   with a new slot \code{object@counts} containing count matrices for RNA-seq and
#'   Ribo-seq data. This slot is a list with elements:
#'   \describe{
#'     \item{\code{rpf}}{Counts and summary info for Ribo-seq data.}
#'     \item{\code{rna}}{Counts and summary info for RNA-seq data.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Assuming 'my_ribotrans_obj' is a valid ribotrans object
#'   # and its BAM files plus GTF path have already been specified:
#'   my_ribotrans_obj <- get_counts(
#'     my_ribotrans_obj,
#'     GTF_attrType.extra = c("gene_name","gene_biotype"),
#'     rna_feature = "exon",
#'     ribo_feature = "CDS",
#'     nThreads = 4
#'   )
#'
#'   # Then, you can inspect the counts:
#'   head(my_ribotrans_obj@counts$rna$counts)
#'   head(my_ribotrans_obj@counts$rpf$counts)
#' }
#'
#' @importFrom Rsamtools idxstatsBam scanBam ScanBamParam scanBamFlag
#' @importFrom dplyr filter select full_join n
#' @importFrom fastplyr f_group_by f_summarise f_left_join f_filter
#' @export
#'
#' @rdname get_counts
setGeneric("get_counts",function(object,...) standardGeneric("get_counts"))




#' @rdname get_counts
#' @export
setMethod("get_counts",
          signature(object = "ribotrans"),
          function(object,
                   GTF_attrType.extra = c("gene_name","gene_biotype"),
                   rna_feature = "exon",
                   rna_PairedEnd = FALSE,
                   ribo_feature = "CDS",
                   ribo_PairedEnd = FALSE,
                   nThreads = 1){
            lib <- object@library
            rnainfo <- subset(lib, type == "rna")
            riboinfo <- subset(lib, type == "ribo")

            # ==================================================================
            # check mapping type
            if(object@mapping_type == "genome"){
              # ==================================================================
              # count rna reads
              if(nrow(rnainfo) > 0){

                if (requireNamespace("Rsubread", quietly = TRUE)) {
                  rna <- Rsubread::featureCounts(files = rnainfo$bam,
                                                 isGTFAnnotationFile = TRUE,
                                                 isPairedEnd = rna_PairedEnd,
                                                 GTF.featureType = rna_feature,
                                                 annot.ext = object@gtf_path,
                                                 GTF.attrType = "gene_id",
                                                 GTF.attrType.extra = GTF_attrType.extra,
                                                 nthreads = nThreads)

                  colnames(rna$counts) <- rnainfo$sample
                } else {
                  warning("Package 'Rsubread' is needed for this function to work.")
                }
              }else{
                rna <- data.frame()
              }

              # ==================================================================
              # count ribo reads
              if(nrow(riboinfo) > 0){

                if (requireNamespace("Rsubread", quietly = TRUE)) {
                  rpf <- Rsubread::featureCounts(files = riboinfo$bam,
                                                 isGTFAnnotationFile = TRUE,
                                                 isPairedEnd = ribo_PairedEnd,
                                                 GTF.featureType = ribo_feature,
                                                 annot.ext = object@gtf_path,
                                                 GTF.attrType = "gene_id",
                                                 GTF.attrType.extra = GTF_attrType.extra,
                                                 nthreads = nThreads)

                  colnames(rpf$counts) <- riboinfo$sample
                } else {
                  warning("Package 'Rsubread' is needed for this function to work.")
                }
              }else{
                rpf <- data.frame()
              }

              # return
              counts_info <- list(rpf = rpf,rna = rna)

              object@counts <- counts_info
            }else{
              # ================================================================
              # mapping type is "transcriptome"

              # count rna reads
              if(nrow(rnainfo) > 0){
                rnabms <- rnainfo$bam
                rnasps <- rnainfo$sample

                # x = 1
                lapply(seq_along(rnabms),function(x){

                  tinfo <- idxstatsBam(rnabms[x]) %>%
                    dplyr::filter(seqnames != "*") %>%
                    dplyr::select(seqnames,mapped)

                  tinfo$seqnames <- sapply(strsplit(as.character(tinfo$seqnames),split = "\\|"),"[",2)

                  tinfo <- tinfo %>%
                    fastplyr::f_group_by(seqnames) %>% fastplyr::f_summarise(mapped = round(mean(mapped)))

                  colnames(tinfo) <- c("gene_name",rnasps[x])

                  return(tinfo)
                }) %>% Reduce(function(x, y) dplyr::full_join(x, y, by = "gene_name"), .) -> rna.count

                rna.count <- rna.count %>%
                  tibble::column_to_rownames(var = "gene_name")
                rna.count[is.na(rna.count)] <- 0
              }else{
                rna.count <- data.frame()
              }

              # ================================================================
              # count ribo reads
              if(nrow(riboinfo) > 0){
                ribobms <- riboinfo$bam
                ribosps <- riboinfo$sample

                # x = 1
                lapply(seq_along(ribobms),function(x){

                  bam_data <- Rsamtools::scanBam(file = ribobms[x],
                                                 nThreads = nThreads,
                                                 param = Rsamtools::ScanBamParam(what = c("rname", "pos"),
                                                                                 flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
                  )

                  tinfo <- data.frame(bam_data[[1]])

                  tinfo.anno <- tinfo %>%
                    fastplyr::f_left_join(y = object@features,by = c("rname" = "idnew")) %>%
                    fastplyr::f_filter(pos >= mstart & pos <= mstop) %>%
                    fastplyr::f_group_by(rname) %>%
                    fastplyr::f_summarise(counts = dplyr::n())

                  tinfo.anno$rname <- sapply(strsplit(as.character(tinfo.anno$rname),split = "\\|"),"[",2)

                  tinfo.anno <- tinfo.anno %>%
                    fastplyr::f_group_by(rname) %>% fastplyr::f_summarise(counts = round(mean(counts)))

                  colnames(tinfo.anno) <- c("gene_name",ribosps[x])

                  return(tinfo.anno)
                }) %>% Reduce(function(x, y) dplyr::full_join(x, y, by = "gene_name"), .) -> ribo.count

                ribo.count <- ribo.count %>%
                  tibble::column_to_rownames(var = "gene_name")
                ribo.count[is.na(ribo.count)] <- 0
              }else{
                ribo.count <- data.frame()
              }

              # return
              counts_info <- list(rpf = ribo.count,rna = rna.count)

              object@counts <- counts_info
            }


            return(object)
          }
)



# ==============================================================================
# differential analysis for rna and rpf
# ==============================================================================

#' Gene Differential Expression Analysis
#'
#' Performs differential expression analysis on gene-level count data for either RNA-seq or Ribo-seq
#' experiments stored in a \code{ribotrans} object. The DESeq2 pipeline is used to identify differentially
#' expressed genes between control and treatment sample groups. Results include fold changes, p-values,
#' and gene annotations with significance classification.
#'
#' @param object A \code{ribotrans} object, which should already include count data generated by \code{get_counts()}.
#' @param type Character string specifying the library type to analyze. Either \code{"rna"} or \code{"ribo"}.
#' @param control_samples Character vector of sample names to serve as the control group.
#' @param treat_samples Character vector of sample names to serve as the treatment group.
#' @param lo2FC A numeric vector of length two specifying lower and upper thresholds for log2 fold change.
#' Default is \code{c(-1, 1)}, meaning fold changes below -1 or above 1 are considered significant.
#' @param pval Numeric. P-value threshold for statistical significance. Default is \code{0.05}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @return A \code{data.frame} containing differential expression results, including:
#' \describe{
#'   \item{\code{gene_id}}{Gene identifier}
#'   \item{\code{log2FoldChange}}{Log2-transformed fold change between treatment and control}
#'   \item{\code{pvalue, padj}}{Raw and adjusted p-values (if available)}
#'   \item{\code{gene_name, gene_biotype}}{Gene annotation information}
#'   \item{\code{type}}{Significance label: "sigUp", "sigDown", or "nonSig"}
#' }
#'
#' @details This method uses the \pkg{DESeq2} package for differential gene expression analysis.
#' Gene annotations are automatically extracted from the stored count metadata.
#'
#'
#' @importFrom dplyr left_join mutate case_when filter
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume `obj` is a ribotrans object with valid count data already obtained using get_counts()
#' # and has RNA samples: "ctrl1", "ctrl2", "treat1", "treat2"
#'
#' result <- gene_differential_analysis(
#'   object = obj,
#'   type = "rna",
#'   control_samples = c("ctrl1", "ctrl2"),
#'   treat_samples = c("treat1", "treat2"),
#'   lo2FC = c(-1.5, 1.5),
#'   pval = 0.01
#' )
#'
#' head(result)
#' table(result$type)
#' }
#'
#' @seealso \code{\link{get_counts}}, \code{\link[DESeq2]{DESeq}}, \code{\link[DESeq2]{results}}
setGeneric("gene_differential_analysis",function(object,...) standardGeneric("gene_differential_analysis"))



#' @rdname gene_differential_analysis
#' @export
setMethod("gene_differential_analysis",
          signature(object = "ribotrans"),
          function(object,
                   type = c("rna","ribo"),
                   control_samples = NULL,
                   treat_samples = NULL,
                   lo2FC = c(-1,1),
                   pval = 0.05){
            type <- match.arg(type,choices = c("rna","ribo"))

            # check counts data
            if(length(object@counts) == 0){
              stop("Please run `get_counts` first!")
            }

            # ==================================================================
            # get counts info
            if(type == "rna"){
              if(object@mapping_type == "genome"){
                counts_info <- object@counts$rna

                count_mtx <- counts_info$counts[,c(control_samples,treat_samples)]
              }else{
                count_mtx <- object@counts$rna
              }

            }else{
              if(object@mapping_type == "genome"){
                counts_info <- object@counts$rpf

                count_mtx <- counts_info$counts[,c(control_samples,treat_samples)]
              }else{
                count_mtx <- object@counts$rpf
              }

            }

            # ==================================================================
            # diff analysis

            # gene annotation
            if(object@mapping_type == "genome"){
              gene_anno <- counts_info$annotation[,c("GeneID","Length", "gene_name","gene_biotype")]
              colnames(gene_anno)[1] <- "gene_id"
            }


            # meta data
            coldata <- data.frame(condition = factor(c(rep('control',length(control_samples)),
                                                       rep('treat', length(treat_samples))),
                                                     levels = c('control', 'treat')))

            if (requireNamespace("DESeq2", quietly = TRUE)) {
              # DESeqDataSet
              dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_mtx, colData = coldata, design= ~condition)

              # diff
              dds <- DESeq2::DESeq(dds)

              # get results
              res <- DESeq2::results(dds, contrast = c('condition', 'treat', 'control'))
            } else {
              warning("Package 'DESeq2' is needed for this function to work.")
            }


            # filter NA
            res_all <- as.data.frame(res) %>% dplyr::filter(log2FoldChange != 'NA')

            # add annotation
            if(object@mapping_type == "genome"){
              res_all$gene_id <- rownames(res_all)
              res_all_anno <- res_all %>%
                dplyr::left_join(y = gene_anno,by = 'gene_id')
            }else{
              res_all$gene_name <- rownames(res_all)
              res_all_anno <- res_all
            }


            # add sig type
            res_all_anno <- res_all_anno |>
              dplyr::mutate(type = dplyr::case_when(log2FoldChange >= lo2FC[2] & pvalue < pval ~ "sigUp",
                                                    log2FoldChange <= lo2FC[1] & pvalue < pval ~ "sigDown",
                                                    .default = "nonSig"))

            return(res_all_anno)
          }
)



# ==============================================================================
# calculation normalized matrix and translation efficiency
# ==============================================================================

#' Get normalized reads from ribotrans object
#'
#' This function extracts and calculates normalized read counts (TPM or RPKM)
#' from a ribotrans object for RNA-seq, ribosome profiling, and translation
#' efficiency data.
#'
#' @param object A \code{ribotrans} object containing RNA and ribosome profiling data
#' @param norm_type Character string specifying the normalization method.
#'   Options are "tpm" (Transcripts Per Million) or "rpkm" (Reads Per Kilobase Million).
#'   Default is "tpm".
#' @param type Character vector specifying which data types to extract.
#'   Options include "rna" (RNA-seq data), "ribo" (ribosome profiling data),
#'   and/or "te" (translation efficiency). Default is all three types.
#' @param ... Additional arguments (currently not used)
#'
#' @return A list containing normalized data based on the specified \code{type} parameter:
#' \describe{
#'   \item{tpm.rna/rpkm.rna}{Data frame with normalized RNA-seq counts and gene annotations}
#'   \item{tpm.ribo/rpkm.ribo}{Data frame with normalized ribosome profiling counts and gene annotations}
#'   \item{te}{Data frame with translation efficiency values (ribosome/RNA ratios) and gene annotations}
#' }
#'
#' @details
#' The function performs normalization accounting for both gene length and sequencing depth:
#'
#' \strong{TPM (Transcripts Per Million):}
#' \enumerate{
#'   \item Calculate RPK (Reads Per Kilobase): count/gene_length_kb
#'   \item Calculate TPM: RPK/sum(RPK) * 10^6
#' }
#'
#' \strong{RPKM (Reads Per Kilobase Million):}
#' \enumerate{
#'   \item Normalize by sequencing depth: count/total_reads * 10^6
#'   \item Normalize by gene length: result/gene_length_kb
#' }
#'
#' Translation efficiency (TE) is calculated as the ratio of normalized
#' ribosome profiling reads to normalized RNA-seq reads for each gene.
#'
#' For genome-based mapping, gene lengths are extracted from annotation data.
#' For transcript-based mapping, average transcript lengths are calculated
#' from the features data.
#'
#' @note
#' \itemize{
#'   \item When calculating TE, only genes present in both RNA and ribosome data are included
#'   \item Sample column order is automatically matched between RNA and ribosome data for TE calculation
#'   \item Gene annotations are added based on the mapping type (genome vs transcript)
#' }
#'
#' @examples
#' \dontrun{
#' # Get TPM normalized data for all types
#' normalized_data <- get_normalized_reads(ribotrans_obj,
#'                                        norm_type = "tpm",
#'                                        type = c("rna", "ribo", "te"))
#'
#' # Get only RPKM normalized RNA data
#' rna_rpkm <- get_normalized_reads(ribotrans_obj,
#'                                 norm_type = "rpkm",
#'                                 type = "rna")
#'
#' # Access specific results
#' rna_tpm <- normalized_data$tpm.rna
#' translation_efficiency <- normalized_data$te
#' }
#'
#' @seealso
#' \code{\link[=ribotrans-class]{ribotrans}} for the ribotrans class structure
#'
#' @export
setGeneric("get_normalized_reads",function(object,...) standardGeneric("get_normalized_reads"))




#' @rdname get_normalized_reads
#' @export
setMethod("get_normalized_reads",
          signature(object = "ribotrans"),
          function(object,
                   norm_type = c("tpm","rpkm"),
                   type = c("rna","ribo","te")){
            feature <- object@features
            norm_type <- match.arg(norm_type,choices = c("tpm","rpkm"))

            # check counts data
            if(length(object@counts) == 0){
              stop("Please run `get_counts` first!")
            }

            # ==================================================================
            # get rna tpm
            if("rna" %in% type){
              if(object@mapping_type == "genome"){
                rna_info <- object@counts$rna
                count_rna <- rna_info$counts

                kb <- rna_info$annotation$Length/1000

                if(norm_type == "tpm"){
                  rpk <- count_rna/kb
                  tpm_rna <- t(t(rpk)/colSums(rpk)*10^6)
                }else{
                  tpm_rna <- t(t(count_rna)/colSums(count_rna)*10^6)
                  tpm_rna <- tpm_rna/kb
                }


                # add gene annotation
                gene_anno <- rna_info$annotation[,c("GeneID", "gene_name","gene_biotype")]
                colnames(gene_anno)[1] <- "gene_id"

                tpm_rna.df <- data.frame(tpm_rna,check.names = F)
                tpm_rna.df$gene_id <- rownames(tpm_rna.df)
                tpm_rna.df <- tpm_rna.df %>%
                  dplyr::left_join(y = gene_anno,by = "gene_id")
              }else{
                count_rna <- object@counts$rna

                tmpf <- count_rna %>%
                  tibble::rownames_to_column(var = "gene") %>%
                  fastplyr::f_left_join(y = feature, by = "gene") %>%
                  fastplyr::f_group_by(gene) %>%
                  fastplyr::f_summarise(len = mean(translen))

                kb <- tmpf$len/1000
                rpk <- count_rna/kb
                tpm_rna <- t(t(rpk)/colSums(rpk)*10^6)

                tpm_rna.df <- data.frame(tpm_rna,check.names = F)
                tpm_rna.df$gene_name <- rownames(tpm_rna.df)
              }

            }else{
              tpm_rna.df <- NULL
            }

            # ==================================================================
            # get ribo tpm

            if("ribo" %in% type){
              if(object@mapping_type == "genome"){
                ribo_info <- object@counts$rpf
                count_ribo <- ribo_info$counts

                kb <- ribo_info$annotation$Length/1000
                rpk <- count_ribo/kb
                tpm_ribo <- t(t(rpk)/colSums(rpk)*10^6)

                # add gene annotation
                gene_anno <- ribo_info$annotation[,c("GeneID", "gene_name","gene_biotype")]
                colnames(gene_anno)[1] <- "gene_id"

                tpm_ribo.df <- data.frame(tpm_ribo,check.names = F)
                tpm_ribo.df$gene_id <- rownames(tpm_ribo.df)
                tpm_ribo.df <- tpm_ribo.df %>%
                  dplyr::left_join(y = gene_anno,by = "gene_id")
              }else{
                count_ribo <- object@counts$rpf

                tmpf <- count_ribo %>%
                  tibble::rownames_to_column(var = "gene") %>%
                  fastplyr::f_left_join(y = feature, by = "gene") %>%
                  fastplyr::f_group_by(gene) %>%
                  fastplyr::f_summarise(len = mean(cds))

                kb <- tmpf$len/1000
                rpk <- count_ribo/kb
                tpm_ribo <- t(t(rpk)/colSums(rpk)*10^6)

                tpm_ribo.df <- data.frame(tpm_ribo,check.names = F)
                tpm_ribo.df$gene_name <- rownames(tpm_ribo.df)
              }

            }else{
              tpm_ribo.df <- NULL
            }

            # ==================================================================
            # get te
            if("te" %in% type){
              if(object@mapping_type == "genome"){
                gene_info <- rbind(rna_info$annotation,ribo_info$annotation)
                gene_info <- gene_info[,c("GeneID","gene_name","gene_biotype")] %>% unique()
              }


              # calculation TE(translation efficiency)
              ids <- intersect(rownames(tpm_rna),rownames(tpm_ribo))

              tpm.rna <- tpm_rna[ids,]
              tpm.ribo <- tpm_ribo[ids,]

              # check colnames order to be same
              if(identical(colnames(tpm.rna),colnames(tpm.ribo))){
                te <- tpm.ribo/tpm.rna
              }else{
                tmp.tpm <- tpm.rna[,colnames(tpm.ribo)]
                te <- tpm.ribo/tmp.tpm
              }

              # add gene annotation
              te.df <- data.frame(te,check.names = F)

              if(object@mapping_type == "genome"){
                te.df$gene_id <- rownames(te.df)
                te.df <- te.df %>%
                  dplyr::left_join(y = gene_anno,by = "gene_id")
              }else{
                te.df$gene_name <- rownames(te.df)
              }

            }else{
              te.df <- NULL
            }

            # return
            if(norm_type == "tpm"){
              res <- list(tpm.rna = tpm_rna.df,tpm.ribo = tpm_ribo.df,te = te.df)
            }else{
              res <- list(rpkm.rna = tpm_rna.df,rpkm.ribo = tpm_ribo.df,te = te.df)
            }


          }
)







# ==============================================================================
# TE differential analysis
# ==============================================================================

#' @title Differential Translation Efficiency (TE) Analysis
#'
#' @description
#' This function performs differential translation efficiency analysis using either
#' the riborex or xtail package. It compares two conditions (control vs. treated) at
#' both RNA-seq and Ribo-seq levels, determining changes in translational efficiency
#' and identifying significantly up- or down-regulated genes.
#'
#' @details
#' 1. If \code{pkg = "riborex"}, the function uses \code{\link[riborex]{riborex}}
#'    to assess differential expression on the matched RNA and RPF (ribosome-protected
#'    fragment) counts. This approach can use multiple back-ends such as DESeq2, edgeR,
#'    or voom-limma (as specified in \code{method}) to detect significantly changed transcripts.
#'
#' 2. If \code{pkg = "xtail"}, the function uses xtail to directly compute
#'    changes in translational efficiency (log2FC_TE_final) and associated p-values for each gene.
#'    Genes are then annotated with a simple significance category (sigUp, sigDown, or nonSig)
#'    based on log fold-change thresholds (\code{lo2FC}) and p-value threshold (\code{pval}).
#'
#' 3. The input requires that the \code{ribotrans} object contains:
#'    • RNA-seq and RPF count matrices (in \code{object@counts$rna} and \code{object@counts$rpf})
#'    • Metadata for genome-mapped counts (e.g. \code{object@counts$rna$annotation} and
#'      \code{object@counts$rpf$annotation}) if \code{object@mapping_type == "genome"}
#'
#' 4. The function returns a list containing raw differential expression results as well as
#'    annotation-augmented results, marking significant up/down changes.
#'
#' @param object A \code{ribotrans} object containing count matrices for RNA-seq and RPF,
#'               plus relevant annotation. Must have \code{object@mapping_type} set to
#'               "genome" or "transcriptome".
#' @param pkg A character string, either "riborex" or "xtail", indicating which package
#'            to use for differential TE analysis.
#' @param method Character string specifying the statistical engine to use when
#'               \code{pkg = "riborex"}. Possible values are "DESeq2", "edgeR", "edgeRD", or "Voom".
#' @param control_samples A character vector of sample names that are considered controls
#'                        (must be column names in \code{rna_mtx} and \code{ribo_mtx}).
#' @param treat_samples A character vector of sample names that are considered treated
#'                      (must be column names in \code{rna_mtx} and \code{ribo_mtx}).
#' @param min_count A numeric threshold indicating the minimum mean count required
#'                  for a gene to be included in the analysis (default is 1).
#' @param lo2FC A numeric vector of length 2 specifying the log2 fold change thresholds
#'              for calling significant up/down (default is c(-1, 1)).
#' @param pval A numeric value specifying the p-value cutoff for significance (default is 0.05).
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with two components:
#'   \item{deg_raw}{Raw differential expression results returned by riborex or xtail.}
#'   \item{deg_anno}{Results data frame with columns for gene annotation (if genome-mapped)
#'                   and a "type" column indicating significance status: "sigUp", "sigDown", or "nonSig".}
#'
#'
#'
#' @examples
#' \dontrun{
#'   # Suppose 'my_ribotrans_obj' has been populated with RNA-seq and RPF count data
#'   # for several control and treated samples.
#'   te_results <- TE_differential_analysis(
#'     object = my_ribotrans_obj,
#'     pkg = "riborex",
#'     method = "DESeq2",
#'     control_samples = c("ctrl_1","ctrl_2"),
#'     treat_samples = c("trt_1","trt_2"),
#'     min_count = 10,
#'     lo2FC = c(-1, 1),
#'     pval = 0.01
#'   )
#'
#'   # Check significant upregulated genes
#'   subset(te_results$deg_anno, type == "sigUp")
#' }
#'
#' @export
setGeneric("TE_differential_analysis",function(object,...) standardGeneric("TE_differential_analysis"))




#' @rdname TE_differential_analysis
#' @export
setMethod("TE_differential_analysis",
          signature(object = "ribotrans"),
          function(object,
                   pkg = c("riborex","xtail"),
                   method = c("DESeq2","edgeR","edgeRD","Voom"),
                   control_samples = NULL,
                   treat_samples = NULL,
                   min_count = 1,
                   lo2FC = c(-1,1),
                   pval = 0.05){
            method <- match.arg(method,choices = c("DESeq2","edgeR","edgeRD","Voom"))
            pkg <- match.arg(pkg, c("riborex","xtail"))

            # check counts data
            if(length(object@counts) == 0){
              stop("Please run `get_counts` first!")
            }

            # ==================================================================
            # get counts info
            if(object@mapping_type == "genome"){
              rna_info <- object@counts$rna
              rna_mtx <- rna_info$counts[,c(control_samples,treat_samples)]

              ribo_info <- object@counts$rpf
              ribo_mtx <- ribo_info$counts[,c(control_samples,treat_samples)]

              gene_info <- rbind(rna_info$annotation,ribo_info$annotation)
              gene_info <- gene_info[,c("GeneID","gene_name","gene_biotype")] %>% unique()
              colnames(gene_info)[1] <- "gene_id"
            }else{
              rna_mtx <- object@counts$rna
              ribo_mtx <- object@counts$rpf
            }


            # make row order to be same
            ids <- intersect(rownames(rna_mtx),rownames(ribo_mtx))

            rna_mtx <- rna_mtx[ids,]
            ribo_mtx <- ribo_mtx[ids,]
            # ==================================================================
            # diff analysis
            # make groups
            rnaCond <- c(rep("control",length(control_samples)), rep("treated",length(treat_samples)))
            riboCond <- c(rep("control",length(control_samples)), rep("treated",length(treat_samples)))


            # check which package to use
            if(pkg == "riborex"){
              # deg
              if (requireNamespace("riborex", quietly = TRUE)) {
                deg.res <- riborex::riborex(rnaCntTable = rna_mtx,
                                            riboCntTable = ribo_mtx,
                                            rnaCond = rnaCond,
                                            riboCond = riboCond,
                                            engine = method,
                                            minMeanCount = min_count)
              } else {
                warning("Package 'riborex' is needed for this function to work.")
              }

              # add annotation
              diff_df <- data.frame(deg.res)

              if(object@mapping_type == "genome"){
                diff_df$gene_id <- rownames(diff_df)

                # add gene_name
                diff_df_anno <- diff_df %>% dplyr::inner_join(gene_info,by = 'gene_id')
              }else{
                diff_df$gene_name <- rownames(diff_df)
                diff_df_anno <- diff_df
              }


              # add sig type
              if(method == "DESeq2"){
                res_all_anno <- diff_df_anno |>
                  dplyr::mutate(type = dplyr::case_when(log2FoldChange >= lo2FC[2] & pvalue < pval ~ "sigUp",
                                                        log2FoldChange <= lo2FC[1] & pvalue < pval ~ "sigDown",
                                                        .default = "nonSig"))
              }else if(method %in% c("edgeR","edgeRD")){
                res_all_anno <- diff_df_anno |>
                  dplyr::mutate(type = dplyr::case_when(logFC >= lo2FC[2] & PValue < pval ~ "sigUp",
                                                        logFC <= lo2FC[1] & PValue < pval ~ "sigDown",
                                                        .default = "nonSig"))
              }else if(method == "Voom"){
                res_all_anno <- diff_df_anno |>
                  dplyr::mutate(type = dplyr::case_when(logFC >= lo2FC[2] & P.Value < pval ~ "sigUp",
                                                        logFC <= lo2FC[1] & P.Value < pval ~ "sigDown",
                                                        .default = "nonSig"))
              }
            }else{
              # ================================================================
              # xtail for diff TE analysis
              ft <- apply(rna_mtx,1,function(x){
                if(sum(x>0) >= min_count){
                  return(TRUE)
                }else{
                  return(FALSE)
                }
              })

              mrna <- rna_mtx[ft,]
              rpf <- ribo_mtx[ft,]

              # diff test
              if (requireNamespace("xtail", quietly = TRUE)) {
                deg.res <- xtail::xtail(mrna = mrna,
                                        rpf = rpf,
                                        condition = riboCond)

                res <- data.frame(xtail::resultsTable(deg.res))
              } else {
                warning("Package 'xtail' is needed for this function to work.")
              }


              res <- res |>
                dplyr::mutate(type =  dplyr::case_when(log2FC_TE_final >= lo2FC[2] & pvalue_final < pval ~ "sigUp",
                                                       log2FC_TE_final <= lo2FC[1] & pvalue_final < pval ~ "sigDown",
                                                       .default = "nonSig"))

              if(object@mapping_type == "genome"){
                res$gene_id <- rownames(res)

                # add gene_name
                res_all_anno <- res %>% dplyr::inner_join(gene_info,by = 'gene_id')
              }else{
                res$gene_name <- rownames(res)
                res_all_anno <- res
              }

            }


            # return
            diff.res <- list(deg_raw = deg.res, deg_anno = res_all_anno)
            return(diff.res)
          }
)






