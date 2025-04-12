# ==============================================================================
# enrichment function
# ==============================================================================


#' Perform GO/KEGG Enrichment or GSEA Analysis
#'
#' This function conducts Gene Ontology (GO) or Gene Set Enrichment Analysis (GSEA) on input differential data
#' using the clusterProfiler package. It can also perform KEGG pathway enrichment or GSEA. The function returns
#' both the enrichment result objects and filtered data frames containing significant terms/pathways.
#'
#' @param diff_data A data frame containing differential expression results with at least the columns \code{gene_name},
#' \code{log2fc_col}, and \code{pval_col}.
#' @param log2fc_col A character string specifying the name of the column in \code{diff_data} that contains log2 fold changes.
#' @param log2fc_threshold A numeric value specifying the log2 fold-change threshold for filtering differentially expressed genes.
#' Default is 1.
#' @param pval_col A character string specifying the name of the column in \code{diff_data} that contains p-values or adjusted p-values.
#' @param pval_threshold A numeric value specifying the threshold for filtering by p-values.
#' Default is 0.05.
#' @param OrgDb An \strong{OrgDb} object (e.g., \code{org.Hs.eg.db}) for gene annotation.
#' @param organism A character string specifying the KEGG organism code (e.g., \code{"hsa"} for Homo sapiens).
#' @param type A character string indicating the analysis type. Can be either \code{"go"} for GO/KEGG enrichment
#' or \code{"gsea"} for Gene Set Enrichment Analysis. Default is \code{c("go","gsea")} (match-arg).
#' @param trend A character string specifying whether to focus on up- or down-regulated genes for GO/KEGG enrichment.
#' Can be \code{"up"} or \code{"down"}. Default is \code{c("up","down")} (match-arg).
#' @param ont A character string specifying the GO sub-ontology. Can be \code{"BP"}, \code{"CC"}, \code{"MF"}, or \code{"ALL"}.
#' Default is \code{"BP"}.
#' @param pvalueCutoff A numeric value specifying the significance cutoff threshold for terms or pathways in the enrichment/GSEA results.
#' Default is 0.05.
#'
#' @return A named list with the following elements:
#' \itemize{
#' \item \code{enrich.go.obj}: The GO enrichment or GSEA result (depending on \code{type}).
#' \item \code{enrich.go.df}: A data frame filtered by \code{pvalueCutoff}, containing the GO terms or GSEA results.
#' \item \code{enrich.kegg.obj}: The KEGG enrichment or GSEA result.
#' \item \code{enrich.kegg.df}: A data frame filtered by \code{pvalueCutoff}, containing the KEGG pathways or GSEA results.
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#'
#' # Assume diff_data is a data frame with columns: gene_name, log2FoldChange, pvalue
#' res <- enrich_analysis(
#' diff_data = diff_data,
#' log2fc_col = "log2FoldChange",
#' pval_col = "pvalue",
#' OrgDb = org.Hs.eg.db,
#' organism = "hsa",
#' type = "go",
#' trend = "up",
#' ont = "BP"
#' )
#'
#' # Inspect the returned objects:
#' head(resenrich.go.df) #' head(resenrich.kegg.df)
#' }
#'
#' @export
enrich_analysis <- function(diff_data = NULL,
                            log2fc_col = NULL,
                            log2fc_threshold = 1,
                            pval_col = NULL,
                            pval_threshold = 0.05,
                            OrgDb = NULL,
                            organism = NULL,
                            type = c("go","gsea"),
                            trend = c("up","down"),
                            ont = c("BP","CC","MF","ALL"),
                            pvalueCutoff = 0.05){
  type <- match.arg(type,choices = c("go","gsea"))
  trend <- match.arg(trend,choices = c("up","down"))
  ont <- match.arg(ont,choices = c("BP","CC","MF","ALL"))
  # ============================================================================
  # process data
  tmp <- diff_data

  if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    # id trans
    id <- clusterProfiler::bitr(geneID = tmp$gene_name,
                                fromType = "SYMBOL",
                                toType = "ENTREZID",
                                OrgDb = OrgDb)

    # merge
    tmp.anno <- id %>%
      dplyr::left_join(y = tmp,by = c("SYMBOL" = "gene_name")) %>%
      dplyr::arrange(dplyr::desc(.data[[log2fc_col]]))

    glist <- tmp.anno[,log2fc_col]
    names(glist) <- tmp.anno$ENTREZID

    if(type == "go"){
      # check trend
      if(trend == "up"){
        glist <- tmp.anno %>%
          dplyr::filter(.data[[log2fc_col]] >= log2fc_threshold & .data[[pval_col]] < pval_threshold)
      }else{
        glist <- tmp.anno %>%
          dplyr::filter(.data[[log2fc_col]] <= -log2fc_threshold & .data[[pval_col]] < pval_threshold)
      }
    }

    # ============================================================================
    # enrichment type
    if(type == "go"){
      # go
      go.res <- clusterProfiler::enrichGO(gene = glist$ENTREZID,
                                          OrgDb = OrgDb,
                                          keyType = "ENTREZID",
                                          ont = ont,
                                          pvalueCutoff = 1,
                                          readable = T)

      go.res.df <- data.frame(go.res) %>% dplyr::filter(pvalue < pvalueCutoff)

      # kegg
      keg.res <- clusterProfiler::enrichKEGG(gene = glist$ENTREZID,
                                             organism = organism,
                                             pvalueCutoff = 1)

      keg.res <- clusterProfiler::setReadable(x = keg.res,OrgDb = OrgDb,keyType = "ENTREZID")

      keg.res.df <- data.frame(keg.res) %>% dplyr::filter(pvalue < pvalueCutoff)
    }else{
      go.res <- clusterProfiler::gseGO(geneList = glist,
                                       ont = ont,
                                       OrgDb = OrgDb,
                                       keyType = "ENTREZID",
                                       pvalueCutoff = 1)

      go.res <- clusterProfiler::setReadable(x = go.res,OrgDb = OrgDb,keyType = "ENTREZID")

      go.res.df <- data.frame(go.res) %>% dplyr::filter(pvalue < pvalueCutoff)

      keg.res <- clusterProfiler::gseKEGG(geneList = glist,
                                          organism = organism,
                                          pvalueCutoff = 1)

      keg.res <- clusterProfiler::setReadable(x = keg.res,OrgDb = OrgDb,keyType = "ENTREZID")

      keg.res.df <- data.frame(keg.res) %>% dplyr::filter(pvalue < pvalueCutoff)
    }

  } else {
    warning("Package 'clusterProfiler' is needed for this function to work.")
  }



  # return
  res <- list(enrich.go.obj = go.res,
              enrich.go.df = go.res.df,
              enrich.kegg.obj = keg.res,
              enrich.kegg.df = keg.res.df)

  return(res)
}
