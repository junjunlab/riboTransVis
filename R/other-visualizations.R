#' Scatter Plot of Log-Transformed Values with Correlation Coefficient
#'
#' This function generates a scatter plot comparing two numeric variables (e.g., gene expression values)
#' after log2 transformation with pseudocount. It adds a dashed y = x reference line and overlays
#' a Pearson correlation coefficient using "ggpubr::stat_cor()".
#'
#' @param data A data frame containing the variables to be plotted.
#' @param x A character string specifying the name of the column to be used for the x-axis.
#' @param y A character string specifying the name of the column to be used for the y-axis.
#' @param color A string specifying the color of the points. Default is "#074799".
#' @param point_size A numeric value indicating the size of the points. Default is 0.5.
#'
#' @return A ggplot2 object representing the scatter plot.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(A = rnorm(100, 10, 2), B = rnorm(100, 10, 2))
#' cor_plot(df, x = "A", y = "B", color = "steelblue", point_size = 1)
#' }
#'
#' @import ggplot2
#'
#' @export
cor_plot <- function(data = NULL,
                     x = NULL, y = NULL,
                     color = "#074799", point_size = 0.5){
  p <-
    ggplot(data,
           aes(x = .data[[x]] + 1,y = .data[[y]] + 1)) +
    geom_point(color = color,size = point_size) +
    geom_abline(intercept = 0,slope = 1,lty = "dashed") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black")) +
    coord_equal() +
    scale_x_continuous(transform = "log2", guide = "axis_logticks",
                       labels = scales::label_log(base = 2,digits = 1)) +
    scale_y_continuous(transform = "log2", guide = "axis_logticks",
                       labels = scales::label_log(base = 2,digits = 1)) +
    xlab(x) + ylab(y)

  if(requireNamespace("ggpubr", quietly = TRUE)){
    p + ggpubr::stat_cor()
  }else{
    warning("Package 'ggpubr' is needed for this function to work.")
  }

}





#' Volcano Plot for Differential Gene Expression Results
#'
#' Generate a volcano plot to visualize differentially expressed genes based on log2 fold change and p-value. Highlight significantly upregulated and downregulated genes, optionally label top genes or specific marker genes, and display counts of significant genes on the plot.
#'
#' @param diff_data A data frame containing the differential expression results. Must include columns for log2 fold change, p-value, and a categorical 'type' column (e.g., "sigUp", "sigDown", "nonSig") indicating significance status. A 'gene_name' column is also required for labeling.
#' @param log2FC String. The column name representing log2 fold changes. Default is `"log2FoldChange"`.
#' @param pval String. The column name representing the p-value. Default is `"pvalue"`.
#' @param log2FC_threshold Numeric vector of length 2. Thresholds for determining significant up/down regulation based on log2 fold change. Default is `c(-1, 1)`.
#' @param top_gene Integer. Number of most significant upregulated/downregulated genes (based on p-value) to label if `marker_gene` is not provided. Default is `10`.
#' @param marker_gene Character vector. Optional list of gene names to highlight on the plot. If NULL, top genes are selected based on significance and p-value. Default is `NULL`.
#' @param gene_number_size Numeric. Text size for displaying significant gene counts. Default is `4`.
#' @param gene_number_label_pos Numeric vector of length 2. Position (npc coordinates) to place gene count labels. Default is `c(0.95, 0.95)`.
#' @param gene_label_size Numeric. Text size for gene name labels. Default is `3`.
#' @param color Character vector of length 3. Colors used for sigUp, nonSig, and sigDown categories. Default is `c("#AF1740", "grey", "#074799")`.
#'
#' @details
#' The function draws a volcano plot using ggplot2. Points are colored by significance category (`type`), with cutoff thresholds visualized as dashed lines. It also annotates the number of significant genes and optionally labels key genes. Top significant genes are determined based on lowest p-value unless `marker_gene` is specified.
#'
#' @return A `ggplot` object representing the volcano plot.
#'
#' @import ggplot2
#' @importFrom dplyr slice_max
#'
#' @examples
#' \dontrun{
#' vocalno_plot(diff_data = my_diff_data,
#'              top_gene = 10,
#'              marker_gene = c("RPL3", "ACT1"),
#'              gene_label_size = 3.5)
#' }
#'
#' @export
vocalno_plot <- function(diff_data = NULL,
                         log2FC = "log2FoldChange", pval = "pvalue",
                         log2FC_threshold = c(-1,1),
                         top_gene = 10,
                         marker_gene = NULL,
                         gene_number_size = 4,
                         gene_number_label_pos = c(0.95,0.95),
                         gene_label_size = 3,
                         color = c("#AF1740", "grey", "#074799")){
  # ============================================================================
  if(is.null(marker_gene)){
    # filter top genes
    up <- subset(diff_data,type == "sigUp") %>%
      dplyr::slice_max(order_by = -log10(.data[[pval]]),n = top_gene)

    dn <- subset(diff_data,type == "sigDown") %>%
      dplyr::slice_max(order_by = -log10(.data[[pval]]),n = top_gene)
  }else{
    mk <- subset(diff_data,gene_name %in% marker_gene)

    up <- subset(mk,type == "sigUp")
    dn <- subset(mk,type == "sigDown")
  }


  # ============================================================================
  # count sig genes
  upn <- paste("sigUp:",nrow(subset(diff_data,type == "sigUp")))
  dnn <- paste("sigDown:",nrow(subset(diff_data,type == "sigDown")))

  # ============================================================================
  # label pos
  posx <- round(range(diff_data[,log2FC]),digits = 0)
  # ============================================================================
  # plot
  p <-
    ggplot(diff_data) +
    geom_point(aes(x = .data[[log2FC]],y = -log10(.data[[pval]]), color = type)) +
    geom_hline(yintercept = -log10(0.05),lty = "dashed") +
    geom_vline(xintercept = log2FC_threshold,lty = "dashed") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black")) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    scale_color_manual(values = c(sigUp = color[1],nonSig = color[2],sigDown = color[3]))

  if(requireNamespace(c("ggrepel","ggpp"), quietly = TRUE)){
    p +
      # ================ sig gene numbers
      ggpp::geom_text_npc(data = data.frame(),
                          aes(npcx = gene_number_label_pos[1],npcy = gene_number_label_pos[2],label = upn),
                          color = color[1],size = gene_number_size) +
      ggpp::geom_text_npc(data = data.frame(),
                          aes(npcx = 1-gene_number_label_pos[1],npcy = gene_number_label_pos[2],label = dnn),
                          color = color[3],size = gene_number_size) +
      # ================ sig gene labels
      ggrepel::geom_text_repel(data = up,
                               aes(x = .data[[log2FC]],y = -log10(.data[[pval]]),label = gene_name),
                               direction = "y",
                               nudge_x = posx[2] - up[,log2FC],
                               hjust = 1,
                               size = gene_label_size,
                               fontface = "italic",
                               segment.linetype = "dashed",
                               segment.colour = "grey",
                               max.overlaps = Inf) +
      ggrepel::geom_text_repel(data = dn,
                               aes(x = .data[[log2FC]],y = -log10(.data[[pval]]),label = gene_name),
                               direction = "y",
                               nudge_x = -posx[2] - dn[,log2FC],
                               hjust = 0,
                               size = gene_label_size,
                               fontface = "italic",
                               segment.linetype = "dashed",
                               segment.colour = "grey",
                               max.overlaps = Inf)
  }else{
    warning("Package 'ggrepel', 'ggpp' is needed for this function to work.")
  }

}




#' Create ECDF Plot with Optional Boxplot Insert
#'
#' This function creates an empirical cumulative distribution function (ECDF) plot
#' for RNA-seq, Ribosome profiling, or Translation Efficiency (TE) data. It can
#' optionally include an inserted boxplot and performs statistical comparisons
#' between samples.
#'
#' @param normed_data A list containing normalized data with elements:
#'   \itemize{
#'     \item \code{tpm.rna}: TPM values for RNA-seq data
#'     \item \code{tpm.ribo}: TPM values for Ribosome profiling data
#'     \item \code{te}: Translation efficiency values
#'   }
#'   Each element should be a data frame with columns: gene_id, gene_name,
#'   gene_biotype, and sample columns.
#' @param type Character string specifying data type. One of "rna", "ribo", or "te".
#'   Default is c("rna", "ribo", "te").
#' @param samples_selected Character vector of sample names to include in the plot.
#'   If NULL, all samples will be used.
#' @param un_selected_cols Character vector of column names to exclude from
#'   pivoting (typically gene annotation columns like gene_id, gene_name).
#'   Default: c("gene_id", "gene_name", "gene_biotype").
#' @param ref_group Character string specifying the reference group for statistical
#'   comparisons. Should match one of the sample names.
#' @param colors Named character vector of colors for each sample. If NULL,
#'   default ggplot2 colors will be used.
#' @param ecdf_xlim Numeric vector of length 2 specifying x-axis limits for the
#'   ECDF plot. If NULL, automatic limits will be used.
#' @param box_ylim Numeric vector of length 2 specifying y-axis limits for the
#'   boxplot insert. If NULL, automatic limits will be used.
#' @param insert_box Logical indicating whether to insert a boxplot. Default is FALSE.
#' @param x_pos Numeric value specifying the x position of the boxplot insert
#'   (as proportion of plot width). Default is 0.45.
#' @param y_pos Numeric value specifying the y position of the boxplot insert
#'   (as proportion of plot height). Default is 0.4.
#' @param width Numeric value specifying the width of the boxplot insert
#'   (as proportion of plot width). Default is 0.2.
#' @param height Numeric value specifying the height of the boxplot insert
#'   (as proportion of plot height). Default is 0.65.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{plot}: A ggplot2 object containing the ECDF plot (with optional boxplot insert)
#'     \item \code{statistics}: A data frame with statistical comparison results from Wilcoxon tests
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Selects the appropriate data based on the \code{type} parameter
#'   \item Converts data from wide to long format
#'   \item Filters out NA and infinite values
#'   \item Performs pairwise Wilcoxon tests using the specified reference group
#'   \item Creates an ECDF plot with log2-transformed values
#'   \item Optionally adds a boxplot insert using cowplot
#' }
#'
#' @note
#' This function requires the following packages:
#' \itemize{
#'   \item \code{ggplot2}
#'   \item \code{reshape2}
#'   \item \code{dplyr}
#'   \item \code{ggpubr}
#'   \item \code{cowplot} (only if \code{insert_box = TRUE})
#' }
#'
#' @examples
#' \dontrun{
#' # Basic ECDF plot for TE data
#' result <- ecdf_plot(normed_data = my_data,
#'                     type = "te",
#'                     samples_selected = c("sample1", "sample2", "sample3"),
#'                     ref_group = "sample1")
#'
#' # ECDF plot with boxplot insert and custom colors
#' result <- ecdf_plot(normed_data = my_data,
#'                     type = "te",
#'                     samples_selected = c("sample1", "sample2", "sample3"),
#'                     ref_group = "sample1",
#'                     colors = c("sample1" = "blue", "sample2" = "red", "sample3" = "green"),
#'                     insert_box = TRUE,
#'                     x_pos = 0.5, y_pos = 0.3)
#'
#' # Access the plot and statistics
#' print(result$plot)
#' print(result$statistics)
#' }
#'
#' @seealso
#' \code{\link[ggplot2]{stat_ecdf}}, \code{\link[ggpubr]{compare_means}},
#' \code{\link[cowplot]{ggdraw}}
#'
#' @export
ecdf_plot <- function(normed_data = NULL,
                      type = c("rna", "ribo", "te"),
                      samples_selected = NULL,
                      un_selected_cols = c("gene_id", "gene_name", "gene_biotype"),
                      ref_group = NULL,
                      colors = NULL,
                      ecdf_xlim = NULL,
                      box_ylim = NULL,
                      insert_box = FALSE,
                      x_pos = 0.45, y_pos = 0.4,
                      width = 0.2, height = 0.65){
  # check data
  if(type == "rna"){
    df <- normed_data$tpm.rna
    xlabel <- "log2(RNA abundence)"
  }else if(type == "ribo"){
    df <- normed_data$tpm.ribo
    xlabel <- "log2(RPF abundence)"
  }else{
    df <- normed_data$te
    xlabel <- "log2(Translation efficiency)"
  }

  #wide format to long
  df.long <- df %>%
    tidyr::pivot_longer(cols = -dplyr::all_of(un_selected_cols),
                        names_to = "sample",
                        values_to = "value") %>%
    stats::na.omit() %>%
    dplyr::filter(!is.infinite(value))

  # filter samples
  if(!is.null(samples_selected)){
    df.long <- subset(df.long, sample %in% samples_selected)
  }


  # ===========================================================================
  # stastics
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required. Please install it.")
  }

  stcs <- ggpubr::compare_means(data = df.long,
                                formula = value ~ sample,
                                ref.group = ref_group,
                                method = "wilcox.test")

  # ===========================================================================
  # plot

  if(is.null(colors)){
    cols <- NULL
    cols2 <- NULL
  }else{
    cols <- scale_color_manual(values = colors,name = "")
    cols2 <- scale_fill_manual(values = colors,name = "")
  }

  if(is.null(ecdf_xlim)){
    xlims <- NULL
  }else{
    xlims <- xlim(ecdf_xlim)
  }

  if(is.null(box_ylim)){
    ylims <- NULL
  }else{
    ylims <- ylim(box_ylim)
  }

  ecdf <-
    ggplot(df.long) +
    stat_ecdf(aes(x = log2(value), color = sample),linewidth = 0.3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          strip.text = element_text(face = "bold",size = rel(1)),
          axis.text = element_text(colour = "black")) +
    cols + xlims +
    ylab("Cumulative Fraction") +
    xlab(xlabel)

  # ============================================================================
  # insert boxplot
  if(insert_box == TRUE){
    ecdf.box <-
      ggplot(df.long) +
      geom_boxplot(aes(x = sample,y = log2(value),fill = sample),
                   width = 0.6,notch = T,
                   outliers = F,show.legend = F) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            plot.background = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_text(face = "bold",size = rel(1))) +
      cols2 + ylims +
      xlab("") + ylab("log2(TE)")

    # INSERT
    if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Package 'cowplot' is required. Please install it.")
    }

    p <- cowplot::ggdraw(ecdf) +
      cowplot::draw_plot(plot = ecdf.box,
                         x = x_pos,y = y_pos,
                         width = width,height = height,
                         vjust = 0.5,hjust = 0)

    return(list(plot = p,statistics = stcs))
  }else{
    return(list(plot = ecdf,statistics = stcs))
  }
}
