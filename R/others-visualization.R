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
           aes(x = log2(.data[[x]] + 1),y = log2(.data[[y]] + 1))) +
    geom_point(color = color,size = point_size) +
    geom_abline(intercept = 0,slope = 1,lty = "dashed") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black")) +
    coord_equal()

  if(requireNamespace("ggpubr", quietly = TRUE)){
    p + ggpubr::stat_cor()
  }else{
    warning("Package 'ggpubr' is needed for this function to work.")
  }

}
