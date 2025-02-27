#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")
  if (requireNamespace("cli", quietly = TRUE)) {
    start_up_mg <- cli::cat_boxx("Welcome to use riboTransVis package for Ribo-seq analysis.",
                                 col = "#8B1874")
  } else {
    warning("Package 'cli' is needed for this function to work.")
  }

  packageStartupMessage(start_up_mg)
  packageStartupMessage(paste("The version of riboTransVis:",
                              pkgVersion,
                              "\nAny advice or suggestions please contact with me: 3219030654@stu.cpu.edu.cn.",
                              sep = " "))
}
