#' @importFrom utils packageDescription
#' @importFrom cli cat_boxx
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")
  start_up_mg <- cli::cat_boxx("Welcome to use riboTransVis package for Ribo-seq analysis.",
                               col = "#8B1874")
  packageStartupMessage(start_up_mg)
  packageStartupMessage(paste("The version of riboTransVis:",
                              pkgVersion,
                              "\nAny advice or suggestions please contact with me: 3219030654@stu.cpu.edu.cn.",
                              sep = " "))
}
