## Messages to be displayed when the user loads configural:
.onAttach <- function(libname, pkgname) {
    crayon_enabled <- requireNamespace("crayon", quietly = TRUE)
    pkg_version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname), fields = "Version")
    if (crayon_enabled) {
        packageStartupMessage(crayon::white("----------------------------------------------------- ", crayon::bold(paste(pkgname, "version", pkg_version)), " --"))
        packageStartupMessage("\nPlease report any bugs to ", crayon::italic("github.com/bwiernik/configural/issues"), "\nor ", crayon::italic("brenton@wiernik.org"))
        packageStartupMessage("\nDevelopers work hard to produce these open-source tools for the R community,\n",
                              "please cite configural when you use it in your research:\n",
                              "  Wiernik, B. M. (2024). \n  configural: Criterion profile analysis (Version ", pkg_version, ") [R package].\n",
                              "  https://cran.r-project.org/package=configural (Original work published 2019)\n\n",
                              "  Wiernik, B. M., Wilmot, M. P., Davison, M. L., & Ones, D. S. (2020).\n",
                              "  Meta-analytic criterion profile analysis.\n  ",
                              crayon::italic("Psychological Methods. "), "https://doi.org/10.1037/met0000305")
    } else {
        packageStartupMessage(paste("----------------------------------------------------- ", pkgname, "version", pkg_version, " --"))
        packageStartupMessage("\nPlease report any bugs to github.com/bwiernik/configural/issues \nor brenton@wiernik.org")
        packageStartupMessage("\nDevelopers work hard to produce these open-source tools for the R community,\n",
                              "please cite configural when you use it in your research:\n",
                              "  Wiernik, B. M. (2020). \n  configural: Multivariate profile analysis (Version ", pkg_version, ") [R package].\n",
                              "  https://cran.r-project.org/package=configural (Original work published 2019)\n\n",
                              "  Wiernik, B. M., Wilmot, M. P., Davison, M. L., & Ones, D. S. (2020).\n",
                              "  Meta-analytic criterion profile analysis.\n",
                              "  Psychological Methods. https://doi.org/10.1037/met0000305")
    }

    packageStartupMessage("\nFind info about configural on the web at ",
                          if (crayon_enabled) crayon::italic("https://wiernik.org") else "https://wiernik.org")

}


#' Retrieve the NEWS file for the configural package
#'
#' @description
#' This function gives a shortcut to the `utils::news(package = "configural")` function and displays configural's NEWS file, which contains version information, outlines additions and changes to the package, and describes other updates.
#'
#' @export
#'
#' @importFrom utils news
#'
#' @encoding UTF-8
#'
#' @examples
#' configural_news()
configural_news <- function(){
     news(package = "configural")
}

