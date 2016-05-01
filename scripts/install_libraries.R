# adapted from http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
# adapted from @roryk https://github.com/roryk/bcbio.rnaseq/blob/master/resources/scripts/install_libraries.R
library(methods)
mirror = "http://cran.at.r-project.org"
update.packages(checkBuilt = TRUE, ask = FALSE, repos=mirror)
is_package_installed = function(package) {
    available = suppressMessages(suppressWarnings(sapply(package, require,
        quietly = TRUE, character.only = TRUE, warn.conflicts = FALSE)))
    missing = package[!available]
    if (length(missing) > 0) return(FALSE)
    return(TRUE)
}

install.cran = function(packages, default_mirror="http://cran.at.r-project.org") {
  for(i in packages) {
    #  require returns TRUE invisibly if it was able to load package
    if(!is_package_installed(i)) {
      options(repos = c(CRAN = default_mirror))
      suppressMessages(suppressWarnings(install.packages(i, dependencies = TRUE, repo=default_mirror)))
    }
  }
}

install.bioconductor = function(packages, default_mirror="http://cran.at.r-project.org") {
  for(i in packages) {
    #  require returns TRUE invisibly if it was able to load package
    if(!is_package_installed(i)) {
      suppressMessages(suppressWarnings(biocLite(i, ask=FALSE)))
    }
  }
}

install.github = function(package,  branch="master") {
    library(devtools)
    for(i in packages) {
        suppressMessages(suppressWarnings(install_github(package, ref=branch)))
    }
}

cran_packages = c("devtools", "ggplot2", "reshape", "dplyr", "knitr", "gplots", "gridExtra", "pheatmap", "gtools")
install.cran(cran_packages)

bioconductor_packages = c("edgeR", "HTSFilter", "DESeq2", "DEGreport)
source("http://bioconductor.org/biocLite.R")
install.bioconductor(bioconductor_packages)

# novel packages
install.github("hbc/CHBUtils")
install.github("lpantano/isomiRs", "master")
install_github('rstudio/rmarkdown')


