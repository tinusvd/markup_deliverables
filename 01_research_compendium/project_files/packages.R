# Packages
packages <- c("bain", "mvtnorm", "renv")

# Install above-mentioned packages if not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load the packages
invisible(lapply(packages, library, character.only = TRUE))

renv::restore("renv.lock")

rm(packages, installed_packages)

