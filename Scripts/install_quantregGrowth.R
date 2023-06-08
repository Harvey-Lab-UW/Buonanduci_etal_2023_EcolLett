# Run this code to install the appropriate version of the quantregGrowth package

# First, install current version of the quantreg package
install.packages("quantreg")

# Next, install devtools (if not already installed)
install.packages("devtools")

# Finally, install version 1.4-0 of the quantregGrowth package
devtools::install_version("quantregGrowth", version = "1.4-0", repos = "https://cran.r-project.org")
