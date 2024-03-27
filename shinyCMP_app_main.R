#######################################################################
# Define the list of packages to be installed and loaded
packages <- c("latex2exp", "tseries", "oce", "ggplot2", 
              "spatstat", "zoo", "akima", "imager", 
              "jpeg", "dbscan", "OPI", "data.table", 
              "DT", "shinycssloaders", "shinyWidgets", 
              "shiny", "shinyjs", "shinyBS", "shinydashboard", 
              "useful", "shinyFiles", "crch", "survival", 
              "later", "pracma", "matrixStats", "RNiftyReg","this.path")

# Function to check, install, and load packages
check_and_load_packages <- function(pkg){
  # Check if the package is installed
  if (!require(pkg, character.only = TRUE)) {
    # Install the package if it is not installed
    install.packages(pkg, dependencies = TRUE)
    # Load the package after installation
    library(pkg, character.only = TRUE)
  }
}

# Loop through the packages vector to check, install, and load each package
for (pkg in packages) {
  check_and_load_packages(pkg)
}
#######################################################################

#Sets working directory to current file path
setwd(this.path::here())

### THIS WILL CLEAR ALL VARIABLES IN THE WORKSPACE! ### 
rm(list=ls())
# Loads response times from OPI package
data(RtDbUnits)
# Runs the app
runApp("CMP_App/")
