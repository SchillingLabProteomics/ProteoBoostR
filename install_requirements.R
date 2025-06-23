# install CRAN packages
install.packages(c('shiny','shinydashboard','shinyjs','DT','dplyr','tibble','pROC','ROSE','ggplot2','fontawesome'),
 repos='https://cran.rstudio.com/')

if (!requireNamespace('remotes', quietly = TRUE))
  install.packages('remotes', repos='https://cran.rstudio.com/')

# force specific versions using binary packages if available
remotes::install_version('caret', version='7.0-1', repos='https://cran.rstudio.com/')
remotes::install_version('xgboost', version='1.7.8.1', repos='https://cran.rstudio.com/', type="source")
remotes::install_version('rBayesianOptimization', version='1.2.1', repos='https://cran.rstudio.com/', type="source")
