if (!file.exists('r_pkg_install_success.txt')) {
  options(repos = c(CRAN = 'http://cran.us.r-project.org'))
  dir.create('~/R_library', showWarnings = F)
  .libPaths(c('~/R_library', .libPaths()))

  install.packages('devtools')
  install.packages('rmarkdown')
  install.packages('dplyr')
  install.packages('cowplot')
  install.packages('openxlsx')

  source("https://bioconductor.org/biocLite.R")
  biocLite()
  biocLite(c('DESeq2', 'limma', 'SRAdb', 'RUVSeq'))

  ## Install modified repos
  devtools::install_github('warrenmcg/ALDEx2', ref = 'RELEASE_3_6')
  devtools::install_github('warrenmcg/sleuth', ref = 'speedy_fit')
  devtools::install_github('warrenmcg/polyester')

  ## Install other github packages
  devtools::install_github('pimentel/mamabear')
  devtools::install_github('warrenmcg/sleuthALR')
  devtools::install_github('warrenmcg/absSimSeq')

  message <- "all packages have been installed correctly"
  sink('r_pkg_install_success.txt')
  cat(message)
  sink()
} else {
  message('It appears that all the R packages have been installed')
  message('There is nothing to do')
}
