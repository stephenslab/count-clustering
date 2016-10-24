

######   BA data run  ################

rm(list=ls())
data <- get(load("../external_data/GTEX_V6/ba_data_gtex.rda"))

devtools::install_github("kkdey/maptpx")
library(maptpx)

topic_fit <- CountClust::FitGoM(data, K=2, tol=0.1)


