## code to prepare `example_shipdata` dataset goes here
library(usethis)
example_lengthdata <- read.csv("data-raw\\Shipdata_Example.csv",sep=";",dec=",")
use_data(example_lengthdata, overwrite = TRUE)
