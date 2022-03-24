## code to prepare `example_catchdata` and Mediterranean Stocks dataset goes here
library(usethis)
example_catchdata <- read.csv("data-raw\\Catchdata_Example.csv",sep=";",dec=",")
use_data(example_catchdata, overwrite = TRUE)

MED_stocks <- read.csv("data-raw\\MED_Stocks.csv")
use_data(MED_stocks, overwrite = TRUE)
