## code to prepare `example_shipdata` dataset goes here
library(usethis)
example_lengthdata <- read.csv("data-raw\\Shipdata_Example.csv",sep=";",dec=",")
use_data(example_lengthdata, overwrite = TRUE)

assemblage <- read.csv("data-raw\\assemblage.csv",sep = ";",dec = ",")
use_data(assemblage, overwrite = TRUE)

example_fleetdata <- read.csv("data-raw\\Fleetdata_Example.csv")
use_data(example_fleetdata, overwrite = TRUE)
