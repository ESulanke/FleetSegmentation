#'@docType data
#'
#'@title Comprehensive data frame of catch, length and main gear of example ships from the FleetSegmentation package
#'
#'@description A dataset containing the ship length over all in meters of the 75 example ships from the FleetSegmentation package.
#' Created by Erik Sulanke, Thuenen-Institute for Sea Fisheries, Bremerhaven, Germany, in March 2022.
#'
#'@keywords datasets
#'
#'@format A data frame with 296 rows and 6 variables:
#' \describe{
#'   \item{ship_ID}{ship identification number}
#'   \item{shiplength}{length of the ship [m]}
#'   \item{gear}{main gear of the ship}
#'   \item{species}{three-digit code for caught species}
#'   \item{fao_area}{area code where species was caught}
#'   \item{landings}{total catch per vessel, species, and area [kg]}

#' }
"example_fleetdata"
