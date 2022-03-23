#'@docType data
#'
#'@title Assemblage table of Northeast-Atlantic Species
#'
#'@description A dataset containing the stock identifiers and corresponding assemblages for fish and shellfish species of the Northeast Atlantic
#' Created by Erik Sulanke, Thuenen-Institute for Sea Fisheries, Bremerhaven, Germany, in March 2021. Based on species information of the German Federal Administration for Fisheries (BLE) and ecological information of fishbase.org and marinespecies.org
#'
#'@keywords datasets
#'
#'@format A data frame with 2840 rows and 8 variables:
#' \describe{
#'  \item{BLE_fischart}{(German three-letter ICES species code}
#'  \item{species_code}{German three-letter ICES species code}
#'  \item{german_name}{German trivial name of species}
#'  \item{species_name}{English trivial name of species}
#'  \item{species_namesc}{Scientific name of species}
#'  \item{target_assemblage_code}{Three-letter code of corresponding assemblage}
#'  \item{grouped_species}{Comprehensive species grouping}
#'  \item{target_assemblage}{Target assemblage written out}
#' }
"assemblage"

