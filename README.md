# FleetSegmentation
This is an R-package created for the segmentation of fishing fleets. It contains a standardized, fishery-based approach to fleet segmentation based on multivariate statistics. It was created and is maintained by Erik Sulanke for the Thünen Insitute of Sea Fisheries in Bremerhaven, Germany. To install the package, you will need the most recent version of R (at least 4.1.1) and R Studio (at least 1.4.1717), a recent version of Rtools (rtools40), and the most recent version of the devtools-package (at least 2.4.2). 
The best way to install the FleetSegmentation package is by running:
library(devtools)
and then
install_github("ESulanke/FleetSegmentation", build_manual = T, build_vignettes = T)
from your console. This will not only install the package but also build the function manuals and vignette in your library. The vignette will give you guidance on how to apply the package to your data set and can be found in the 'doc'-folder of the FleetSegmentation-package in your R-library. If you have all the prerequisites set and the installation still doesn't work, don't hesitate to contact me via erik.sulanke@thuenen.de so we can find a way to fix this.

All the best, 
Erik
