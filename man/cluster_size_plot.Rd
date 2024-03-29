% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FleetSegmentation.R
\name{cluster_size_plot}
\alias{cluster_size_plot}
\title{Plot number of ships in clusters}
\usage{
cluster_size_plot(clustering, subset = NULL)
}
\arguments{
\item{clustering}{The result of the clustering procedure, stored as a data frame.}

\item{subset}{Display only a subset of clusters in plot. Can be a single number or a vector of numbers.}
}
\description{
This is function creates an overview plot of the number of ships in each cluster.
}
\examples{
library(tidyverse)
data <- example_catchdata
stockdata <- assign_stocks(data=data)
catchdata <- catchdata_transformation(data = stockdata)
clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
cluster_size_plot(clustering = clustering)
}
\keyword{clustering}
\keyword{plot}
\keyword{size}
