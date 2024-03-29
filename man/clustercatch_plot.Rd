% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FleetSegmentation.R
\name{clustercatch_plot}
\alias{clustercatch_plot}
\title{Plot of total catch of ships in clusters}
\usage{
clustercatch_plot(data, clustering, subset = NULL)
}
\arguments{
\item{data}{The original, untransformed data that was used for the clustering.}

\item{clustering}{The result of the clustering procedure, stored as a data frame.}

\item{subset}{Display only a subset of clusters in plot. Can be a single number or a vector of numbers.}
}
\description{
This is function creates an overview barplot of the total catch of all ships in each cluster.
}
\examples{
data <- example_catchdata
stockdata <- assign_stocks(data=data)
catchdata <- catchdata_transformation(data = stockdata)
clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
clustercatch_plot(data = stockdata,clustering = clustering)
}
\keyword{catch}
\keyword{clustering}
\keyword{plot}
