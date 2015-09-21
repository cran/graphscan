\name{france_two_clusters}
\alias{france_two_clusters}
\alias{france}
\docType{data}
\title{A 2D cluster example: France with two clusters.}
\description{This data set is a 2D example with 1000 points ditributed on the French territory: 26 cases and 974 controls. The outline of France is included in object "france".}
\usage{data(france_two_clusters)}
\format{"france_two_clusters" and "france" are of class "SpatialPointsDataFrame"}

\source{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.}

\seealso{graphscan_nd}

\examples{
data(france_two_clusters)
sp::plot(france)
points(france_two_clusters[france_two_clusters$cases==1,],pch=16,col="red")
points(france_two_clusters[france_two_clusters$cases==0,],pch=16,cex=0.2,col="green")
}
\keyword{datasets}
