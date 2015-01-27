\name{sample3d}
\alias{sample3d}
\docType{data}
\title{A 3D cluster example.}
\description{This data set is a 3D example with 1060 points ditributed into a cube. These points are devided into 60 cases and 1000 controls.}
\usage{data(sample3d)}
\format{"sample3d" is of class "SpatialPointsDataFrame"}

\source{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.}

\seealso{graphscan_nd}

\examples{
data(sample3d)
g1<-graphscan_nd(data=sample3d)
}
\keyword{datasets}
