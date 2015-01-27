\name{events_series}
\alias{events_series}
\alias{normalisation_factor}
\docType{data}
\title{A 1D cluster example: events series example.}
\description{A set of 9 series of 1000 events generated for the graphscan package.
}
\usage{data(events_series)}
\format{
 Data is a list of 9 vectors of reals and a list named 'normalisation_factor' containing the normalisation factor of each events serie.
}

\source{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.}

\seealso{graphscan_1d}

\examples{
data(events_series)
g1<-graphscan_1d(data=events_series,normalisation_factor=normalisation_factor)
}
\keyword{datasets}

