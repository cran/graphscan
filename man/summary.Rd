\name{summary}

\alias{summary}
\alias{summary,graphscan-method}

\title{Summarize a graphscan object.}

\description{This function give a summary of a graphscan object.}

\usage{
summary(object,...)
}

\arguments{
    \item{object}{a graphscan object containing or not a cluster analysis.}
    \item{...}{arguments to be passed to methods.}
}

\references{
Cucala, L. 2008. A hypothesis-free multiple scan statistic with variable window, Biometrical Journal, 2, p. 299-310.

Cucala, L. 2009. A flexible spatial scan test for case event data, Computational Statistics and Data Analysis, 53, p. 2843-2850.}

\author{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.
Maintainer: David Abrial <graphscan@clermont.inra.fr>}

\examples{
\dontrun{
dna_file<-list.files(path=system.file("extdata",package="graphscan"),
 pattern="fna",full.names=TRUE)
g1<-graphscan_1d(data=dna_file)
summary(g1)
}
}
\keyword{spatial}
\keyword{cluster}
\keyword{DNA mutation}
