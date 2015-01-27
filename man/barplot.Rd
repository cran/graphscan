\name{barplot}

\alias{barplot}
\alias{barplot,graphscan-method}

\title{Barplot of the 1D clusters lengths.}

\description{This function draws a barplot of the total length of the 1D clusters for each events serie among several events series of same length and aligned. The barplot is expressed in frequencies of different percent ranges of the event series length detected as positive, resp. negative, cluster.}

\usage{
barplot(height,...)
}

\arguments{
    \item{height}{a graphscan object containing 1D cluster analysis.}
    \item{...}{arguments to be passed to methods.}
}

\references{
Cucala, L. 2008. A hypothesis-free multiple scan statistic with variable window, Biometrical Journal, 2, p. 299-310.

Cucala, L. 2009. A flexible spatial scan test for case event data, Computational Statistics and Data Analysis, 53, p. 2843-2850.}

\author{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.
Maintainer: David Abrial <graphscan@clermont.inra.fr>}

\seealso{graphscan_plot}

\examples{
\dontrun{
dna_file<-list.files(path=system.file("extdata",package="graphscan"),
 pattern="fna",full.names=TRUE)
g1<-graphscan_1d(data=dna_file)
g1<-cluster(g1)
barplot(g1)
}
}

\keyword{spatial}
\keyword{cluster}
\keyword{DNA mutation}
