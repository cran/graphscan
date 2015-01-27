\name{graphscan-class}
\Rdversion{1.1}
\docType{class}

\alias{graphscan-class}
\alias{cluster,graphscan-method}
\alias{print,graphscan-method}
\alias{show,graphscan-method}

\title{Class \code{"graphscan"}}
\description{Class of graphscan objects used to store data, parameters and results of cluster analysis with the method of Cucala.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{"graphscan_1d(data, ...)"} or \code{"graphscan_nd(data, ...)"}.}

\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"list"}, containing the data as events series.}
    \item{\code{param}:}{Object of class \code{"list"}, containing the parameters used in the analysis.}
    \item{\code{cluster}:}{Object of class \code{"list"}, containing the results.}
  }
}

\section{Methods}{
  \describe{
    \item{cluster}{\code{signature(gr = "graphscan")}: performed the cluster analysis. }
    \item{barplot}{\code{signature(x = "graphscan")}: draw an barplot of the length of clusters in multiple events series analysis. }
    \item{graphscan_plot}{\code{signature(x = "graphscan")}: plot clusters localisations or the 1D events distributions.. }
    \item{print}{\code{signature(x = "graphscan")}: print informations about a graphscan object.}
    \item{show}{\code{signature(object = "graphscan")}: print informations about a graphscan object. }
    \item{summary}{\code{signature(object = "graphscan")}: summary for graphscan objects. }
	 }
}
\references{
Kulldorff, M. 1997. A spatial scan statistic, Commun. Statist. - Theory Meth., 26, 6, p. 1481-1496.

Cucala, L. 2008. A hypothesis-free multiple scan statistic with variable window, Biometrical Journal, 2, p. 299-310.

Cucala, L. 2009. A flexible spatial scan test for case event data, Computational Statistics and Data Analysis, 53, p. 2843-2850.}

\author{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
Maintainer: David Abrial <graphscan@clermont.inra.fr>}


\seealso{graphscan_1d, graphscan_nd}
\examples{
\dontrun{
data(events_series)
g4<-graphscan_1d(events_series,normalisation_factor=normalisation_factor)
g4<-cluster(g4)
graphscan_plot(g4)
}
}
\keyword{spatial}
\keyword{cluster}
\keyword{DNA mutation}
