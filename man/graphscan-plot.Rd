\name{graphscan_plot}

\alias{graphscan_plot}
\alias{graphscan_plot,graphscan-method}

\title{Plot clusters localisations or the 1D events distributions.}

\description{plot clusters localisations in a 1D, 2D or 3D space. In 1D, the distribution (in frequencies) of number of events on each position can be ploted.}

\usage{
graphscan_plot(x,events_series=1,map=NULL,indice="cucala",sphere=TRUE,
	projection=FALSE,\dots)
}

\arguments{
  \item{x}{a graphscan object containing a cluster analysis.}
  \item{events_series}{a numeric or character vector containing 1D cluster identifiants of the events series to draw. If 'events_series="all"' the distribution of number of clusters on each position is plotted.}
  \item{map}{a 'SpatialPolygons' or 'SpatialPolygonsDataFrame' object to add to the 2D graph of clusters localisations (generally the outline of the studied region).}
  \item{indice}{a character string used in nD, to define the type of index to draw. Possible values are "cucala" (default) and "kulldorff".}
  \item{sphere}{a boolean ("TRUE" by default) to define if the spheres used to represent the 3D envelope of the cluster are drawn.}
  \item{projection}{a boolean ("FALSE" by default) to draw a projection in 2D of a 3D cluster. Three plots are drawn respectively for 'y vs x', 'z vs x' and 'z vs y'. }
  \item{...}{further arguments passed to or from other methods.} 
}
\details{To draw the distribution of number of events on each position in 1D, the events series must be of same length and aligned.
The 3D representation of cluster use 'OpenGL' (http://www.opengl.org) an environment for interactive 2D and 3D graphics. If the number of cases points is very important old computers will display graphics quite slowly. In this case, use the option "projection=TRUE".}

\references{
Cucala, L. 2008. A hypothesis-free multiple scan statistic with variable window, Biometrical Journal, 2, p. 299-310.

Cucala, L. 2009. A flexible spatial scan test for case event data, Computational Statistics and Data Analysis, 53, p. 2843-2850.}

\author{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.
Maintainer: David Abrial <graphscan@clermont.inra.fr>}

\seealso{barplot}

\examples{# 1D example:
	\dontrun{
		data(woodmouse)
		g1<-graphscan_1d(data=woodmouse)
		g1<-cluster(g1)
		graphscan_plot(g1,events_series=3)
		dev.new()
		graphscan_plot(g1,events_series="all")

		# 2D example:
		data(france_two_clusters)
		g2<-graphscan_nd(data=france_two_clusters)
		g2<-cluster(g2)
		graphscan_plot(g2,map=france)

		# 3D example:
		data(sample3d)
		g3<-graphscan_nd(data=sample3d)
		g3<-cluster(g3)
		graphscan_plot(g3,projection=TRUE) # 2D plot
		graphscan_plot(g3) # 3D plot
	}
}

\keyword{spatial}
\keyword{cluster}
\keyword{DNA mutation}
