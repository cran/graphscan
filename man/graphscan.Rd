\name{graphscan-package}
\alias{graphscan-package}
\alias{graphscan}
\docType{package}
\title{Cluster detection with hypothesis free scan statistic}
\description{
Multiple scan statistic with variable window for one dimension data and scan statistic based on connected components in 2D or 3D.
}
\details{
\tabular{ll}{
Package: \tab graphscan\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2016-10-11\cr
License: \tab GPL-2 | GPL-3\cr
Depends: \tab R (>= 3.0.2), snowfall (>= 1.84-6)\cr
Imports: ape, sp, methods, rgl, utils, graphics, grDevices, stats\cr
NeedsCompilation: \tab yes\cr
}

Index:
\tabular{ll}{
cluster \tab Performs cluster analysis on 'graphscan' class object.\cr
events_series \tab A 1D cluster example: events series example.\cr
france_two_clusters \tab A 2D cluster example: France with two clusters.\cr
graphscan-class \tab Class '"graphscan"'\cr
graphscan-package \tab Cluster detection with hypothesis free scan statistic.\cr
graphscan_1d \tab Creates objects of class 'graphscan' using 1D data.\cr
graphscan_nd \tab Creates objects of class 'graphscan' using 2D or 3D data.\cr
barplot \tab Barplot of the 1D clusters lengths.\cr
graphscan_plot \tab Plot clusters localisations or the 1D events distributions.\cr
sample3d \tab A 3D cluster example.\cr
summary \tab Summary for graphscan objects.\cr
}

This package implements a statistical method for detecting clusters in dimensions 1, 2 and 3 proposed by Cucala (2008,2009). 
In 1D, this hypothesis multiple scan statistic with variable window in 1D can detect positive clusters only, negative clusters only, or simultaneously positive or negative clusters. Positive clusters correspond to a particularly high concentration in events, while negative clusters correspond to a particularly low concentration in events. The concentration index of Cucala is based on the properties of the distances between order statistics under the hypothesis of uniform distribution of the events. The 1D functions are adapted to study the repartition of mutations along the genome and detect zones with numerous or few mutations when comparing two aligned DNA sequences. When studying a population, i.e. several DNA sequences, pair comparisons can be created and analysed.
In nD (2D or 3D), the flexible spatial scan test for case event data is implemented to detect only the most significant positive cluster. The candidate clusters are build from connected components. The number of candidate clusters is equal to the number of cases minus one. The concentration index of Cucala is based on the properties of the distances between order statistics under the hypothesis of Poisson distribution of the cases. The analysis for nD data returns both the Cucala and the Kulldorff concentration index.
}

\author{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.
Maintainer: David Abrial <graphscan@clermont.inra.fr>}

\references{
Kulldorff, M. 1997. A spatial scan statistic, Commun. Statist. - Theory Meth., 26, 6, p. 1481-1496.

Cucala, L. 2008. A hypothesis-free multiple scan statistic with variable window, Biometrical Journal, 2, p. 299-310.

Cucala, L. 2009. A flexible spatial scan test for case event data, Computational Statistics and Data Analysis, 53, p. 2843-2850.}

\keyword{spatial}
\keyword{cluster}
\keyword{DNA mutation}

