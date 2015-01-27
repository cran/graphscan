\name{graphscan_nd}

\alias{graphscan_nd}
\alias{graphscan_nd,SpatialPoints-method}

\title{Creates objects of class 'graphscan' using 2D or 3D data.}

\description{This function produces objects of class 'graphscan' using 'SpatialPointsDataFrame' class objects for 2D and 3D analysis. }

\usage{
graphscan_nd(data,field_cases = NULL,field_controls = NULL,
n_simulation = 199, alpha = 0.05)
}

\arguments{ 
  \item{data}{a 'SpatialPointsDataFrame' object containing coordinates of cases and controls points. The 'data.frame' of 'data' must contain two numeric fields: one for the number of cases for each point (named by default 'cases') and one for the number of controls for each point (named by default 'controls'). The minimal number of cases to perform the analysis is 2 and 50 for controls.}
  
  \item{field_cases}{a character string to define the name of the field containing the number of cases per points. If a field named 'cases' exists in the 'data.frame' this argument is optional.}
  
  \item{field_controls}{a character string to define the name of the field containing the number of controls per points. If a field named 'controls' exists in the 'data.frame' this argument is optional.}
  
  \item{n_simulation}{number of simulations (default value 199) used to compute the significance of the clusters i.e. the p-values computed with a Monte-Carlo process. The value of 'n_simulation' is stored in the slot 'param' and can be modified by the function 'cluster'.}
  
  \item{alpha}{the threshold of significance (p-value) to keep the candidate clusters. The value of 'alpha' is stored in the slot 'param'.}
 
}

\details{
This function implements a statistical method for detecting clusters in nD (2D or 3D) proposed by Cucala (2009).
This flexible spatial scan test for case event data is implemented to detect the most significant positive cluster. The candidate clusters are build from connected components. There are the number of case minus one candidate clusters. The concentration index of Cucala is based on the properties of the distances between order statistics under the hypothesis of Poisson distribution of the cases. The analysis for nD data returns both the Cucala and the Kulldorff concentration index.
}

\value{'graphscan_nd' returns an object of class 'graphscan' with 3 slots:
\item{param}{this slot contains the informations about the data and the parameters used to perform the analysis.}
\item{data}{this slot contains a list with one item 'x': the 'SpatialPointsDataFrame' of the 'data' argument. }

\item{cluster}{this slot contains the results of the analysis: a list with three items two 'SpatialPointsDataFrame' named 'cluster_nd_cucala' and 'cluster_nd_kulldorff' and a vector of string named 'cluster_nd_description' (see 'cluster' function for more details).}}

\references{
Kulldorff, M. 1997. A spatial scan statistic, Commun. Statist. - Theory Meth., 26, 6, p. 1481-1496.

Cucala, L. 2008. A hypothesis-free multiple scan statistic with variable window, Biometrical Journal, 2, p. 299-310.

Cucala, L. 2009. A flexible spatial scan test for case event data, Computational Statistics and Data Analysis, 53, p. 2843-2850.}

\author{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.
Maintainer: David Abrial <graphscan@clermont.inra.fr>}

\seealso{graphscan_1d, cluster}

\examples{
data(france_two_clusters)
g3<-graphscan_nd(data=france_two_clusters)

}
\keyword{spatial}
\keyword{cluster}
