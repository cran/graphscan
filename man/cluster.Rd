\name{cluster}
\alias{cluster}
\title{Performs cluster analysis on 'graphscan' class object.}
\description{
This function performs cluster detections for both 'graphscan_1d' and 'graphscan_nd' objects. For 'graphscan_nd' objects both Kulldorff and Cucala indices are computed.
}
\usage{cluster(gr, n_simulation = NULL, cluster_analysis = NULL, memory_size = 2000)}

\arguments{
  \item{gr}{an object of class graphscan.}
  
  \item{n_simulation}{number of simulations (default value 199) to compute the p-value. This value is set during the generation of the 'graphscan' object or redefined by this argument.}
  
  \item{cluster_analysis}{type of cluster detection. Possible values are "positive", "negative" or "both". For nD detection only "positive" is possible.}
  \item{memory_size}{memory size (default value 2000) to use for the simulation, in mega-bytes(mb, 10^6 bytes). Possible values are integers > 10, and you should try to adapt it to the memory currently available (in the RAM) of your computer.
  }
}
\details{
This is the main function to run the cluster detection analysis from data and parameters contained in an object of class graphscan. The data come from the slot 'data' and the parameters from the slot 'param' of the 'gr' object. The results of the analysis are saved in the slot 'cluster' of the returned object (see value). The analysis for 1D data returns only the concentration index of Cucala. Analysis for nD data returns both the Cucala and the Kulldorff concentration index.
}
\value{
Returns a 'graphscan' object containing the analysis in the slot 'cluster'.
For 1D data, the slot 'cluster' is a list with three matrices 'cluster_1d_raw', 'cluster_1d' and 'cluster_1d_description'. 'cluster_1d_raw' is a matrix with 7 columns containing the raw results of the analysis. 'xleft' and 'xright' columns indicate the boundaries of each detected cluster. 'index' is the concentration index of Cucala and 'pvalue' is the significance of the cluster. 'positivity' is a boolean indicating if the cluster is positive or negative. 'id_segment' and 'id_serie' are the identifiers respectively of the clusters and the events series. 
'cluster_1d' is a matrix with 9 columns containing treated results. Indeed, some clusters are cut into several pieces called segments because they contain one or more included cluster. This matrix indicates for each cluster the start ('xleft') and the end ('xright') of the non-overlapping segments (in 'cluster_1d_raw' matrix clusters are composed by only one segment potentialy overlapping with other detected clusters). Columns 'index','pvalue','positivity', 'id_segment' and 'id_serie' are the same than in 'cluster_1d_raw' matrix. The column 'n_segment' indicate how many segments compose the cluster and 'length' is size of each segment. 'cluster_1d_description' is a matrix with 4 columns giving general informations on all the clusters. 'n_pos' and 'n_neg' give respectively the number of positive and negative clusters for each events series. 'l_pos' and 'l_neg' are respectively the ratio (in percent) of positive and negative clusters total length.
For nD data, the slot 'cluster' is a list with two 'SpatialPointsDataFrame' named 'cluster_nd_cucala' and 'cluster_nd_kulldorff' and a vector of characters named 'cluster_nd_description'. The two 'SpatialPointsDataFrame' objects contain the points of the significant cluster, the 'index' of concentration of Cucala or Kulldorff, the 'radius' of circles to draw the cluster area, the 'pvalue', the number of cases and controls present in the cluster. The vector 'cluster_nd_description' gives a brief description. 'memory_size' is used to limit the number of parallel simulations, due to the possible big memory consumption of a simulation. Indeed, without limitation it could use swap memory in the hard drive and highly decrease the performance of the algorithm.The nd algorithm uses kd-tree for the calculation, and is parallelized, so it has a theorical O(nlog\enc{²}{2}n) complexity and can work with millions of points up to 15 dimensions on a desktop computer in raisonnable time (for example: 200 simulations on a 2D dataset of 1 million of control points and 30'000 case points has been done in 41 minutes using 3 threads at 2 GHz each).
}
\references{
Cucala, L. 2008. A hypothesis-free multiple scan statistic with variable window, Biometrical Journal, 2, p. 299-310.

Cucala, L. 2009. A flexible spatial scan test for case event data, Computational Statistics and Data Analysis, 53, p. 2843-2850.}

\author{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.
Maintainer: David Abrial <graphscan@clermont.inra.fr>}

\seealso{graphscan_1d, graphscan_nd, plot, barplot}

\examples{
\dontrun{
# 1d example with 2 fasta format files containing each 2 DNA aligned sequences.
dna_file<-list.files(path=system.file("extdata",package="graphscan"),
                     pattern="fna",full.names=TRUE)
g1<-graphscan_1d(data=dna_file)
g1<-cluster(g1)

# 2d example
data(france_two_clusters)
g3<-graphscan_nd(data=france_two_clusters)
g3<-cluster(g3)
graphscan_plot(g3,map=france)
}
}
\keyword{spatial}
\keyword{cluster}
\keyword{DNA mutation}
