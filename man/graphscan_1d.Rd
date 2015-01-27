\name{graphscan_1d}
\docType{methods}
\alias{graphscan_1d}
\alias{graphscan_1d,character-method}
\alias{graphscan_1d,DNAbin-method}
\alias{graphscan_1d,list-method}
\title{Creates objects of class 'graphscan' using 1D data.}
\description{This function produces objects of class graphscan_1d.}
\usage{
graphscan_1d(data, format = "fasta", events_series = "all", id = NULL, 
n_simulation = 199, cluster_analysis = "both", 
normalisation_factor = NULL, alpha = 0.05, cluster_user_choice = "positive")
}
\arguments{
  \item{data}{
main argument to define the format of the data. 'data' can be a vector of character string corresponding to files names of aligned DNA sequences. In this case, the format can be precised with argument 'format'. 'data' can be also a list of class 'DNAbin' produced by 'read.dna' function of 'ape' package. In all cases, DNA sequences must be aligned.
Finally 'data' can be a 'list' of numeric vector containing the positions of events. This list is called a series of events. These events are not be necessarily on the same segment and not also necessarily on a [0,1] segment. The argument 'normalisation_factor' allows to fix the upper and lower bounds of each events series.
}
  \item{format}{a character string corresponding to the format of the DNA
          sequences contained in files of argument 'data'. This is used by 'read.dna' function of 'ape' package. Possibles values are "interleaved","sequential","clustal" or "fasta" (default).
}
  \item{events_series}{used if 'data' is a set of files names of aligned DNA sequences or list of class 'DNAbin'. 'events_series' can be a list of the form 'list(A,B)' where 'A' and 'B' corresponding to 2 vectors of sequences identifiants. The crossing AxB product is made to obtain a list of the series of events corresponding to the comparison between each sequence from 'A' to each sequence from 'B'. 'events_series' can be also a character string containing "all", in this case all possible comparison between sequences is made.}

  \item{id}{a character string corresponding to the prefix used to create the names of the events series.}
  
  \item{n_simulation}{number of simulations (default value 199) used to compute the p-values of clusters in a Monte-Carlo process. The value of 'n_simulation' is stored in the slot 'param' and can be modified by the function 'cluster'.}

  \item{cluster_analysis}{a character string corresponding to "positive", "negative" or "both" (default value) to detect respectively only the positives clusters, only the negatives clusters or both  positives and negatives clusters. The value of 'cluster_analysis' is stored in the slot 'param' and can be modified by the function 'cluster'.}
  
  \item{normalisation_factor}{a list of vectors with a size equal to the number of events series. Each vector contains 2 integers: the minimum and the maximum for the events positions of series of events. 
  
  The maximum is the length of the DNA sequences if 'data' argument is a vector of character or an object of class "DNAbin". In these cases, the 'normalisation_factor' is automatically computed by the function 'graphscan_1d'.
  
  If 'data' is a 'list' of numeric vector containing the positions of events the 'normalisation_factor' must be specified as a 'list' containing the upper and lower bounds of each events series. The values of 'normalisation_factor' are stored in the slot 'param'.
}
  \item{alpha}{the threshold of significance (p-value) for keeping the candidate clusters. The value of 'alpha' is stored in the slot 'param'.}
  
  \item{cluster_user_choice}{use if 'cluster_analysis="both"'. 'cluster_user_choice' is a string character corresponding to "positive" (default value), "negative" or "random". If two candidates clusters one positive and one negative have the same p-value this argument indicates how to choose between these 2 clusters. The value of 'cluster_user_choice' is stored in the slot 'param'.}
}
\details{
This function implements a statistical method for detecting clusters in 1D proposed by Cucala (2008). This hypothesis multiple scan statistic with variable window in 1D can detect positive clusters only, negative clusters only, or simultaneously positive or negative clusters. Positive clusters correspond to a particularly high concentration in events, while negative clusters correspond to a particularly low concentration in events. The concentration index of Cucala is based on the properties of the distances between oder statistics under the hypothesis of uniform distribution of the events. The 1D functions are adapted to study the repartition of mutations along the genome and detect zones with numerous or few mutations when comparing two aligned DNA sequences. When studying a population, i.e. several DNA sequences, pair comparisons can be created and analysed.
}
\value{
'graphscan_1d' returns an object of class 'graphscan' with 3 slots:
\item{param}{this slot contains the informations about data and the parameters used to perform the analysis.}
\item{data}{this slot contains data of events series as a list 'x' of numeric vectors. }

\item{cluster}{this slot contains the results of the analysis. For 1d, three matrices 'cluster_1d_raw', 'cluster_1d' and 'cluster_1d_description' (see 'cluster' function for more details).}
}

\references{
Cucala, L. 2008. A hypothesis-free multiple scan statistic with variable window, Biometrical Journal, 2, p. 299-310.

Cucala, L. 2009. A flexible spatial scan test for case event data, Computational Statistics and Data Analysis, 53, p. 2843-2850.}

\author{Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer.
Maintainer: David Abrial <graphscan@clermont.inra.fr>}

\seealso{cluster,graphscan_2d}

\examples{
# example with 2 fasta format files containing each 2 DNA aligned sequences.
# the output object contain 6 events series.
dna_file<-list.files(path=system.file("extdata",package="graphscan"),
 pattern="fna",full.names=TRUE)
g1<-graphscan_1d(data=dna_file)

# to perform only 4 comparisons between DNA sequences 
# 1 vs 3, 1 vs 4, 2 vs 3 and 2 vs 4.
g2<-graphscan_1d(data=dna_file,events_series=list(1:2,3:4))

# example with 'DNABin' class object :
data(woodmouse)
g3<-graphscan_1d(data=woodmouse)

# example with a list of 9 events series 
data(events_series)
g4<-graphscan_1d(events_series,normalisation_factor=normalisation_factor)
}
\keyword{cluster}
\keyword{DNA mutation}
