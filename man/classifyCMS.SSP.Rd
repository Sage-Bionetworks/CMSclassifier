\name{classifyCMS.SSP}
\alias{classifyCMS.SSP}
\title{classifyCMS.SSP}
\description{nearest-centroid Single Sample Predictor (SSP) of the Consensus Molecular Subtype (CMS) of colorectal cancer samples, based on log2_scaled GEP   }
\usage{classifyCMS.SSP(Exp,minCor=.15,minDelta=.06)}
\arguments{
\item{Exp}{a dataframe with log2_scaled Gene Expression Profiles (GEP) data values, samples in columns, genes in rows, rownames corresponding to Entrez IDs  }
\item{minCor}{minimal correlation (sample x CMS centroid) for a sample to be classified in a CMS (NB : for a correlation below this threshold the sample will remain unclassified)}
\item{minDelta}{minimal difference between the correlation of a sample to its two nearest centroids, to be classified in a CMS (NB: for a delta of correlation below this threshold the sample will remain unclassified)}
}
\value{a dataframe, samples in rows, columns = correlation of the GEP to the four CMS centroids, nearest CMS, correlation to the nearest CMS, delta of correlation between the 2 nearest CMS, predicted CMS (NA value in case of uncertainty)
}
\author{Aurelien de Reynies}
\keyword{internal}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
