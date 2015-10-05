\name{classifyCMS.RF}
\alias{classifyCMS.RF}
\title{classifyCMS.RF}
\description{Random Forest (RF) predictor of the Consensus Molecular Subtype (CMS) of colorectal cancer samples, based on log2_scaled GEP   }
\usage{classifyCMS.RF(Exp,center=TRUE,minPosterior=.5)}
\arguments{
\item{Exp}{a dataframe with log2_scaled Gene Expression Profiles (GEP) data values, samples in columns (n>1 mandatory), genes in rows, rownames corresponding to Entrez IDs  }
\item{center}{ boolean : should GEP be row-centered (default : TRUE)}
\item{minPosterior}{ numeric : minimal posterior probablity to classify a sample in a CMS}
}
\value{a dataframe, samples in rows, columns = posterior probability to be classified in each of the four CMS centroids, nearest CMS (ie CMS with highest posterior prob), predicted CMS 
}
\author{Justin Guinney, Aurelien de Reynies}
\keyword{internal}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
