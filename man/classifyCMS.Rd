\name{classifyCMS}
\alias{classifyCMS}
\title{classifyCMS}
\description{predict the Consensus Molecular Subtype (CMS) of colorectal cancer samples, based on log2_scaled GEP   }
\usage{classifyCMS(Exp,method=c("RF","SSP"))}
\arguments{
\item{Exp}{a dataframe with log2_scaled Gene Expression Profiles (GEP) data values, samples in columns, genes in rows, rownames corresponding to Entrez IDs  }
\item{method}{a character vector, accepted values are: RF = random forest  predictor (here data will be automatically row-centered) ; SSP = single sample predictor , based on correlation to centroids (the data won't be row-centered)}
}
\value{a dataframe, samples in rows, columns : output of the predictor(s)
}
\author{Aurelien de Reynies}
\keyword{internal}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
