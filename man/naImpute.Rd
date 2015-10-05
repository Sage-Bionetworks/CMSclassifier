\name{naImpute}
\alias{naImpute}
\title{naImpute}
\description{impute NA values in a dataframe}
\usage{naImpute(Exp,G=NULL)}
\arguments{
\item{Exp}{ a dataframe of GEP measures, samples in columns, genes in rows, rownames corresponding to Entrez IDs}
\item{G}{a vector of character strings (optionnal) : entrez gene ids that must be found in the rownames of Exp (will be added if absent, imputed value = global mean)  }
}
\value{the imputed datataframe
}
\author{Aurelien de Reynies}
\keyword{internal}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
