\name{listModelGenes}
\alias{listModelGenes}
\title{listModelGenes}
\description{genes used in the RF or SSP CMS predictor}
\usage{listModelGenes(method=c("RF","SSP")[1])}
\arguments{
\item{method}{a character string : either RF (Random Forest predictor) or SSP (Single Sample Predictor)}
}
\value{a vector of character strings : entrez gene ids of the genes used in the RF or SSP predictor (depending on the chosen method)
}
\author{Justin Guinney, Rodrigo Dienstmann, Aurelien de Reynies}
\keyword{internal}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
