\name{probesToEntrez}
\alias{probesToEntrez}
\title{probesToEntrez}
\description{aggregate probe level GEP to gene level (entrez id)  GEP }
\usage{probesToEntrez(Exp,Gpl,entrez="Entrez")}
\arguments{
\item{Exp}{a dataframe with log2_scaled GEP data values, samples in columns, genes in rows, rownames corresponding to Entrez IDs  }
\item{Gpl}{a dataframe giving the mapping between probes and Entrez gene Ids, probes correspond to rownames, Entrez Gene Ids are given in the column entrez}
\item{entrez}{a character string, the name of the column of the Entrez Gene Ids}
}
\value{ a dataframe with GEP data - 1 row per entrez ID 
}
\author{Aurelien de Reynies}
\keyword{internal}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
