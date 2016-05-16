## CMSclassifier

An R package and an example data set are provided to run the CMSclassifier. 

To install the CMSclassifier, install from github using devtools:
```
library(devtools)
install_github("Sage-Bionetworks/CMSclassifier")
```

To run the demo below, you will need a Synapse account and the Synapse R client installed:
```
source('http://depot.sagebase.org/CRAN.R')
pkgInstall("synapseClient")
```

Both classifiers (Random forest and single sample) expect data formatted according to the example data set provided:
```
library(synapseClient)
library(CMSclassifier)

synapseLogin()

sampleData <- read.table(synGet("syn4983432")@filePath, sep="\t",header = TRUE,row.names = 1,check.names=FALSE)

Rfcms <- CMSclassifier::classifyCMS(t(sampleData),method="RF")[[3]]
SScms <- CMSclassifier::classifyCMS(t(sampleData),method="SSP")[[3]]
```


