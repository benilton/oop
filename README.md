oop
===

1. Save bigRMA.R
2. Install rhdf5 - biocLite('rhdf5')
3. Install rhdf5utils - check my repos

# Example

```R
source('bigRMA.R')
library(oligo)

## Load parallel back-end, if wanted
library(doMC)
registerDoMC(4)

## get the CEL files
celFiles = list.celfiles()

## run bigRMA
## note that the 'resultBigRMA.h5' file may be very large
## and its size is linearly associated with the size of the
## dataset (ie, it grows with sample size). Try to use a
## local disk
results = bigRMA(celFiles, pkgname='pd.huex.1.0.st.v2',
                 target='core', h5fname='resultBigRMA.h5')

## The line below assumes you have RAM to keep all the results at once
rmaResults = results$h5container$rmaSummaries[,]
dimnames(rmaResults) = results$dimnamesSummaries
head(rmaResults)

## If RAM is not enough, work with chunks
rows = 1:1000
cols = 1:100
chunk = results$h5container$rmaSummaries[rows, cols]
dmns = mapply('[', results$dimnamesSummaries, list(rows, cols))
dimnames(chunk) = dmns
head(chunk)
```
