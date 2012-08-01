oop
===

1. Save bigRMA.R
2. Install rhdf5 - biocLite('rhdf5')
3. Install rhdf5utils - check my repos

# Example

```R
source('bigRMA.R')
library(oligo)
celFiles = list.celfiles()
results = bigRMA(celFiles, pkgname='pd.huex.1.0.st.v2', target='core')

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
