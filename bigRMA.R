bigRMA <- function(celFiles, target='core', background=TRUE,
                   normalize=TRUE, h5fname='resultBigRMA.h5',
                   pkgname){
    library(oligo)
    library(preprocessCore)
    library(rhdf5utils)
    library(foreach)
    if (missing(pkgname)){
        annot <- annotation(read.celfiles(celFiles[1]))
    }else{
        annot <- pkgname
        library(annot, character.only=TRUE)
    }

    probeInfo <- getProbeInfo(annot, field='fid', target=target)
    probeInfo <- probeInfo[order(probeInfo$man_fsetid,
                                 probeInfo$fid),]
    rownames(probeInfo) <- NULL
    
    celBatches <- splitIndicesByNode(celFiles)
    
    container <- hdf5Container(h5fname)
    hdf5AddArray(container, 'PM', c(nrow(probeInfo), length(celFiles)),
                 storage.mode='double')
    hdf5AddArray(container, 'rmaSummaries',
                 c(length(unique(probeInfo$man_fsetid)), length(celFiles)),
                 storage.mode='double')
    
    ## remember to export it to nodes
    ## also export probeInfo, celFiles
    myReadCELS <- function(cels){
        headdetails <- .Call("ReadHeader", as.character(cels[1]),
                             PACKAGE="affyio")
        .Call("read_abatch", cels, FALSE,
              FALSE, FALSE, headdetails[[1]],
              headdetails[[2]], FALSE, PACKAGE="affyio")
    }
    
    stats <- foreach(cels=celBatches) %dopar% {
        fid <- probeInfo$fid
        subBatches <- splitIndicesByLength(cels, ocSamples())
        outStats <- vector('list', length(subBatches))
        i <- 1L
        for (cels2read in subBatches){
            colsInFile <- match(cels2read, celFiles)
            pms <- myReadCELS(cels2read)[fid,]
            if (background)
                pms <- backgroundCorrect(pms, copy=FALSE, verbose=FALSE)
            container$PM[, colsInFile] <- pms
            ## prep normalization
            if (normalize){
                outStats[[i]] <- list(n=ncol(pms),
                                      sum=normalize.quantiles.determine.target(pms)*ncol(pms))
            }else{
                outStats[[i]] <- list(n=NA, sum=NA)
            }
            i <- i+1L
            rm(pms)
        }
        stats <- list(n=sum(sapply(outStats, '[[', 'n')),
                      sum=rowSums(do.call(cbind, lapply(outStats, '[[', 'sum'))))
    }
    
    if (normalize){
        targetD <- rowSums(do.call(cbind, lapply(stats, '[[', 'sum')))
        targetD <- targetD/sum(sapply(stats, '[[', 'n'))
        foreach(cels=celBatches) %dopar% {
            subBatches <- splitIndicesByLength(cels, ocSamples())
            for (cels2read in subBatches){
                colsInFile <- match(cels2read, celFiles)
                container$PM[, colsInFile] <- normalize.quantiles.use.target(container$PM[, colsInFile], targetD)
            }
        }
    }
    
    psets <- unique(probeInfo$man_fsetid)
    ## create rmaSummaries matrix length(psets) x length(cels)
    psetBatches <- splitIndicesByNode(psets)
    foreach(psetsOnNode=psetBatches) %dopar% {
        subBatches <- splitIndicesByLength(psetsOnNode, ocProbesets())
        for (psets2read in subBatches){
            rout <- match(psets2read, psets)
            rin <- which(probeInfo$man_fsetid %in% psets2read)
            container$rmaSummaries[rout,] <- basicRMA(container$PM[rin,],
                                                      probeInfo$man_fsetid[rin], normalize=FALSE,
                                                      background=FALSE, verbose=FALSE)
        }
    }

    dnmsBGNM <- list(probeInfo$man_fsetid, celFiles)
    dnmsRMA <- list(psets, celFiles)
    list(h5container=container, dimnamesProbeLevel=dnmsBGNM, dimnamesSummaries=dnmsRMA)
}
