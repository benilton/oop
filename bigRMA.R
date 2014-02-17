bigRMA <- function(celFiles, target='core', background=TRUE,
                   normalize=TRUE, h5fname='resultBigRMA.h5',
                   pkgname, verbose=TRUE){
    library(oligo)
    library(preprocessCore)
    library(rhdf5utils)
    library(foreach)
    if (verbose)
        message('Loading annotation package.')
    if (missing(pkgname)){
        annot <- annotation(read.celfiles(celFiles[1]))
    }else{
        annot <- pkgname
        library(annot, character.only=TRUE)
    }

    if (verbose)
        message('Obtaining probe information.')
    probeInfo <- getProbeInfo(annot, field='fid', target=target)
    probeInfo <- probeInfo[order(probeInfo$man_fsetid,
                                 probeInfo$fid),]
    rownames(probeInfo) <- NULL
    
    celBatches <- splitIndicesByNode(celFiles)
    nBatches <- length(celBatches)

    if (verbose)
        message('Creating HDF5 container: ', h5fname)
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
    
    stats <- foreach(batchID=1:nBatches) %dopar% {
        if (verbose)
            message('Background correcting (if requested) and preparing target distribution for Batch ', batchID)
        cels <- celBatches[[batchID]]
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
        if (verbose)
            message('Background correcting (if requested) and preparing target distribution for Batch ', batchID, '. Done.')
    }
    
    if (normalize){
        if (verbose)
            message('Determining final target distribution. ', appendLF=FALSE)
        targetD <- rowSums(do.call(cbind, lapply(stats, '[[', 'sum')))
        targetD <- targetD/sum(sapply(stats, '[[', 'n'))
        if (verbose)
            message('Done.')
        foreach(batchID=1:nBatches) %dopar% {
            if (verbose)
                message('Normalizing Batch ', batchID)
            cels <- celBatches[[batchID]]
            subBatches <- splitIndicesByLength(cels, ocSamples())
            for (cels2read in subBatches){
                colsInFile <- match(cels2read, celFiles)
                container$PM[, colsInFile] <- normalize.quantiles.use.target(container$PM[, colsInFile], targetD)
            }
            if (verbose)
                message('Normalizing Batch ', batchID, '. Done.')
        }
    }
    
    psets <- unique(probeInfo$man_fsetid)
    ## create rmaSummaries matrix length(psets) x length(cels)
    psetBatches <- splitIndicesByNode(psets)
    nPsetBatches <- length(psetBatches)
    foreach(psetBID=1:nPsetBatches) %dopar% {
        if (verbose)
            message('Running median-polish on Batch ', psetBID)
        psetsOnNode <- psetBatches[[psetBID]]
        subBatches <- splitIndicesByLength(psetsOnNode, ocProbesets())
        for (psets2read in subBatches){
            rout <- match(psets2read, psets)
            rin <- which(probeInfo$man_fsetid %in% psets2read)
            container$rmaSummaries[rout,] <- basicRMA(container$PM[rin,],
                                                      probeInfo$man_fsetid[rin], normalize=FALSE,
                                                      background=FALSE, verbose=FALSE)
        }
        if (verbose)
            message('Running median-polish on Batch ', psetBID, '. Done.')
    }

    dnmsBGNM <- list(probeInfo$man_fsetid, celFiles)
    dnmsRMA <- list(psets, celFiles)
    list(h5container=container, dimnamesProbeLevel=dnmsBGNM, dimnamesSummaries=dnmsRMA)
}
