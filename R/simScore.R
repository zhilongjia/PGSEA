simScore <- function(obj, upgene, downgene, is.rank=TRUE, ...) {
    
    switch(class(obj),
           matrix = mat <- obj,
           ExpressionSet = mat <- Biobase::exprs(obj),
           data.frame = {
               if(any(!sapply(obj,class)%in%c("numeric","integer")))
                   stop("data frame 'obj' contains non-numeric data")
               mat <- as.matrix(obj)
           },
           stop("argument 'obj' must be a matrix, data.frame, or ExpressionSet 
            object")
    )
    
    upsmc=new("smc",ids=upgene)
    downsmc=new("smc",ids=downgene)
    
    uppgscore <- PGSEA(mat, cl=list(upsmc))
    downpgscore <- PGSEA(mat, cl=list(downsmc))
    
    uppgscore <- PGSEA(mat, cl=list(upsmc), ...)
    downpgscore <- PGSEA(mat, cl=list(downsmc), ...)
    
    # Normalisation factor
    sim.max <- simMax(mat, upgene, downgene, is.rank)
    
    if (is.rank) {
        score <- (downpgscore/sim.max[1] - uppgscore/sim.max[1])/2
    } else {
        score <- (uppgscore/sim.max[1,] - downpgscore/sim.max[2,])/2
    }

    score
}

# Called by simMax
getMax <- function(x, mat, upgene, downgene) {
    upsmc=new("smc",ids=tail(names(sort(mat[,x]) ), length(upgene)) )
    downsmc=new("smc",ids=head(names(sort(mat[,x]) ), length(downgene)) )
    uppgscore <- PGSEA(mat, cl=list(upsmc))
    downpgscore <- PGSEA(mat, cl=list(downsmc))
    c(abs(uppgscore[x]), abs(downpgscore[x]) )
}

# Get the max similar score
simMax <- function(mat, upgene, downgene, is.rank=TRUE) {
    if (isTRUE(is.rank)) {
        # only the legnth of genes affects sim.max
        sim.max <- getMax( 1, mat, downgene, upgene)
    } else {
        sim.max <- sapply(1:ncol(mat),  getMax, mat, upgene, downgene)
    }
    return (sim.max)
    
}



makeExpressionSet <- function(dat, state=colnames(dat)){
    
    dat <- data.matrix(dat)
    pdata <- as.data.frame(state)
    rownames(pdata) <- colnames(dat)
    
    metadata <- data.frame(labelDescription=c("state"), row.names=c("state"))
    phenoData <- new("AnnotatedDataFrame", data=pdata, varMetadata=metadata) 
    dataExp <- ExpressionSet(assayData=dat,  phenoData=phenoData)
    dataExp
}




