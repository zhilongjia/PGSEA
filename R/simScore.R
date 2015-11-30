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
    
    uppgscore <- PGSEA(mat, cl=list(upsmc), ...)
    downpgscore <- PGSEA(mat, cl=list(downsmc), ...)
    
    if (is.rank) {
        score <- downpgscore - uppgscore
    } else {
        score <- uppgscore - downpgscore
    }
    score <- score/ifelse(max(abs(score), na.rm=TRUE)>0, max(abs(score), na.rm=TRUE), 1)
    
}
