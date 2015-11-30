\name{simScore}
\alias{simScore}

\title{ Similarity Score based on PGSEA }
\description{
    Similarity Score based on PGSEA. Either ratio or rank is accepted.
}
\usage{
    simScore(obj, upgene, downgene, is.rank=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{obj}{ a numeric matrix, a data.frame, or an ExpressionSet objec}
    \item{upgene}{ up genes. Names should based on obj}
    \item{downgene}{ down genes. Names should based on obj }
    \item{is.rank}{ obj is rank (TRUE) or ratio (FALSE) }
    \item{...}{ Other parameter to PGSEA }
}

\value{
    a simScore vector
}
\author{ Zhilong Jia}
\examples{
    
data(nbEset)
scorePGSEA <- simScore(nbEset, 
    upgene=rownames(exprs(nbEset))[1:250], 
    downgene=rownames(exprs(nbEset))[1:250], 
    is.rank=FALSE, p.value=1)
    
    
}
\keyword{ simScore }
