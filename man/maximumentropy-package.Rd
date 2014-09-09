\name{maximumentropy-package}
\alias{maximumentropy}
\docType{package}
\title{
Multinomial Logistic Regression
}
\description{
This R Package allows to use a multinomial logistic regression with binaries and real-valued features.
}
\usage{
maximumentropy(x, y=NULL, data=NULL, iteration=NULL, verbose=TRUE, 
method=c("L-BFGS-B", "GIS", "CG", "BFGS"), 
addslack=TRUE, normalize=FALSE)
}

\details{
\tabular{ll}{
Package: \tab maximumentropy\cr
Type: \tab Package\cr
Version: \tab 0.4\cr
Date: \tab 2014-09-06\cr
License: \tab GNU GPL\cr
}
}
\author{
Jean-Philippe.Fauconnier <jean.philippe.fauconnier@gmail.com>
}
\references{

}

\keyword{ package }

\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
x <- subset(iris,select=-Species)
maximumentropy(x,y)
}
