% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BootMackChainLadder.R
\name{mean.MultiBootMackChainLadder}
\alias{mean.MultiBootMackChainLadder}
\title{mean.MultiBootMackChainLadder}
\usage{
\method{mean}{MultiBootMackChainLadder}(x, ...)
}
\arguments{
\item{x}{A MultiBootMackChainLadder Object}
}
\value{
Data frames containing mean informations.
}
\description{
mean.MultiBootMackChainLadder
}
\details{
Return mean informations about the MBMCL computed model. Gives the mean values of estimate quantities in a MBMCL bootstrap.
}
\examples{
data(ABC)
triangles <- list(tri1 = ABC, tri2 = ABC, tri3 = ABC)
MBMCL <- MultiBootmackChainLadder(triangles,100)
mean(MBMCL)
}
\seealso{
MultiBootMackChainLadder, mean.BootMackChainLadder
}
