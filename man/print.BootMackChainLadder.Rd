% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BootMackChainLadder.R
\name{print.BootMackChainLadder}
\alias{print.BootMackChainLadder}
\title{print.BootMackChainLadder}
\usage{
\method{print}{BootMackChainLadder}(x, ...)
}
\arguments{
\item{BMCL}{}
}
\value{
the BMCL object
}
\description{
print.BootMackChainLadder
}
\examples{
data(ABC)
BMCL <- BootMackChainLader(Triangle = ABC, B = 100, distNy = "residuals", seuil = 2)
print(BMCL)
}
