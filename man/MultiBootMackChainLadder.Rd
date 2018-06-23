% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BootMackChainLadder.R
\name{MultiBootMackChainLadder}
\alias{MultiBootMackChainLadder}
\title{MultiBootMackChainLadder}
\usage{
MultiBootMackChainLadder(triangles, B = 100, distNy = "normal",
  names = NULL, seuil = NA)
}
\arguments{
\item{triangles}{A List of Triangles objects of the same size.}

\item{B}{The numebr of boostrap replicates}

\item{distNy}{The process distribution for next year increments. Either "normal" (default), "residuals.bycolumn","residuals.global" or "residuals". See details.}

\item{names}{enventual names of the different triangles. IF set to NULL, the procedure will try to get names from the triangles list.}

\item{seuil}{Eventual exclusions limit for residuals. Set to NA (default) to avoid excluding anything.}
}
\value{
a MBMCL object containing a list of BMCL objects and a little more.
}
\description{
The multi boot mack chain ladder algorythme computes todays and next-year common estimates on a portefolio of several triangles, following closely a synchronised version of BootMackChainLadder.
}
\details{
This model uses the fact that the Mack model can be seen as a quasi-glm to found nice residuals. Bootstrap thoose residuals on the upper-left triangle allows to get bootstrap distribution of today's estimatins ( reserves, ultimates, ...).

In each bootstrap resample, the function use the specified process distrubiton to simulate next-year payments. If a normal law is used, it follows boumezoued et all and converges to the Mer-wuthrich fromula in the Braun model. If set to residuals, the convergence is there but not to the same result since no specific law is supposed for residuals.
}
\examples{
data(ABC)
triangles <- list(tri1 = ABC, tri2 = ABC, tri3 = ABC)
MultiBootmackChainLadder(triangles,100)
}
\seealso{
BootMackChainLadder
}