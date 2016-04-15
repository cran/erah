\name{setAlPar}
\alias{setAlPar}
\title{Set Alignment Parameters}
\description{Setting alignment parameters for eRah.}
\usage{
setAlPar(min.spectra.cor, max.time.dist,
mz.range=c(70:600))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{min.spectra.cor}{Minimum spectral correlation value. From 0 (non similar) to 1 (very similar). This value sets how similar two or more compounds have be to be considered for alignment between them.}
  \item{max.time.dist}{Maximum retention time distance. This value (in seconds) sets how far two or more compounds can be to be considered for alignment between them.}
  \item{mz.range}{The range of masses that is considered when comparing spectra.}
}
\references{
[1] eRah: an R package for spectral deconvolution, alignment, and metabolite identification in GC/MS-based untargeted metabolomics. Xavier Domingo-Almenara, Alexandre Perera, Maria Vinaixa, Sara Samino, Xavier Correig, Jesus Brezmes, Oscar Yanes. (2016) Article in Press.

}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{newExp}}, \code{\link{setDecPar}}, \code{\link{alignComp}}
}
\examples{

# The following will set eRah for aligning compounds which are
# at least 90 (per cent) similar, and which peaks are at a 
# maximum distance of 2 seconds. All the masses are considered when
# computing the spectral similarity.

# ex.al.par <- setAlPar(min.spectra.cor=0.90, max.time.dist=2,
# mz.range=1:600)

}


