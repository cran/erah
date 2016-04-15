\name{setDecPar}
\alias{setDecPar}
\title{Set Software Parameters}
\description{Sets Software Parameters for eRah.}
\usage{
setDecPar(min.peak.width, min.peak.height=500, 
noise.threshold=500,
avoid.processing.mz=c(73:75,147:149), 
compression.coef=2, analysis.time=0)
}
\arguments{
  \item{min.peak.width}{Minimum compound peak widht (in seconds). This is a critical parameter that conditions the efficiency of eRah. Typically, this should be the half of the mean compound width.}
  \item{min.peak.height}{Minimum compound peak height}
  \item{noise.threshold}{Data above this threshold will be considered as noise}
  \item{avoid.processing.mz}{The masses that do not want to be considered for processing. Typically, in GC-MS those masses are 73,74,75,147,148 and 149, since they are they are ubiquitous mass fragments typically generated from compounds carrying a trimethylsilyl moiety.}
  \item{compression.coef}{Data will be compressed when using OSD algorithm according to this value}
  \item{analysis.time}{The chromatographic retention time window to process. If 0, all the chromatogram is processed.}
 
 }
\details{
See eRah vignette for more details. To open the vignette, execute the following code in R:
vignette("eRahManual", package="erah")

(Available also in http://metabolomicsplatform.com/applications).
} 
 
\references{
[1] eRah: an R package for spectral deconvolution, alignment, and metabolite identification in GC/MS-based untargeted metabolomics. Xavier Domingo-Almenara, Alexandre Perera, Maria Vinaixa, Sara Samino, Xavier Correig, Jesus Brezmes, Oscar Yanes. (2016) Article in Press.
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{newExp}}, \code{\link{deconvolveComp}}, \code{\link{alignComp}}, \code{\link{setAlPar}}
}
\examples{

# The following will set eRah for analyzing the chromatograms
# from minutes 5 to 15, and withouth taking into account the masses
# 35:69,73:75,147:149, widht a minimum peak widht of 0.7 seconds.

# ex.dec.par <- setDecPar(min.peak.width=0.7, min.peak.height=5000, 
# noise.threshold=500, avoid.processing.mz=c(35:69,73:75,147:149), 
# analysis.time=c(5,15))

}


