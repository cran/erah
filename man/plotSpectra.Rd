% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotSpectra}
\alias{plotSpectra}
\alias{plotSpectra,MetaboSet-method}
\title{Plotting Spectra}
\usage{
plotSpectra(Experiment, AlignId, n.putative = 1,
compare = T, id.database = mslib, comp.db = NULL, 
return.spectra = F, draw.color = "purple", xlim = NULL)

\S4method{plotSpectra}{MetaboSet}(
  Experiment,
  AlignId,
  n.putative = 1,
  compare = T,
  id.database = mslib,
  comp.db = NULL,
  return.spectra = F,
  draw.color = "purple",
  xlim = NULL
)
}
\arguments{
\item{Experiment}{A 'MetaboSet' S4 object containing the experiment after being deconolved, aligned and (optionally) identified.}

\item{AlignId}{the Id identificator for the compound to be shown.}

\item{n.putative}{The hit number (position) to be returned when comparing the empirical spectrum with the reference. See details}

\item{compare}{logical. If TRUE, then the reference spectrum from the library is shown for comparison.}

\item{id.database}{The mass-spectra library to be compared with the empirical spectra. By default, the MassBank-[2] - Mass Bank of North America (MoNa) database are employed.}

\item{comp.db}{If you want to compare the empirical spectrum with another spectrum from the database, select the comp.db number from the database.}

\item{return.spectra}{logical. If TRUE, the function returns the empirical spectrum for the selected compound}

\item{draw.color}{Selects the color for the reference spectrum (see \code{\link{colors}}).}

\item{xlim}{x axsis (mass - m/z) limits (see \code{\link{plot.default}}).}
}
\value{
\code{plotSpectra} returns an vector when return.spectra=TRUE.

     \item{x}{vector. Containts the empirical spectrum.}
}
\description{
Plots the empirical spectra found by eRah, and allows comparing it with the reference spectra.
}
\details{
When identification is applied (see \code{\link{identifyComp}}), the number of hits to be returned (n.putative) has to be selected. Therefore, here you can compare the empirical spectrum (found by eRah) with each n.putative hit returned (1, 2, ...) by (see \code{\link{identifyComp}}).
}
\references{
[1] eRah: an R package for spectral deconvolution, alignment, and metabolite identification in GC/MS-based untargeted metabolomics. Xavier Domingo-Almenara, Alexandre Perera, Maria Vinaixa, Sara Samino, Xavier Correig, Jesus Brezmes, Oscar Yanes. (2016) Article in Press.

[2] MassBank: A public repository for sharing mass spectral data for life sciences, H. Horai, M. Arita, S. Kanaya, Y. Nihei, T. Ikeda, K. Suwa. Y. Ojima, K. Tanaka, S. Tanaka, K. Aoshima, Y. Oda, Y. Kakazu, M. Kusano, T. Tohge, F. Matsuda, Y. Sawada, M. Yokota Hirai, H. Nakanishi, K. Ikeda, N. Akimoto, T. Maoka, H. Takahashi, T. Ara, N. Sakurai, H. Suzuki, D. Shibata, S. Neumann, T. Iida, K. Tanaka, K. Funatsu, F. Matsuura, T. Soga, R. Taguchi, K. Saito and T. Nishioka, J. Mass Spectrom., 45, 703-714 (2010)
}
\seealso{
\code{\link{plotProfile}} \code{\link{plotAlign}}
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}
