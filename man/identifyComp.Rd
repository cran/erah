\name{identifyComp}
\alias{identifyComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identification of compounds}
\description{Identification of compounds. Each empirical spectrum is compared against a ms library.}
\usage{
identifyComp(Experiment, id.database=mslib,
mz.range=NULL, n.putative=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data previously created by newExp, deconvolved by deconvolveComp and optionally aligned by alignComp.}
  \item{id.database}{The mass-spectra library to be compared with the empirical spectra. By default, the MassBank-[2] - Mass Bank of North America (MoNa) database are employed.}
  \item{mz.range}{The same as in alignComp. If specified already in alignComp, then there is no need to especify it again. If not, it has to be specified.}
  \item{n.putative}{The number of hits (compound candidate names) to be returned for each spectrum found.}
}

\value{
The function returns an updated S4 'MetaboSet' class, where the GC-MS samples have been now aligned.
}

\references{
[1] eRah: an R package for spectral deconvolution, alignment, and metabolite identification in GC/MS-based untargeted metabolomics. Xavier Domingo-Almenara, Alexandre Perera, Maria Vinaixa, Sara Samino, Xavier Correig, Jesus Brezmes, Oscar Yanes. (2016) Article in Press.

[2] MassBank: A public repository for sharing mass spectral data for life sciences, H. Horai, M. Arita, S. Kanaya, Y. Nihei, T. Ikeda, K. Suwa. Y. Ojima, K. Tanaka, S. Tanaka, K. Aoshima, Y. Oda, Y. Kakazu, M. Kusano, T. Tohge, F. Matsuda, Y. Sawada, M. Yokota Hirai, H. Nakanishi, K. Ikeda, N. Akimoto, T. Maoka, H. Takahashi, T. Ara, N. Sakurai, H. Suzuki, D. Shibata, S. Neumann, T. Iida, K. Tanaka, K. Funatsu, F. Matsuura, T. Soga, R. Taguchi, K. Saito and T. Nishioka, J. Mass Spectrom., 45, 703-714 (2010)
 
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{newExp}}, \code{\link{alignComp}}, \code{\link{setAlPar}}, \code{\link{setDecPar}}
}



