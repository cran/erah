\name{plotAlign}
\alias{plotAlign}
\title{Plotting chromatophic profile with and without alignement}
\description{Plots the chromatophic profiles of the compounds found by eRah. Similarly to plotProfile, but with two sub-windows, showing the chromatophic profiles before and after alignment.}
\usage{
plotAlign(Experiment,AlignId, per.class=T, xlim=NULL)
}
\arguments{
 \item{Experiment}{A 'MetaboSet' S4 object containing the experiment after being deconolved, aligned and (optionally) identified.}
   \item{AlignId}{the Id identificator for the compound to be shown.}
   \item{per.class}{logical. if TRUE the profiles are shown one color per class, if FALSE one color per sample.}
   \item{xlim}{x axsis (retention time) limits (see \code{\link{plot.default}}).}
}
\author{
Xavier Domingo-Almenara. xavier.domingoa@eurecat.org
}

\seealso{
\code{\link{plotSpectra}}, \code{\link{plotProfile}}
}
