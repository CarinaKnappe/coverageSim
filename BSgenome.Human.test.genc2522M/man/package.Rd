\name{BSgenome.Human.test.genc2522M}
\docType{package}

\alias{BSgenome.Human.test.genc2522M-package}
\alias{BSgenome.Human.test.genc2522M}
\alias{Human.test}

\title{Full genome sequences for Human.test, genc2522M}

\description{
  Full genome sequences for Human.test, genc2522M
}

\details{
  
}

\note{
  This BSgenome data package was made from the following source data files:
  \preformatted{
-- information not available --
  }
}

\author{The Bioconductor Dev Team}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects in the \pkg{BSgenome}
          software package.

    \item The \code{\link[Seqinfo]{seqinfo}} getter and \link[Seqinfo]{Seqinfo}
          objects in the \pkg{Seqinfo} package.

    \item The \code{\link[GenomeInfoDb]{seqlevelsStyle}} getter and
          setter in the \pkg{GenomeInfoDb} package.

    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.

    \item The "A quick introduction to the BSgenomeForge package" vignette
          in the \pkg{BSgenomeForge} package for how to create a BSgenome
          data package. Install the package and access the vignette with
          \code{browseVignettes(package="BSgenomeForge")}.
  }
}

\examples{
BSgenome.Human.test.genc2522M
bsg <- BSgenome.Human.test.genc2522M
head(seqlengths(bsg))
seqinfo(bsg)



## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------

## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
