% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnastring2dnabin.R
\name{dnastring2dnabin}
\alias{dnastring2dnabin}
\title{dnastring2dnabin}
\usage{
dnastring2dnabin(dna)
}
\arguments{
\item{dna}{\code{DNAStringSet} [mandatory]}
}
\value{
An object of class \code{DNAbin}
}
\description{
This function converts a \code{DNAStringSet} into an
\code{ape} \code{DNAbin}.
}
\examples{
## define two cds sequences
cds1 <- Biostrings::DNAString("ATGCAACATTGC")
cds2 <- Biostrings::DNAString("ATG---CATTGC")
cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
    Biostrings::DNAStringSet(cds2))
## convert into DNAbin
#dnastring2dnabin(cds1.cds2.aln)
cds1.cds2.aln |> dnastring2dnabin()
}
\seealso{
\code{\link[seqinr]{as.alignment}}
\code{\link[ape]{as.DNAbin.alignment}}
}
\author{
Kristian K Ullrich
}
