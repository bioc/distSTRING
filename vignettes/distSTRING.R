## ----echo=FALSE, results='hide', warning=FALSE, message=FALSE-----------------
suppressPackageStartupMessages({
    library(distSTRING)
    library(Biostrings)
    library(GenomeInfoDb)
    library(tidyverse)
    library(grDevices)
    })

## -----------------------------------------------------------------------------
# load distSTRING
library(distSTRING)

## -----------------------------------------------------------------------------
## define two cds sequences
cds1 <- Biostrings::DNAString("ATGCAACATTGC")
cds2 <- Biostrings::DNAString("ATG---CATTGC")
cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
    Biostrings::DNAStringSet(cds2))
## define names
names(cds1.cds2.aln) <- c("seq1", "seq2")
## convert into alignment
cds1.cds2.aln |> dnastring2aln()

## -----------------------------------------------------------------------------
## convert back into DNAStringSet
cds1.cds2.aln |> dnastring2aln() |> aln2dnastring()

## -----------------------------------------------------------------------------
## convert into alignment
cds1.cds2.aln |> dnastring2dnabin()

## -----------------------------------------------------------------------------
## convert back into DNAStringSet
cds1.cds2.aln |> dnastring2dnabin() |> dnabin2dnastring()
## use woodmouse data
data(woodmouse, package="ape")
woodmouse |> dnabin2dnastring()

## -----------------------------------------------------------------------------
## translate cds into aa
aa1.aa2.aln <- cds1.cds2.aln |> cds2aa()
## convert into alignment
aa1.aa2.aln |> aastring2aln()

## -----------------------------------------------------------------------------
## convert back into AAStringSet
aa1.aa2.aln |> aastring2aln() |> aln2aastring()

## -----------------------------------------------------------------------------
## convert into AAbin
aa1.aa2.aln |> aastring2aabin()

## -----------------------------------------------------------------------------
## convert back into AAStringSet
aa1.aa2.aln |> aastring2aabin() |> aabin2aastring()

## -----------------------------------------------------------------------------
## define two cds sequences
cds1 <- Biostrings::DNAString("ATGCAACATTGC")
cds2 <- Biostrings::DNAString("ATG---CATTGC")
cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
    Biostrings::DNAStringSet(cds2))
## define names
names(cds1.cds2.aln) <- c("seq1", "seq2")
## translate cds into aa
cds1.cds2.aln |> cds2aa()
aa1.aa2.aln <- cds1.cds2.aln |> cds2aa()

## -----------------------------------------------------------------------------
## translate cds into aa using frame = 2
## result is empty due to not multiple of three
cds1.cds2.aln |> cds2aa(frame=2)
## translate cds into aa using frame = 2 and shorten = TRUE
cds1.cds2.aln |> cds2aa(frame=2, shorten=TRUE)
## translate cds into aa using frame = 3 and shorten = TRUE
cds1.cds2.aln |> cds2aa(frame=3, shorten=TRUE)
## use woodmouse data
data(woodmouse, package="ape")
woodmouse |> dnabin2dnastring() |> cds2aa(shorten=TRUE)

## -----------------------------------------------------------------------------
## alternative genetic code
data(woodmouse, package="ape")
woodmouse |> dnabin2dnastring() |> cds2aa(shorten=TRUE,
                                          genetic.code=Biostrings::getGeneticCode("2"))

## -----------------------------------------------------------------------------
## load example sequence data
data("hiv", package="distSTRING")

## calculate pairwise AA distances based on Grantham's distance
aa.dist <- hiv |> cds2aa() |> aastring2dist(score=granthamMatrix())
## obtain distances
head(aa.dist$distSTRING)
## obtain pairwise sites used
head(aa.dist$sitesUsed)

## ----sessionInfo, echo=TRUE---------------------------------------------------
sessionInfo()

