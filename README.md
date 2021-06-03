[![pipeline status](https://gitlab.gwdg.de/mpievolbio-it/diststring/badges/master/pipeline.svg)](https://gitlab.gwdg.de/mpievolbio-it/diststring/-/commits/master)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)
[![MITlicense](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)


distSTRING - calculates pairwise distances between all sequences of a DNAStringSet or a AAStringSet
=========

# Installation

see also here for the R package pages [https://mpievolbio-it.pages.gwdg.de/diststring/](https://mpievolbio-it.pages.gwdg.de/diststring/)

## R specific installation prerequisites

### install packages from [cran](https://cran.r-project.org/web/packages/index.html)

In most cases you need to first install the following system-wide packages to be able to compile the R dependencies.

Ubuntu/Debian

```
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev libglu1-mesa-dev libgit2-dev
#pkgdown dependencies - pkgdown is used to build R package pages
#sudo apt-get install libssh2-1-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev
```

CentOS

```
sudo yum install libcurl-devel openssl-devel libxml2-devel mesa-libGLU-devel libgit2-devel
#pkgdown dependencies - pkgdown is used to build R package pages
#sudo yum install libssh2-devel fontconfig-devel harfbuzz-devel fribidi-devel
```

- [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
- [RcppThread](https://cran.r-project.org/web/packages/RcppThread/index.html)
- [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
- [testthat](https://cran.r-project.org/web/packages/testthat/index.html)
- [curl](https://cran.r-project.org/web/packages/curl/index.html)
- [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html)
- [ape](https://cran.r-project.org/web/packages/ape/index.html)

```
install.packages("Rcpp")
install.packages("RcppThread")
install.packages("devtools")
install.packages("testthat")
install.packages("curl")
install.packages("seqinr")
install.packages("ape")
```

### install packages from [Bioconductor](https://www.bioconductor.org/)

- [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```

### install [distSTRING](https://gitlab.gwdg.de/mpievolbio-it/diststring)

```
library(devtools)
install_gitlab("mpievolbio-it/diststring", host = "https://gitlab.gwdg.de",
build_vignettes = FALSE, dependencies = FALSE)
#install_github("kullrich/distSTRING", build_vignettes = FALSE, dependencies = FALSE)
```

## Quick-guide

```
library(distSTRING)
## load example sequence data
data("hiv", package="distSTRING")

## calculate pairwise AA distances based on Grantham's distance
aa.dist <- aastring2dist(cds2aa(hiv), score=granthamMatrix())
head(aa.dist$distSTRING)

## create and plot bionj tree
aa.dist.bionj <- ape::bionj(as.dist(aa.dist$distSTRING))
plot(aa.dist.bionj)

## calculate pairwise DNA distances based on IUPAC distance
dna.dist <- dnastring2dist(hiv, score=iupacMatrix())

## create and plot bionj tree
dna.dist.bionj <- ape::bionj(as.dist(dna.dist$distSTRING))
head(dna.dist$distSTRING)

## creation of the association matrix:
association <- cbind(aa.dist.bionj$tip.label, aa.dist.bionj$tip.label)

## cophyloplot
ape::cophyloplot(aa.dist.bionj,
                 dna.dist.bionj,
                 assoc=association,
                 length.line=4,
                 space=28,
                 gap=3,
                 rotate=TRUE)

## calculate pairwise DNA distances based on K80 distance
dna.dist.K80 <- dnastring2dist(hiv, model="K80")

## calculate pairwise AA distances based on getAAMatrix() function from the alakazam package
data("AAMatrix", package="distSTRING")
aa.dist <- aastring2dist(cds2aa(hiv), score=AAMatrix)

## example how to calculate all pairwise kaks values given a MSA
hiv_kaks <- dnastring2kaks(hiv, model="Li")
hiv_kaks <- dnastring2kaks(hiv, model="NG86")

## codon plot - sites under possible positive selection
library(tidyr)
library(dplyr)
library(ggplot2)
hiv.xy <- codonmat2xy(dnastring2codonmat(hiv))
hiv.xy %>% select(Codon,SynMean,NonSynMean,IndelMean) %>%
  gather(variable, values, -Codon) %>% 
  ggplot(aes(x=Codon, y=values)) + 
    geom_line(aes(colour=factor(variable))) + 
    geom_point(aes(colour=factor(variable))) + 
    ggtitle("HIV-1 sample 136 patient 1 from Sweden envelope glycoprotein (env) gene")
```

## License

MIT (see LICENSE)

## Contributing Code

If you would like to contribute to distSTRING, please file an issue so that one can establish a statement of need, avoid redundant work, and track progress on your contribution.

Before you do a pull request, you should always file an issue and make sure that someone from the distSTRING developer team agrees that it’s a problem, and is happy with your basic proposal for fixing it.

Once an issue has been filed and we've identified how to best orient your contribution with package development as a whole, [fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) the [main repo](https://gitlab.gwdg.de/mpievolbio-it/diststring.git), branch off a [feature branch](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-branches) from `master`, [commit](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/committing-and-reviewing-changes-to-your-project) and [push](https://docs.github.com/en/github/using-git/pushing-commits-to-a-remote-repository) your changes to your fork and submit a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests) for `distSTRING:master`.

By contributing to this project, you agree to abide by the Code of Conduct terms.

## Bug reports

Please report any errors or requests regarding [distSTRING](https://gitlab.gwdg.de/mpievolbio-it/diststring) to Kristian Ullrich (ullrich@evolbio.mpg.de)

or use the issue tracker at [https://gitlab.gwdg.de/mpievolbio-it/diststring/issues](https://gitlab.gwdg.de/mpievolbio-it/diststring/issues)

## Code of Conduct - Participation guidelines

This repository adhere to [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://gitlab.gwdg.de/mpievolbio-it/diststring/-/blob/master/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.
