#' @title dnastring2dist
#' @name dnastring2dist
#' @description This function calculates pairwise distances for all
#' combinations of a \code{DNAStringSet}.
#' @param dna \code{DNAStringSet} [mandatory]
#' @param model specify model either "IUPAC" or any model from
#' \code{ape::dist.dna} [default: IUPAC]
#' @param threads number of parallel threads [default: 1]
#' @param score \code{scoringMatrix} use scoring matrix to calculate
#' distances [default: NULL]
#' @param ... other \code{ape::dist.dna} parameters
#' (see \code{\link[ape]{dist.dna}})
#' @return A data.frame of pairwise distance values \code{distSTRING} and
#' sites used \code{sitesUsed}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom ape dist.dna
#' @seealso \code{\link[ape]{dist.dna}}
#' @examples
#' ## load example sequence data
#' data("hiv", package="distSTRING")
#' dnastring2dist(hiv, model="IUPAC")
#' dnastring2dist(hiv, model="K80")
#' data("woodmouse", package="ape")
#' dnastring2dist(dnabin2dnastring(woodmouse), score=iupacMatrix())
#' dnastring2dist(hiv, model = "IUPAC", threads = 2)
#' @export dnastring2dist
#' @author Kristian K Ullrich

dnastring2dist <- function(dna,
    model = "IUPAC",
    threads = 1,
    score = NULL,
    ...){
    if(class(dna) != "DNAStringSet"){
        stop("Error: Input needs to be DNAStringSet")
    }
    if(is.null(score) & !model %in%
        c("IUPAC", "raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81",
        "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin",
        "indel", "indelblock")){
            stop("Error: either choose model 'IUPAC' or '?ape::dist.dna'")
        }
    if(!is.null(score)){
        OUT <- rcpp_distSTRING(dnavector = as.character(dna),
            scoreMatrix = score, ncores = threads)
        OUT$distSTRING <- as.data.frame(OUT$distSTRING)
        OUT$sitesUsed <- as.data.frame(OUT$sitesUsed)
    }
    if(is.null(score) & model == "IUPAC"){
        OUT <- rcpp_distSTRING(dnavector = as.character(dna),
            scoreMatrix=iupacMatrix(), ncores=threads)
        OUT$distSTRING <- as.data.frame(OUT$distSTRING)
        OUT$sitesUsed <- as.data.frame(OUT$sitesUsed)
    }
    if(is.null(score) & model != "IUPAC"){
        distSTRING_ <- as.matrix(ape::dist.dna(x=dnastring2dnabin(dna),
            model=model, ...))
        sitesUsed_ <- rcpp_pairwiseDeletionDNA(dnavector=as.character(dna),
            ncores=threads)
        OUT <- list(distSTRING = as.data.frame(distSTRING_))
        OUT <- append(OUT, sitesUsed_)
        OUT$sitesUsed <- as.data.frame(OUT$sitesUsed)
    }
    return(OUT)
}
