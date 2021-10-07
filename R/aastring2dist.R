#' @title aastring2dist
#' @name aastring2dist
#' @description This function calculates pairwise distances for all combinations
#' of a \code{AAStringSet}.
#' @param aa \code{AAStringSet} [mandatory]
#' @param threads number of parallel threads [default: 1]
#' @param score \code{scoringMatrix} use scoring matrix to calculate distances
#' [default: NULL]
#' @return A data.frame of pairwise distance values \code{distSTRING} and sites
#' used \code{sitesUsed}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' pairwiseAlignment
#' @importFrom methods is
#' @seealso \code{\link[distSTRING]{dnastring2dist}}
#' @examples
#' ## load example sequence data
#' data("hiv", package="distSTRING")
#' #aastring2dist(cds2aa(hiv), score=granthamMatrix())
#' hiv |> cds2aa() |> aastring2dist(score=granthamMatrix())
#' @export aastring2dist
#' @author Kristian K Ullrich

aastring2dist <- function(aa,
    threads = 1,
    score = NULL){
    stopifnot("Error: input needs to be a AAStringSet"=
        methods::is(aa, "AAStringSet"))
    stopifnot("Error: set score matrix e.g 'granthamMatrix()'"= !is.null(score))
    OUT <- distSTRING::rcpp_distSTRING(dnavector=as.character(aa),
        scoreMatrix=score, ncores=threads)
    OUT$distSTRING <- as.data.frame(OUT$distSTRING)
    OUT$sitesUsed <- as.data.frame(OUT$sitesUsed)
    return(OUT)
}
