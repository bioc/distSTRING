#' @title aastring2aabin
#' @name aastring2aabin
#' @description This function converts a \code{AAStringSet} into an \code{ape}
#' \code{DNAbin}.
#' @param aa \code{AAStringSet} [mandatory]
#' @return An object of class \code{DNAbin}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom ape as.AAbin
#' @seealso \code{\link[seqinr]{as.alignment}}
#' \code{\link[ape]{as.DNAbin.alignment}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'  Biostrings::DNAStringSet(cds2))
#' ## convert into alignment
#' aastring2aabin(cds2aa(cds1.cds2.aln))
#' @export aastring2aabin
#' @author Kristian K Ullrich

aastring2aabin <- function(aa){
    if(class(aa)!="AAStringSet"){stop("Error: input needs to be a AAStringSet")}
    return(ape::as.AAbin(aa))
}
