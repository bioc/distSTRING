#' @title addmask2string
#' @name addmask2string
#' @description This function adds mask information as an \code{IRanges} object,
#' \code{START} and \code{END} information, to a
#' \code{DNAStringSet} or an \code{AAStringSet} and puts them into the
#' \code{metadata} information.
#' This information can be used to restrict the distance calculation to
#' specific regions of the \code{DNAStringSet} or the \code{AAStringSet}.
#' @param seq \code{DNAStringSet} or \code{AAStringSet} [mandatory]
#' @param mask \code{IRanges} object [mandatory]
#' @param append indicate if mask should be appended or overwritten
#' [default: TRUE]
#' @return An object of class \code{DNAStringSet} or \code{AAStringSet}
#' @importFrom methods is
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet
#' readDNAStringSet readAAStringSet writeXStringSet width subseq
#' @importFrom IRanges IRanges IRangesList reduce start end findOverlaps
#' disjoin overlapsRanges
#' @seealso \code{\link[distSTRING]{addregion2string}},
#' \code{\link[distSTRING]{addpop2string}},
#' \code{\link[distSTRING]{addpos2string}}
#' @examples
#' data(iupac, package="distSTRING")
#' iupac.aa <- iupac |> cds2aa(shorten = TRUE)
#' ## create mask
#' mask1 <- IRanges::IRanges(start=c(1,41), end=c(20,50))
#' ## add mask
#' iupac.aa <- iupac.aa |> addmask2string(mask=mask1)
#' iupac.aa@metadata$mask
#' ## append mask
#' mask2 <- IRanges::IRanges(start=c(21), end=c(30))
#' iupac.aa <- iupac.aa |> addmask2string(mask=mask2)
#' iupac.aa@metadata$mask
#' ## overwrite mask
#' iupac.aa <- iupac.aa |> addmask2string(mask=mask2, append=FALSE)
#' iupac.aa@metadata$mask
#' ## reduce by mask
#' iupac.aa.region <- iupac.aa |> string2region(mask=iupac.aa@metadata$mask)
#' iupac.aa.region@metadata
#' @export addmask2string
#' @author Kristian K Ullrich

addmask2string <- function(seq, mask = NULL, append = TRUE){
    stopifnot("Error: input needs to be a DNAStringSet or AAStringSet"=
        methods::is(seq, "AAStringSet") | methods::is(seq, "DNAStringSet"))
    stopifnot("Error: set mask"= !is.null(mask))
    stopifnot("Error: mask needs to be an IRanges object"=
        methods::is(mask, "IRanges"))
    if(append){
        if(is.null(seq@metadata$mask)){
            seq@metadata$mask <- IRanges::reduce(mask)
        } else{
            seq@metadata$mask <- IRanges::reduce(
                c(seq@metadata$mask, mask))
        }
    } else{
        seq@metadata$mask <- IRanges::reduce(mask)
    }
    return(seq)
}
