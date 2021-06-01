#' @title cds2aa
#' @name cds2aa
#' @description This function translates a \code{DNAStringSet} into an \code{AAStringSet}.
#' @param cds \code{DNAStringSet} [mandatory]
#' @param shorten shorten all sequences to multiple of three [default: FALSE]
#' @param frame  indicates the first base of a the first codon [default: 1]
#' @param framelist  supply vector of frames for each entry [default: NULL]
#' @return \code{AAStringSet}
#' @importFrom Biostrings DNAString DNAStringSet AAString AAStringSet readDNAStringSet readAAStringSet writeXStringSet width subseq translate
#' @importFrom stringr word
#' @seealso \code{\link[Biostrings]{XStringSet-class}},
#' \code{\link[seqinr]{translate}}
#' @examples
#' ## define two cds sequences
#' cds1 <- Biostrings::DNAString("ATGCAACATTGC")
#' cds2 <- Biostrings::DNAString("ATG---CATTGC")
#' cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
#'  Biostrings::DNAStringSet(cds2))
#' cds2aa(cds1.cds2.aln)
#' @export cds2aa
#' @author Kristian K Ullrich

cds2aa <- function(cds, shorten = FALSE, frame = 1, framelist = NULL){
  if(class(cds)!="DNAStringSet"){stop("Error: input needs to be a DNAStringSet")}
  if(!frame %in% c(1, 2, 3)){stop("Error: frame needs to be 1 or 2 or 3")}
  if(!is.null(framelist)){
    if(length(framelist) != length(cds)){stop("Error: framelist needs to be of equal length as cds")}
  }
  if(!is.null(names(cds))){
    names(cds) <- stringr::word(names(cds), 1)
  }
  if(is.null(framelist)){
    cds <- Biostrings::subseq(cds, frame, unique(Biostrings::width(cds)))
  }
  if(!is.null(framelist)){
    cds <- Biostrings::subseq(cds, framelist, unique(Biostrings::width(cds)))
  }
  if(shorten){
    cds <- Biostrings::subseq(cds, 1, Biostrings::width(cds) - Biostrings::width(cds) %% 3)
  }
  cds_not_multiple_of_three.idx <- which(Biostrings::width(cds) %% 3 != 0)
  if(length(cds_not_multiple_of_three.idx) > 0){
    cds_not_multiple_of_three <- cds[cds_not_multiple_of_three.idx]
    cds <- cds[-cds_not_multiple_of_three.idx]
  }
  #aa <- Biostrings::AAStringSet(unlist(lapply(as.character(cds), function(x) paste0(seqinr::translate(unlist(strsplit(x, ""))), collapse=""))))
  #aa <- Biostrings::AAStringSet(unlist(lapply(as.character(cds), function(x) seqinr::c2s(seqinr::translate(seqinr::s2c(x))))))
  cds <- DNAStringSet(gsub("-", "N", cds))
  cds <- DNAStringSet(gsub("X", "N", cds))
  aa <- Biostrings::translate(cds, if.fuzzy.codon = "X")
  return(aa)
}
