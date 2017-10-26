#' Enumerate genotypic vectors
#'
#' Enumerates all genotypic vectors of length k consistent with each identity vector.
#'
#' @param      k (numeric integer) specifies the number of putatitive replicates associated with 
#'               the proband in question.
#'        
#' @param I.list (list) output from \code{\link{enumerateI}}.        
#'
#' @return The function returns a list of length two: 
#'      (1) `G|I`: (list) enumerates all genotype vectors consistent with identity vector i
#'      (2) `P(G|I)`: (vector) provides the probability of genotype vector g given identity vector i. 
#'                    We compute P(G|I) assuming the uniform probability law (i.e. 1/cardinality of set).
#'
#' @author Ariel W Chan, \email{ac2278@@cornell.edu}
#' 
#' @export
################################################################################
enumerateG_revised = function(k, I.list){
  require(gtools); require(parallel)
  # (1) Supply a list of all possible identity vectors of length k (i.e. I.list).
  #     We index each identity vector using i. 
  I <- do.call(rbind, unlist(I.list, recursive=F))
  rownames(I) <- paste("I=", apply(I, 1, function(i) paste(i, collapse=",")), sep="")
  # (2) Enumerate all possible genotype vectors of length k. (indexed by g)
  G <- permutations(n=3, r=k, v=c("AA","AB","BB"), repeats.allowed=T)
  rownames(G) <- apply(G, 1, function(g) paste(g, collapse=","))
  colnames(G) <- paste("d=", 1:k, sep="")
  # (3) We answer the question "Does g belong to category i?". 
  membership <- apply(G, 1, function(g){ x <- length(unique(g));
                                         putative <- do.call(rbind, unlist(I.list[x:k], recursive=F))
                                         rownames(putative) <- apply(putative, 1, function(i) paste(i, collapse=","))
                                         consistent = apply(putative, 1, function(evaluate){
                                           count = 0 
                                           for (i in 1:k){
                                             test = length(unique(g[evaluate==i]))>1
                                             if(test==T){ break } else count = count + 1
                                           }
                                           return(count)
                                         }) 
                                         return(paste("I=",names(consistent[consistent==k]),sep="")) })
  
  enumerateGgivenI <- sapply(setNames(rownames(I), rownames(I)), function(i) names(grep(T, rapply(membership, function(x) ifelse(length(intersect(x,i))==1, T, F)), value=T)))
  # We return a list of length two: 
  #      (1) `G|I`: (list) enumerates all genotype vectors consistent with identity vector i
  #      (2) `P(G|I)`: (vector) provides the probability of genotype vector g given identity vector i. 
  #                    We compute P(G|I) assuming (for now) the uniform probability law (i.e. 1/cardinality of set).                  
  return(list("G"=G, "G|I"=enumerateGgivenI, "P(G|I)"=rapply(enumerateGgivenI, function(x) setNames(rep(1/length(x), length(x)), x), how="list")))
}










