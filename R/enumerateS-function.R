#' Enumerate source vectors
#'
#' Enumerates all source vectors associated with the k putative replicates of a given proband
#'
#' @param k (numeric integer) specifies the number of putative replicates associated with 
#'        the proband in question.   
#'
#' @return The function returns a named list enumerating all of the source vectors that describe how 
#'         the k putative replicates could potentially be related to one another.
#'
#' @examples 
#' k3 <- enumerateS(k=3) 
#' This function call returns a named list consisting of three elements, which are also themselves 
#' lists: 'i=1', 'i=2', and 'i=3'.
#' 
#' Printing the output of enumerateS(k=3)[['i=2']] allows us to see the source vectors that describe the case 
#' where the k=3 supposed replicates actually originated from i=2 genotypic sources:
#' $`120`
#' [,1] [,2] [,3]
#' [1,]    1    2    2     # corresponds to the source vector S=1,2,2, where samples d=2 and d=3 originate from one genotype source and d=1 is assumed to be the outlier
#' $`210`
#' [,1] [,2] [,3]
#' [1,]    1    1    2     # corresponds to the source vector S=1,1,2, where samples d=1 and d=2 originate from one genotype source and d=3 is assumed to be the outlier
#' [2,]    1    2    1     # corresponds to the source vector S=1,2,1, where samples d=1 and d=3 originate from one genotype source and d=2 is assumed to be the outlier
#'
#' @author Ariel W Chan, \email{ac2278@@cornell.edu}
#' 
#' @export
################################################################################
enumerateS = function(k){
  require(stringr); require(gtools)
  A <- 1:(k-1)
  composition <- k
  input <- setNames(1:k,1:k)
  count <- 1
  while(count<k){
    output = unlist(sapply(input, function(x) structure(x+A,names=A), simplify=F))              #  
    composition = c(composition, names(output[output==k]))                                      # Output outcome if the elements sum to k
    input = output[output<k]                                                                    # Eliminate outcome if the elements sum to >k
    count = count + 1
  }
  composition = gsub("[.]", "", composition)                                                    # Replace all "." with "". 
  composition = str_pad(composition, k, side="right", pad="0")                                  # 
  order = sapply(composition, function(x) matrix(rep(1:k, times=unlist(strsplit(x,""))), ncol=k), simplify=F)
  for (i in 2:(length(composition)-1)){                                                         #  
    composition.permutation = cbind(1, unique(permutations(k-1, k-1, order[[i]][,2:k], set=F)))
    
    # Identify the subset of permutations that satisfy constraint (1): 
    order[[i]] = matrix(composition.permutation[apply(composition.permutation, 1, function(x){
      length(unlist(sapply(2:k, function(i) setdiff(1:(x[i]), x[1:i]))))==0 }),], ncol=k)
  }
  I = rapply(order, function(x) length(table(x))) 
  S = sapply(setNames(1:k,paste("i=",1:k,sep="")), simplify=F, function(i) order[I==i])
  return(S)
}

