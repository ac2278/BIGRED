#' Normalized genotype likelihoods
#'
#' Computes a (normalized) likelihood for each of the three possible genotypes (i.e. AA, AB, and BB) 
#' of a diploid individual at a given biallelic site.
#'
#' @param             x (character vector) specifies the observed counts of the reference allele (allele A) and 
#'                      alternative allele (allele B) at a site for a diploid individual
#'
#' @param         error (numeric float) specifies a fix sequencing error rate (errors that may have been introduced 
#'                      in base calling, alignment and assembly. This value may range from 0 to 1. 
#'                      For an error rate of 1 percent, enter 0.01. 
#' 
#' @param callthreshold (numeric float; optional) only call genotypes with a posterior probability (calculated 
#'                      assuming a uniform prior) of >callthreshold. If no posterior probability is greater than the 
#'                      callthreshold, the function returns NA. If no argument is specified, the function returns 
#'                      the normalized likelihoods (or equivalently, the posterior probabilities calculated assuming 
#'                      a uniform prior on the three genotypes) for all three genotypes.
#' 
#' @param        coding (character) specifies the format of the function output. The user may select either
#'                      'character' or 'numeric'.
#'
#' @examples 
#' normGL(x="1,5", error=0.01, coding="character")         
#'           AA           AB           BB 
#' 3.938746e-09 6.216456e-01 3.783544e-01
#'
#' normGL(x='1,5', error=0.01, coding="numeric")
#'            0            1            2 
#' 3.938746e-09 6.216456e-01 3.783544e-01
#'                          
#' @return Refer to paramemter callthreshold and coding.
#'
#' @author Ariel W Chan, \email{ac2278@@cornell.edu}
#' 
#' @export
################################################################################
normGL = function(x, error, callthreshold, coding){ #, where x is a character vector of length one "REFcount, ALTcount"
  REFcount <- as.numeric(unlist(strsplit(x,","))[1])
  ALTcount <- as.numeric(unlist(strsplit(x,","))[2])
  n <- REFcount + ALTcount
  if(n == 0) return(setNames(rep(1/3, 3), c("AA","AB","BB"))) else
  k <- ALTcount
  e <- error
  likelihood_AA <- choose(n, k) * e^k * (1-e)^(n-k)
  likelihood_AB <- choose(n, k) * (0.50)^k * (0.50)^(n-k)
  likelihood_BB <- choose(n, n-k) * (1-e)^k * e^(n-k)
  normalization_constant <- sum(likelihood_AA, likelihood_AB, likelihood_BB)
  
  GL_AA <- likelihood_AA/normalization_constant # Normalize the likelihoods such that they sum to one.
  GL_AB <- likelihood_AB/normalization_constant
  GL_BB <- likelihood_BB/normalization_constant
  
  if(coding=="character"){GL <- c("AA"=GL_AA, "AB"=GL_AB, "BB"=GL_BB)}
  if(coding=="numeric"){GL <- c("0"=GL_AA, "1"=GL_AB, "2"=GL_BB)}
  
  if(missing(callthreshold)){return(GL)} else
  return(ifelse(max(GL)>callthreshold, names(which.max(GL)), NA))

}
