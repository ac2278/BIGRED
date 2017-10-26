#' Thin sites
#'
#' Thins a set of sites such that no two sites are within a specified distance (in bp) from one another.
#'
#' @param             x (numeric vector) specifies the set of sites the user wishes to thin.
#'
#' @param      distance (numeric) specifies the desired minimum distance (in bp) between any two sites in the 
#'                      resulting set
#'                          
#' @return The function returns a set of thinned sites where no two sites are within a specified distance 
#'         (in bp) of each other.
#'
#' @author Ariel W Chan, \email{ac2278@@cornell.edu}
#' 
#' @export
################################################################################
thin <- function(x, distance){
    if(length(grep("_", x))>0){
        pos <- do.call(rbind, strsplit(x, "_"));
        mode(pos) <- "numeric"
        chrom <- pos[1,1]; pos <- pos[,2]
    } else pos <- x
    
    sites <- start <- pos[1]
    while(length(pos)>1){
        pass <- ifelse(pos>=(start+distance), T, F)
        pos <- pos[pass]
        start <- pos[1]
        sites <- c(sites, start)
    }
    sites <- sites[!is.na(sites)]
    if(length(grep("_", x))>0){ sites <- paste(chrom, sites, sep="_") }
    return(sites)
}