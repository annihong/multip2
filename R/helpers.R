#' Get pairs given total number of actors/layers
#' 
#' @param n numeric: number of actors (or number of layers)
#'
#' @return list of vectors of size 2: indices of the pairs of actors (layers) in the networks (dyads) 
get_dyads <- function(n){
  res = list()
  counter = 1
  for (i in 1:n) {
    for (j in i:n) {
      if (i != j) {
        res[[counter]] = c(i,j)
        counter = counter + 1
      }
    }
  }
  if (n < 2) {
    res = list(c(0,0))
  }
  return(res)
}

#' Obtain dyad outcomes from network matrices
#' 
#' @param Ms list of adj matrices representing the multiplex network
#'
#' @return return a vector y of size N, such that y in {1,...,2^(2*t)}^N representing the outcomes on a pair of dyads.

matrix_to_dyads_nd <- function(Ms) {
  #given a vector of dyad-wise outcome for dyads i
  #return the multiplex outcome y \in {1,...,2^(2*t)}
  helper <- function(ns){
    res <- sapply(2:t, function(i) (ns[i] - 1)*4^(i - 1))
    res = sum(res) + ns[1]
    return(res)
  }
  t = length(Ms)
  y_1d <- do.call(cbind, lapply(Ms, matrix_to_dyads_1d))
  y = apply(y_1d, 1, helper)
  return(y)
}


#' Obtain dyad outcomes from a single network layer
#' 
#' @param M a single adj matrices representing one layer of the multiplex network
#'
#' @return return a vector y of size N the dyad-wise outcome y in {1,2,3,4}^N

matrix_to_dyads_1d <- function(M){
  n = nrow(M)
  y <- c()
  for (i in 1:n){
    for (j in i:n){
      i_to_j = M[i,j]
      j_to_i = M[j,i]
      if (i != j) {
        outcome = 1
        if (M[i,j] & ! M[j,i]) {
          outcome = 2
        } else if (!M[i,j] & M[j,i]) {
          outcome = 3
        } else if (M[i,j] & M[j,i]) {
          outcome =4
        }
        y <- c(y, outcome)
      }
    }
  }
  return(y)
}


