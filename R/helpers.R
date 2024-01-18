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

get_pair_names <- function(outcome) {
    t = length(outcome)
    pairs <- get_dyads(t)
    lapply(pairs, function(x) paste0(outcome[x[1]],":", outcome[x[2]]))
}


#' Convert Dyadic Data to List of Adjacency Matrices
#'
#' This function takes a dataframe of dyadic data and converts it to a list of adjacency matrices.
#' Each row of the dataframe is a network in dyad form (1 x N) with B total such networks. N is the number of dyads in the network.
#' The function returns a list of B lists of t adjacency matrices representing the B multiplex networks.
#'
#' @param dyad_df A B x N dataframe where each row is a network in dyad form.
#' @param n Number of actors per network.
#' @param t Number of layers in the networks.
#' @param network_type Type of network to return. Options are "adj" for adjacency matrix (default), "igraph" for igraph object, or "network" for network object.
#'
#' @return A list of B lists of t adjacency matrices.
#' @export
dyads_to_matrix_list <- function(dyad_df, n, t, network_type = "adj") {
  B <- nrow(dyad_df)
  matrices <- list()
  if (t > 1) {
    for (i in 1:B) {
      res <- dyads_to_matrix_nd(dyad_df[i,], t, n)
      if (network_type == "igraph") {
        matrices[[i]] <- lapply(res, igraph::graph_from_adjacency_matrix)
      } else if (network_type == "network") {
        matrices[[i]] <- lapply(res, network::network, directed=T)
      } else if(network_type == "adj") {
        matrices[[i]] <- res      
      }
    }
  } else {
    #back to the uniplex case
    for (i in 1:B) {
      res <- dyads_to_matrix_1d(dyad_df[i,], n)
      if (network_type == "igraph") {
        matrices[[i]] <- igraph::graph_from_adjacency_matrix(res)
      } else if (network_type == "network") {
        matrices[[i]] <- network::network(res, directed=T)
      } else if(network_type == "adj") {
        matrices[[i]] <- res      
      }
    }
  }
  return(matrices)
}
