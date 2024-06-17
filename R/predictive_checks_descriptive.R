
# comparing two networks
graph_similarity <- function(matrix1, matrix2){
  both <- sum((matrix1 == matrix2) & (matrix1 == 1), na.rm = T) 
  total_1 <- sum(matrix1 == 1, na.rm = T)
  only_1 <- sum((matrix1 != matrix2) & (matrix1 == 1), na.rm = T)
  total_2 <- sum(matrix2 == 1, na.rm = T)
  only_2 <- sum((matrix1 != matrix2) & (matrix2 == 1), na.rm = T)
  total <- (both + only_1 + only_2)
  jaccard <- both / total
  
  return(list(both=both, total_1=total_1, only_1=only_1, total_2 = total_2, only_2=only_2, total=total, jaccard=jaccard))
}

# given a list of networks, and the function (e.g., jaccard_cross_density)
# return a list of H jaccard indices
jaccard_list <- function(networks, func_J){
  pairs <- get_dyads(length(networks))
  res <- c()
  for(pair in pairs) {
    a <- pair[1]
    b <- pair[2]
    x1 <- networks[[a]]
    x2 <- networks[[b]]
    res <- c(res, func_J(x1, x2))
  }
  return(res)
}

#give two directed binary networks
#return the jaccard index between the first network and the transpose of the second network as a measurement of the crossnetwork reciprocity
jaccard_cross_reciprocity <- function(x1, x2){
  M1 <- as.matrix(igraph::as_adjacency_matrix(x1))
  M2 <- as.matrix(igraph::as_adjacency_matrix(x2))
  sim <- graph_similarity(M1, t(M2))
  return(sim$jaccard)
}

#give two directed binary networks
#return the jaccard index between the two as a measurement of the crossnetwork density
jaccard_cross_density <- function(x1, x2){
  M1 <- as.matrix(igraph::as_adjacency_matrix(x1))
  M2 <- as.matrix(igraph::as_adjacency_matrix(x2))
  sim <- graph_similarity(M1, M2)
  return(sim$jaccard)
}

descriptive_stats <- function(networks){
    network_names <- names(networks)
    networks <- lapply(networks, igraph::graph_from_adjacency_matrix)
    pair_names <- get_pair_names(network_names, ".")
    mu <- sapply(networks, igraph::edge_density)
    names(mu) <- paste(network_names, "density", sep = "_")
    rho <- sapply(networks, igraph::reciprocity)
    names(rho) <- paste(network_names, "reciprocity", sep = "_")
    cross_mu <- jaccard_list(networks, jaccard_cross_density)
    names(cross_mu) <- paste(pair_names, "Jaccard_cross_mu", sep = "_")
    cross_rho <- jaccard_list(networks, jaccard_cross_reciprocity)
    names(cross_rho) <- paste(pair_names, "Jaccard_cross_rho", sep = "_")

  res <- c(mu, rho, cross_mu, cross_rho)
  return(res)
}

descriptive_stats_Sigma_t3 <- function(networks, names_net = names(networks)){
  concat_strings <- function(x, y) {
    paste(x, y, sep = "_X_")
  }
    networks <- lapply(networks, igraph::graph_from_adjacency_matrix)
    indegree <- lapply(networks, igraph::degree, mode = "in", normalized=F)
    outdegree <- lapply(networks, igraph::degree, mode = "out", normalized=F)
    degrees <-list()
    for (i in 1:length(dep_net)) {
      degrees <- append(degrees, outdegree[i])
      degrees <- append(degrees, indegree[i])
    }
    # degrees <- c(outdegree[1], indegree[1], outdegree[2], indegree[2], outdegree[3], indegree[3])
    # #degrees <- c(outdegree[1], indegree[1], outdegree[2], indegree[2])
    var_names <- paste(rep(names_net, each = 2), c("od", "id"), sep = "_")
    max_length <- max(sapply(degrees, length))
    degrees <- do.call(cbind, lapply(degrees, function(x) c(x, rep(0,max_length - length(x)))))
    covariance <- cov(degrees)
    keep = !lower.tri(covariance)
    res = covariance[keep]
    var_names  <- outer(var_names, var_names, concat_strings)[keep]
    names(res) <- var_names
    
    return(res)
}
