
#' Format covariates data for stan estimation
#' 
#' @param outcome vector: names of the layers of networks in the network_data as the outcome multiplex network
#' @param network_data list of p + t matrices of dim n x n : network outcomes and network covariates, p is the number of network covariates and t is the number of layers in the multiplex network outcome. 
#' @param actor_data n x p data frame: p covariates of length of the number of actors (n) 
#' @return a MultiP2Fit object 
#' @export
#'
#' @examples
#' networks = replicate(4, matrix(data=1:9, nrow=3), simplify = FALSE)
#' names(networks) <-  c("network1", "network2", "network_covariate1",      "network_covariate2")
#' MultiP2Fit(outcome = c("network1", "network2"), 
#' network_data = networks, 
#' actor_data = data.frame(actor_attr1=1:3, actor_attr2=1:3))
#constructor
MultiP2Fit <- function(outcome = "", network_data=list(), actor_data=data.frame()) {
    stopifnot(is.character(outcome))
    stopifnot(is.list(network_data))
    stopifnot(is.data.frame(actor_data))

    newMultiP2Fit = structure(
                            list(),
                            outcome=outcome,
                            network_data=network_data,
                            actor_data=actor_data,
                            class = "MultiP2Fit"
                            )
    return(newMultiP2Fit)
}
