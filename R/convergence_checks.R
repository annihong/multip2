#' Check convergence of a multiplex network model
#'
#' This function checks the convergence of a fitted Stan model for a multiplex network using the ggmcmc package.
#'
#' @param multip2fit the fitted MultiP2Fit object.
#' @param params A character string indicating which parameters to plot. Must be one of "fixed" (all the fixed parameters), "random"(actor random effects), or "custom" (use the custom_params argument)
#' @param custom_params A character string indicating a custom regular expression to match parameter names to plot, as given by a character vector or a regular expression, i.e., the "family" parameter in ggmcmc::ggs.
#' @param plot A character vector containing the names of the desired plots.m Default is c("traceplot", "running", or "geweke")
#' @param plot_all A logical indicating whether to plot all types of plots. If TRUE, the plot argument is ignored, instead all plots are generated.
#' @param file A character string indicating the output file name for the plots, default = "ggmcmc-output.pdf"
#' @return NULL
#' @export
convergence_check <- function(multip2fit, params, custom_params="", plot=c("traceplot", "running", "geweke", "density"), plot_all = FALSE, file = "ggmcmc-output.pdf") {
    fit <- multip2fit$stan_fit
    par_labels <- multip2fit$par_labels

    if (params == "fixed") {
        S <- ggmcmc::ggs(fit, family = "^mu|^rho|^cross|fixed", par_labels = par_labels, keep_original_order=TRUE, sort=FALSE)
    } else if (params == "random") {
        S <- ggmcmc::ggs(fit, family = "^Sigma", par_labels = par_labels, keep_original_order=TRUE)
    } else if (params == "custom") {
        S <- ggmcmc::ggs(fit, family = custom_params, par_labels = par_labels, keep_original_order=TRUE)
    } else {
        stop("params must be one of 'fixed', 'random', or 'custom'")
    }

    if (plot_all) {
        plot = NULL
    } 

    ggmcmc::ggmcmc(S, file = file, param_page=6, plot=plot)
}
