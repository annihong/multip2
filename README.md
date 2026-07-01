Sunbelt 2024 workshop: multip2 tutorial
================

## multip2

The multip2 package provides a Bayesian implementation of the multiplex
extension (Hong & Niezink, 2025) to the p2 (van Duijn et al., 2004)
model which includes cross-layer dyadic dependencies as fixed effects
and actor-specific dependencies as random effects. The influence of
covariates are also considered in the analysis of cross-sectional,
directed binary multiplex network data.

You can read more about the multiplex p2 model in the following paper:
[The Multiplex p2 Model: Mixed-Effects Modeling for Multiplex Social
Networks](https://arxiv.org/abs/2405.17707)

## Loading the package + visualization functions

``` r
#remotes::install_github("annihong/multip2") #uncomment if you have not installed the package already.
library(multip2)
library(igraph)
library(ggplot2)
```

``` r
plot_network_binary <- function(M, title, l) {
  colors <- c( "#264653", "#e76f51")
  types_full <- c("Not Governmental", "Government")
  net <- graph_from_adjacency_matrix(M)
  V(net)$govt <- actor_covar$govt
 plot(net, edge.arrow.size=.3,edge.color="#5b5b5b", vertex.label=NA, vertex.size=20,vertex.frame.color="#ffffff",
vertex.color=colors[actor_covar$govt + 1], main = title, layout = l)
  legend("bottomright", legend=types_full, col = colors , bty = "n", pch=20 ,
 pt.cex = 2, cex = 1)
}

plot_network_weighted <- function(M, title, l) {
  colors <- c( "#264653", "#e76f51")
  types_full <- c("Not Governmental", "Government")
  net <- graph_from_adjacency_matrix(M,  mode = "undirected", weighted = TRUE)
  V(net)$govt <- actor_covar$govt
 plot(net, edge.arrow.size=.3,edge.color="#5b5b5b", vertex.label=NA, vertex.size=20,vertex.frame.color="#ffffff",
vertex.color=colors[actor_covar$govt + 1], main = title, layout = l, edge.width=E(net)$weight)
  legend("bottomright", legend=types_full, col = colors , bty = "n", pch=20 ,
 pt.cex = 2, cex = 1)
}
```

## Load the data

ChemNet data is already stored in the package, now we are just going to
load it and construct lists and dataframe for inputting to the multiplex
p2 model.

``` r
data("ChemNet", package = "multip2")
dep_net <- list("political" = ChemNet$pol, "scientific" = ChemNet$sci)
dyad_covar <- list("prefsim" = ChemNet$prefsim, "influPercep" = ChemNet$infrep)
actor_covar <- data.frame(as.numeric(ChemNet[["govt"]]))
colnames(actor_covar) <- "govt"
```

## Visualize the data

Note that self loops will be ignored by the model, but it might be a
good practice to remove them.

``` r
l <- layout_nicely(graph_from_adjacency_matrix(dep_net[[1]])) 
plot_network_binary(dep_net$political, "Political Information", l)
plot_network_binary(dep_net$scientific, "Scientific Information", l)
plot_network_binary(dyad_covar$influPercep, "Perception of Influence", l)
plot_network_weighted(dyad_covar$prefsim, "Preference Similarity", l)
```

## Create an empty model with no covariates

``` r
m_empty <- Mp2Model(dep_net, dyad_covar = dyad_covar, actor_covar = actor_covar)
print(m_empty)
```

## Obtain samples from the prior distributions

Before fitting the model, we want to check if the default prior options
make sense. We will use the fit function to obtain draws from the priors
using stan.

``` r
m_empty_prior <- fit(m_empty, iter = 100, chains = 5, prior_sim = TRUE)
```

## Prior Predictive Checks

After obtaining samples of the parameter values, we will first convert
them to samples of the networks (using the generative model: multiplex
p2 model), then visualize it for inspection.

``` r
prior_draws <- extract_network_draws(m_empty_prior, 1000)
prior_simulated_network_checks(prior_draws)
```

## How to change priors

$\mu = 0$ means that
$\frac{P(i \to j|j \not \to i)}{P(i \not \to j|j \not \to i)} = \exp(0) = 1$,
that is it equally as likely for i to send a tie to j or not, *given
that there is no tie from j to i*. For example if we want the prior to
be such that it is *10 times* as likely for i to *not* send scientific
information to j than for i to send such information, *given that there
is no scientific tie from j to i*, in the *baseline*,  
$\frac{P(i \to j | j \not \to i}{P(i \not \to j| j \not \to i} = \frac{1}{10}$,
means that we need $\mu \sim \mathcal{N}(\log(\frac{1}{10}), \sigma)$.
Thus we change the mean of the prior distribution to be
$\log(\frac{1}{10}) = -2.30$.

But for practice we are going to set the mean to be -10 so there is a
bigger difference in the distributions.

``` r
m_empty_prior <- update_prior(m_empty, param="reciprocity", type="baseline", layer_lab = "scientific", mean=-10)
m_empty_prior <- fit(m_empty, iter = 50, chains = 5, prior_sim = TRUE, refresh=0)
#m_empty$model$prior
prior_draws <- extract_network_draws(m_empty_prior, 1000)
prior_simulated_network_checks(prior_draws)
```

``` r
rm(prior_draws)
rm(m_empty_prior)
```

## Fitting the empty model

Once we are happy with the priors for the model, we are ready to first
fit the empty model.

``` r
# m_fit <- fit(m_empty, iter=100)
load(file = "fittedMp2Models/m_empty_fit_1000.rda")
m_fit <- m_empty_fit_1000
rm(m_empty_fit_1000)

#summary.Mp2Model(m_fit)
knitr::kable(m_fit$fit_res$summary$fixed[,c(1,2,4,8,9,10)], digits = 3)
knitr::kable(m_fit$fit_res$summary$random[,c(1,2,4,8,9,10)], digits = 3)
```

## Convergence check

``` r
convergence_check(m_fit, params = "fixed", file = "empty_model_fixed.pdf")
convergence_check(m_fit, params = "random", file = "empty_model_random.pdf")
```

## Goodness-of-fit

After we assessed model convergence, we can assess how well this model
fits to the observed data. First we look at the uniplex GOF measures
from RSiena:

``` r
sim_nets <- extract_network_draws(m_fit, 1000)
simulated_network_checks(sim_nets, dep_net, "Triad_census")
simulated_network_checks(sim_nets, dep_net, "Indegree_distribution")
simulated_network_checks(sim_nets, dep_net, "Outdegree_distribution")
```

And then the multiplex GOFs.

``` r
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_baseline")
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_random")
```

## Missing values

We first create an example of networks with missing values. Make sure
the missing ties are coded as NA, then the package will take care of
encoding the ties as missing.

``` r
dep_net_missing <- dep_net
dep_net_missing$political[c(3,5,6),c(1,4,18)] <- NA
dep_net_missing$scientific[c(3,5,6),c(1,4,20)] <- NA

m_empty_missing <- Mp2Model(dep_net_missing, dyad_covar, actor_covar)
m_empty_missing_fit <- fit(m_empty_missing, iter=100, refresh=0)
```

The N_obs argument indicates the number of observed dyads in the data,
we can compare the number of observed dyads between the complete model
and the missing model.

``` r
m_empty_missing_fit$fit_res$stan_data$N_obs
m_fit$fit_res$stan_data$N_obs
```

Note that the GOF plots will look worse with models containing
missingness. This will be improved soon.

``` r
sim_nets_missing <- extract_network_draws(m_empty_missing_fit, 100)
simulated_network_checks(sim_nets_missing, dep_net, "Indegree_distribution")
simulated_network_checks(sim_nets_missing, dep_net, "multiplex_gof_baseline")
simulated_network_checks(sim_nets_missing, dep_net, "multiplex_gof_random")
rm(m_empty_missing)
rm(m_empty_missing_fit)
rm(sim_nets_missing)
```

# Adding in covariates

### Model 1: Baseline model with preference similarity:

``` r
m_fit <- m_empty
m_fit <- update_covar(m_fit, layer_1 = "political", density = "prefsim")
m_fit <- update_covar(m_fit, layer_1 = "scientific", density = "prefsim")
names(m_fit$model$covar)
# m_covar1_fit <- fit(m_covar1, iter=100)
load(file = "fittedMp2Models/m_covar1_3000.rda")
m_fit <- m_covar1_fit
rm(m_covar1_fit)
knitr::kable(m_fit$fit_res$summary$fixed[,c(1,2,4,8,9,10)], digits = 3)
knitr::kable(m_fit$fit_res$summary$random[,c(1,2,4,8,9,10)], digits = 3)
```

``` r
# knitr::kable(m_fit$fit_res$summary$fixed[,c(1,2,4,8,9,10)], digits = 3)
# knitr::kable(m_fit$fit_res$summary$random[,c(1,2,4,8,9,10)], digits = 3)
# convergence_check(m_fit, params = "fixed", file = "covar1_fixed.pdf")
# convergence_check(m_fit, params = "random", file = "covar1_random.pdf")
sim_nets <- extract_network_draws(m_fit, 1000)
simulated_network_checks(sim_nets, dep_net, "Triad_census")
simulated_network_checks(sim_nets, dep_net, "Indegree_distribution")
simulated_network_checks(sim_nets, dep_net, "Outdegree_distribution")
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_baseline")
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_random")
```

### Model 2: Adding perception of influence and govt:

Define the model

``` r
m_fit <- m_empty
m_fit <- update_covar(m_fit, layer_1 = "political", density = c("prefsim","influPercep"), receiver = "govt")
m_fit <- update_covar(m_fit, layer_1 = "scientific", density = c("prefsim","influPercep"), receiver = "govt")
names(m_fit$model$covar)
```

Fit the model

``` r
# m_covar2_fit <- fit(m_fit, iter=100)
load(file = "fittedMp2Models/m_covar2_3000.rda")
m_fit <- m_covar2_fit
rm(m_covar2_fit)
```

Check the model

``` r
knitr::kable(m_fit$fit_res$summary$fixed[,c(1,2,4,8,9,10)], digits = 3)
knitr::kable(m_fit$fit_res$summary$random[,c(1,2,4,8,9,10)], digits = 3)
# convergence_check(m_fit, params = "fixed", file = "covar2_fixed.pdf")
# convergence_check(m_fit, params = "random", file = "covar2_random.pdf")
sim_nets <- extract_network_draws(m_fit, 1000)
simulated_network_checks(sim_nets, dep_net, "Triad_census")
simulated_network_checks(sim_nets, dep_net, "Indegree_distribution")
simulated_network_checks(sim_nets, dep_net, "Outdegree_distribution")
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_baseline")
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_random")
```

### Model 3

How it’s your turn!

``` r
m_fit <- m_empty
#m_fit <- update_covar(m_fit, layer_1 = ..., layer_1 = ..., )
```

## Acknowledgement
This material is based upon work supported by the National Science Foundation under Grant Number SES-2020276. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
