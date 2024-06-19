n = 5
networks = replicate(4, matrix(data=sample(c(0,1), n^2, replace=TRUE, prob=c(0.7,0.3)), nrow=n), simplify = FALSE) #self loops are ignored
names(networks) <-  c("network1", "network2", "network_covariate1",      "network_covariate2")
actor_data = data.frame(actor_attr1=1:n, actor_attr2=1:n)

dep_net_uni <- networks[c("network1")]
dep_net<- networks[c("network1", "network2")]
dyad_covar <- networks[c("network_covariate1", "network_covariate2")]
actor_covar <- actor_data[c("actor_attr1", "actor_attr2")]

m_empty_uni <- Mp2Model(dep_net_uni, dyad_covar = dyad_covar, actor_covar = actor_data)
m_empty <- Mp2Model(dep_net, dyad_covar = dyad_covar, actor_covar = actor_data)


test_that("Uniplex empty model works as expected", {
  # Create test data
  dep_net <- list(matrix(0, nrow = 10, ncol = 10))
  dyad_covar <- list(matrix(0, nrow = 10, ncol = 10))
  actor_covar <- data.frame(matrix(0, nrow = 10, ncol = 5))

  # Call the function
  model <- Mp2Model(dep_net, dyad_covar, actor_covar)
  names(dep_net) <- "network1" 

  # Check the attributes of the model
  expect_equal(model$t, 1)
  expect_equal(model$H, 0)
  expect_equal(model$n, 10)
  expect_equal(model$dep_lab, "network1")
  expect_equal(model$pair_lab, NULL)
  expect_equal(model$dyad_covar_lab, names(dyad_covar))
  expect_equal(model$actor_covar_lab, names(actor_covar))
  expect_equal(model$data$dep_net, dep_net)
  expect_equal(model$data$dyad_covar, dyad_covar)
  expect_equal(model$data$actor_covar, actor_covar)
  expect_equal(model$model$covar, list())
  expect_equal(model$model$prior, default_prior_empty_mdl(names(dep_net), NULL))
  expect_equal(model$model$full_prior, NULL)
  expect_equal(model$fit_res$stan_data, NULL)
  expect_equal(model$fit_res$stan_fit, NULL)
  expect_equal(model$fit_res$summary, NULL)
  expect_equal(model$fit_res$par_labels, NULL)
})

test_that("Multiplex empty model works as expected", {
  # Create test data
  dep_net <- list("net1" = matrix(0, nrow = 10, ncol = 10), "net2"=matrix(0, nrow = 10, ncol = 10))
  dyad_covar <- list(matrix(0, nrow = 10, ncol = 10))
  actor_covar <- data.frame(matrix(0, nrow = 10, ncol = 5))

  # Call the function
  model <- Mp2Model(dep_net, dyad_covar, actor_covar)

  # Check the attributes of the model
  expect_equal(model$t, 2)
  expect_equal(model$H, 1)
  expect_equal(model$n, 10)
  expect_equal(model$dep_lab, names(dep_net))
  expect_equal(model$pair_lab, get_pair_names(names(dep_net)))
  expect_equal(model$dyad_covar_lab, names(dyad_covar))
  expect_equal(model$actor_covar_lab, names(actor_covar))
  expect_equal(model$data$dep_net, dep_net)
  expect_equal(model$data$dyad_covar, dyad_covar)
  expect_equal(model$data$actor_covar, actor_covar)
  expect_equal(model$model$covar, list())
  expect_equal(model$model$prior, default_prior_empty_mdl(names(dep_net), get_pair_names(names(dep_net))))
  expect_equal(model$model$full_prior, NULL)
  expect_equal(model$fit_res$stan_data, NULL)
  expect_equal(model$fit_res$stan_fit, NULL)
  expect_equal(model$fit_res$summary, NULL)
  expect_equal(model$fit_res$par_labels, NULL)
})


test_that("Mp2Model throws an error with invalid input", {
  # Create test data with invalid input
  dep_net <- list(matrix(0, nrow = 10, ncol = 10))
  dyad_covar <- list(matrix(0, nrow = 5, ncol = 5))  # Invalid dimensions
  actor_covar <- data.frame(matrix(0, nrow = 10, ncol = 5))

  # Call the function and expect an error
  expect_error(Mp2Model(dep_net, dyad_covar, actor_covar), "is_valid")
})



test_that("update_covar correctly updates within-layer covariates uniplex and multiplex", {
  model_obj <- m_empty_uni
  # Update covariates for layer1
  expect_error(update_covar(model_obj, layer_1 = "network1", density = c("network_covariate1", "network_covariate2"), reciprocity = c("actor_attr2")))
  updated_model_obj <- update_covar(model_obj, layer_1 = "network1", density = c("network_covariate1","network_covariate2"), receiver = c("actor_attr2"))
  # Check if the covariates are correctly updated
  expect_equal(names(updated_model_obj$model$covar), c("density_network_covariate1_network1", "density_network_covariate2_network1", "receiver_actor_attr2_network1"))

  model_obj <- m_empty
  # Update covariates for layer1
  updated_model_obj <- update_covar(model_obj, layer_1 = "network1", density = c("network_covariate1", "network_covariate2"), receiver = c("actor_attr2"))
  
  # Check if the covariates are correctly updated
  expect_equal(names(updated_model_obj$model$covar), c("density_network_covariate1_network1", "density_network_covariate2_network1", "receiver_actor_attr2_network1"))

  updated_model_obj <- update_covar(model_obj, layer_1 = "network2", density = c("network_covariate1", "network_covariate2"), receiver = c("actor_attr2"))
  
  # Check if the covariates are correctly updated
  expect_equal(names(updated_model_obj$model$covar), c("density_network_covariate1_network2", "density_network_covariate2_network2", "receiver_actor_attr2_network2"))
})

test_that("update_covar correctly updates cross-layer covariates", {
  model_obj <- m_empty_uni
  # Update cross-layer covariates for uniplex network should return an error
  
  expect_error(update_covar(model_obj, layer_1 = "network1", layer_2 = "network2",  cross_density = c("network_covariate1","network_covariate2"), receiver = c("actor_attr2")))

  model_obj <- m_empty

   expect_error(update_covar(model_obj, layer_1 = "network1", cross_density = c("network_covariate1","network_covariate2"), receiver = c("actor_attr2")))
  # multiplex cross density
  updated_model_obj <- update_covar(model_obj, layer_1 = "network1", layer_2 = "network2", cross_density = c("network_covariate1","network_covariate2"), receiver = c("actor_attr2"))
  
  # Check if the covariates are correctly updated
  expect_equal(names(updated_model_obj$model$covar), c("cross_density_network_covariate1_network1:network2", "cross_density_network_covariate2_network1:network2"))

  updated_model_obj <- update_covar(model_obj, layer_1 = "network2", layer_2 = "network1", cross_reciprocity = c("network_covariate1","network_covariate2"), receiver = c("actor_attr2"))

    expect_equal(names(updated_model_obj$model$covar), c("cross_reciprocity_network_covariate1_network1:network2", "cross_reciprocity_network_covariate2_network1:network2"))
})

test_that("update_covar returns an error when layer is not found", {
  # Create a dummy Mp2Model object
  model_obj <- m_empty
  # Update covariates for a non-existent layer
  expect_error(update_covar(model_obj, layer_1 = "layer4", density = c("covar1")))
})

test_that("update_covar returns an error when covariate is not found", {
    model_obj <- m_empty
  # Update covariates with a non-existent covariate
  expect_error(update_covar(model_obj, layer_1 = "network1", density = c("covar4")))
})



# Test for baseline parameter with mean and sd
test_that("update_prior sets baseline parameter with mean and sd", {
    model_obj <- m_empty_uni

    updated_model_obj <- model_obj
    expect_equal(as.vector(updated_model_obj$model$prior$rho_mean_prior), c(0))
    updated_model_obj <- update_prior(updated_model_obj, "density", "baseline", mean = 10, sd = 0.5)
    updated_model_obj <- update_prior(updated_model_obj, "reciprocity", "baseline", mean = 10, sd = 0.5)
    expect_equal(as.vector(updated_model_obj$model$prior$rho_mean_prior), c(10))
    expect_equal(as.vector(updated_model_obj$model$prior$mu_mean_prior), c(10))
    expect_equal(as.vector(updated_model_obj$model$prior$rho_sd_prior), c(0.5))
    expect_equal(as.vector(updated_model_obj$model$prior$mu_sd_prior), c(0.5))
    
    
    model_obj <- m_empty

    updated_model_obj <- model_obj
    expect_equal(as.vector(updated_model_obj$model$prior$rho_mean_prior), c(0, 0))
    updated_model_obj <- update_prior(updated_model_obj, "density", "baseline", layer_lab="network1", mean = 10, sd = 0.5)
    updated_model_obj <- update_prior(updated_model_obj, "reciprocity", "baseline", mean = 10, sd = 0.5)
    expect_equal(as.vector(updated_model_obj$model$prior$rho_mean_prior), c(10, 10))
    expect_equal(as.vector(updated_model_obj$model$prior$mu_mean_prior), c(10, 0))
    expect_equal(as.vector(updated_model_obj$model$prior$rho_sd_prior), c(0.5, 0.5))
    expect_equal(as.vector(updated_model_obj$model$prior$mu_sd_prior), c(0.5, 10))
})
  

# Test for random parameter with eta, alpha, and beta
test_that("update_prior sets random parameter with eta, alpha, and beta", {
    model_obj <- m_empty_uni

    updated_model_obj <- model_obj
    expect_equal(updated_model_obj$model$prior$LJK_eta_prior, 2)
    expect_equal(updated_model_obj$model$prior$scale_alpha_prior, 3)
    expect_equal(updated_model_obj$model$prior$scale_beta_prior, 50)
    updated_model_obj <- update_prior(updated_model_obj, "LKJ", "random", eta = 0.5, beta = 0.9)
    expect_equal(updated_model_obj$model$prior$LJK_eta_prior, 0.5)
    expect_equal(updated_model_obj$model$prior$scale_alpha_prior, 3)
    expect_equal(updated_model_obj$model$prior$scale_beta_prior, 0.9)

    model_obj <- m_empty

    updated_model_obj <- model_obj
    expect_equal(updated_model_obj$model$prior$LJK_eta_prior, 2)
    expect_equal(updated_model_obj$model$prior$scale_alpha_prior, 3)
    expect_equal(updated_model_obj$model$prior$scale_beta_prior, 50)
    updated_model_obj <- update_prior(updated_model_obj, "LKJ", "random", eta = 0.5, beta = 0.9)
    expect_equal(updated_model_obj$model$prior$LJK_eta_prior, 0.5)
    expect_equal(updated_model_obj$model$prior$scale_alpha_prior, 3)
    expect_equal(updated_model_obj$model$prior$scale_beta_prior, 0.9)
})

# Test for covariate parameter with mean and sd
test_that("update_prior sets covariate parameter with mean and sd", {
    model_obj <- m_empty_uni
    model_obj <- update_covar(model_obj, layer_1 = "network1", density = c("network_covariate1", "network_covariate2"), receiver = c("actor_attr2"))
    updated_model_obj <- model_obj
    updated_model_obj <- update_prior(model_obj, "density", "covariate", covar_lab = "network_covariate1", mean = 10, sd = 0.5)
    expect_equal(updated_model_obj$model$covar$density_network_covariate1_network1$mean_prior, 10)
    expect_equal(updated_model_obj$model$covar$density_network_covariate1_network1$sd_prior, 0.5)

    model_obj <- m_empty
    model_obj <- update_covar(model_obj, layer_1 = "network1", density = c("network_covariate1", "network_covariate2"), receiver = c("actor_attr2"))
    model_obj <- update_covar(model_obj, layer_1 = "network1", layer_2="network2", cross_density = c("network_covariate1", "network_covariate2"))
    updated_model_obj <- model_obj
    updated_model_obj <- update_prior(updated_model_obj, "density", "covariate", covar_lab = "network_covariate1", layer_lab="network1", mean = 10, sd = 0.5)
    expect_equal(updated_model_obj$model$covar$density_network_covariate1_network1$mean_prior, 10)
    expect_equal(updated_model_obj$model$covar$density_network_covariate1_network1$sd_prior, 0.5)
    fitted_model <- fit(updated_model_obj, iter = 10, chains = 1, refresh=0)

})




## tests for fit.Mp2Model
# Test the fit function
test_that("fit function performs estimation using Stan for empty uniplex and multiplex models", {
  model_obj <- m_empty_uni
  fitted_model <- fit(model_obj, iter=10, chains=1, refresh=0)
  expect_true(class(fitted_model) == "Mp2Model")
  
  # Check if the fitted_model contains the fitted Stan model and parameter labels
  expect_true(!is.null(fitted_model$fit_res$stan_fit))
  expect_true(!is.null(fitted_model$fit_res$par_labels))

  model_obj <- m_empty
  fitted_model <- fit(model_obj, iter=10, chains=1, refresh=0)
  fitted_model <- fit(model_obj, iter=10, chains=1, refresh=0)
  fitted_model <- fit(model_obj, iter=10, chains=1, refresh=0)
  expect_true(class(fitted_model) == "Mp2Model")




  
  # Check if the fitted_model contains the fitted Stan model and parameter labels
  expect_true(!is.null(fitted_model$fit_res$stan_fit))
  expect_true(!is.null(fitted_model$fit_res$par_labels))
})

# Test for covariate parameter with mean and sd
test_that("update_prior sets covariate parameter with mean and sd after fitting", {
  
    model_obj <- m_empty
  
  
    model_obj <- update_covar(model_obj, layer_1 = "network1", density = c("network_covariate1", "network_covariate2"), receiver = c("actor_attr2"))
    model_obj <- update_covar(model_obj, layer_1 = "network1", layer_2="network2", cross_density = c("network_covariate1", "network_covariate2"))
    updated_model_obj <- model_obj
    updated_model_obj <- update_prior(updated_model_obj, "density", "covariate", covar_lab = "network_covariate1", layer_lab="network1", mean = 10, sd = 0.5)
    expect_equal(updated_model_obj$model$covar$density_network_covariate1_network1$mean_prior, 10)
    expect_equal(updated_model_obj$model$covar$density_network_covariate1_network1$sd_prior, 0.5)
    fitted_model <- fit(updated_model_obj, iter=10, chains=1, refresh=0)
    expect_true(!is.null(fitted_model$model$full_prior))
    expect_true(fitted_model$model$full_prior$mu_covariates_mean_prior["density_network_covariate1_network1"] == 10)
    expect_true(fitted_model$model$full_prior$mu_covariates_sd_prior["density_network_covariate1_network1"] == 0.5)
    expect_true(fitted_model$model$full_prior$beta_covariates_mean_prior["receiver_actor_attr2_network1"] == 0)
    expect_true(fitted_model$model$full_prior$beta_covariates_sd_prior["receiver_actor_attr2_network1"] == 10/sd(fitted_model$data$actor_covar$actor_attr2))

})

