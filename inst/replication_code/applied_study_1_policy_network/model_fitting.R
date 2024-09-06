library(multip2)
file_path <- "/home/annihong/projects/simres/fitted_models/"

WARMUP = 2000
CHAINS = 4

# load data
data("ChemNet", package = "multip2") # the data is stored in the package, see data_raw/ChemNet.R for the script to clean the raw data
dep_net <- list("PO" = ChemNet$pol, "SC" = ChemNet$sci, "PE" = ChemNet$infrep)
dyad_covar <- list("prefsim" = ChemNet$prefsim)
actor_covar <- data.frame(as.numeric(ChemNet[["govt"]]))
colnames(actor_covar) <- "govt"

# fit the model, empty model is fitted first:
m_1 <- Mp2Model(dep_net, dyad_covar, actor_covar)
m_1 <- fit(m_1, chains = CHAINS, warmup = WARMUP, thin = 1, iter = WARMUP*2)

# model 2 includes the covariates info
# configure the data
m_2 <- Mp2Model(dep_net, dyad_covar, actor_covar)
# specify the covariates
m_2 <- update_covar(m_2, layer_1 = "PO", density = "prefsim", receiver = "govt", sender = "govt")
m_2 <- update_covar(m_2, layer_1 = "SC", density = "prefsim", receiver = "govt", sender = "govt")
m_2 <- update_covar(m_2, layer_1 = "PE", density = "prefsim", receiver = "govt", sender = "govt")
# fit the model
m_2 <- fit(m_2, chains = CHAINS, warmup = WARMUP, thin = 1, iter = WARMUP*2)

results <- list(m_1 = m_1, m_2 = m_2)
# save the fitted models
saveRDS(results, file = paste0(file_path, "study_1_models.Rds"))