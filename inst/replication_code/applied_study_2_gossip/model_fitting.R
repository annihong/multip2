library(multip2)
library(doSNOW)
# library(tcltk)
file_path <- "/home/annihong/projects/simres/fitted_models/"
fitted_ids <- sub("^([0-9]+).*", "\\1", grep("*_study_2_models.Rds", list.files(file_path), value=TRUE))
#reading in data description to determine which classes we will use
path_to_package <- "/home/annihong/projects/multip2/"

data_path <- "/home/annihong/projects/multiplex-social-networks/stan/coding/gossip_perception/Data/Raw/GossipBully/"
file_path <- "/home/annihong/projects/simres/fitted_models/"

filtered_df <- read.csv(paste0(path_to_package, "data-raw/gossip/filtered_df_no_missing_wave_4.csv")) 
to_fit_ids <- filtered_df$class_id
class_ids <- setdiff(to_fit_ids, fitted_ids)


#output path, change to your own 

WARMUP = 2000
CHAINS = 4
WAVE=4
NUM_CORE = CHAINS * length(class_ids)
#class_begin = 19 #doing it in chuncks
#class_end = 39

summary <- read.csv(paste0(path_to_package, "data-raw/gossip/summary.csv"))

# pb <- tkProgressBar(max=length(class_ids))
# progress <- function(n) setTkProgressBar(pb, n)
# opts <- list(progress=progress)

cl <- parallel::makeCluster(NUM_CORE, outfile="/home/annihong/projects/LogGossip.txt")  
registerDoSNOW(cl)  # Register the cluster for use with foreach

res <- foreach(i=1:length(class_ids), .packages="multip2") %dopar% {
    class_id = class_ids[[i]]
    load(file.path(data_path, paste0(class_id,"_wave_",WAVE,"_cleaned.RData")))

    dep_net = list(sts =cleaned_data$sts, tst = cleaned_data$tst)
    is_female = cleaned_data$is_female
    actor_covar <- data.frame(is_female = is_female[rownames(dep_net$sts)])
    dyad_covar <- list(same_gender = array(outer(actor_covar$is_female, actor_covar$is_female, "=="), dim = c(nrow(actor_covar), nrow(actor_covar))))

    # fit the model, empty model is fitted first:
    m_1 <- Mp2Model(dep_net, dyad_covar, actor_covar)
    m_1 <- fit(m_1, chains = CHAINS, warmup = WARMUP, thin = 1, iter = WARMUP*2, network_sim = FALSE)

    # model 2 includes the covariates info
    # configure the data
    m_2 <- Mp2Model(dep_net, dyad_covar, actor_covar)
    # specify the covariates
    m_2 <- update_covar(m_2, layer_1 = "sts", receiver = "is_female", sender = "is_female")
    m_2 <- update_covar(m_2, layer_1 = "tst", receiver = "is_female", sender = "is_female")
    m_2 <- update_covar(m_2, layer_1 = "sts", layer_2 = "tst", cross_density = "same_gender", cross_reciprocity = "same_gender")
    # fit the model
    m_2 <- fit(m_2, chains = CHAINS, warmup = WARMUP, thin = 1, iter = WARMUP*2, network_sim = FALSE)

    results <- list(m_1 = m_1, m_2 = m_2)
    # save the fitted models
    saveRDS(results, file = paste0(file_path, class_id, "_study_2_models.Rds"))
    rm(m_1)
    rm(m_2)
    rm(results)
}

stopCluster(cl)



