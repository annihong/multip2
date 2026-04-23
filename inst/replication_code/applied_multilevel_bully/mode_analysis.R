##### obtain class-size  ######
setwd("./multip2")
devtools::load_all()
library(rstan)
options(mc.cores = 10)

library(ggplot2)
library(dplyr)
library(patchwork)  # for combining plots
library(rstan)
library(tidytext)

# Pretty theme for consistency
custom_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(10, 15, 10, 10),
    legend.position = "none",
  )


rstan::rstan_options(auto_write = TRUE)
source("./inst/multilevel_extension/helper.R")
data_path = "/home/annihong/projects/multiplex-social-networks/stan/coding/gossip_perception/Data/Raw/GossipBully"

######gossip data##########

filtered_df <- read.csv("data-raw/bully/filtered_df_no_missing_wave_4_bully1.csv")
f_gossip <- read.csv("/home/annihong/projects/multip2/data-raw/gossip/filtered_df_no_missing_wave_4.csv")
class_ids = filtered_df$class_id
summary_df <- read.csv("/home/annihong/projects/multip2/data-raw/gossip/summary.csv")
summary_df <- summary_df %>% rename(class_id = class) %>% filter(class_id %in% class_ids, wave == 4) 
summary_df <- summary_df %>% left_join(filtered_df, by = c("class_id"))
WAVE=4
#######extension to L networks ######
L = length(class_ids)
t = 2
H = 1
layer_covar = TRUE
group_covar = TRUE

#########group covar###########
roma_percent <- sapply(class_ids, create_classlevel_covar, var_name = "roma", summary_df = summary_df)
ses_percent <- sapply(class_ids, create_classlevel_covar, var_name = "ses", summary_df = summary_df)
class_size <- sapply(class_ids, create_classlevel_covar, var_name = "num_present", summary_df = summary_df)

class_size_id <- paste(class_ids, class_size, sep = " : ")

filtered_df$roma_percent <- roma_percent
filtered_df$ses_percent <- ses_percent
filtered_df$class_size <- class_size
cor_plot <- filtered_df %>% select(bully1_density, bully1_rev_density, bully1_reciprocity, bully1_rev_reciprocity, sender_psyche_jaccard, target_percep_jaccard,class_size, roma_percent, ses_percent, prop_female) %>% 
cor() %>% 
melt() %>% 
ggplot( aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)

# ggsave(glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/plot/bully1_cor_plot.png"), plot = cor_plot, width = 12, height = 10, dpi = 300)  
################################

file_name <- "bully1_2000_full"
res <- readRDS(file = glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/code/{file_name}.Rds"))

s <-  summary(res)$summary


params <- c("mu_PS","mu_fixed_coef", "rho_PS", "alpha_fixed_coef", "beta_fixed_coef", "cross_mu_PS", "cross_rho_PS", "cross_rho_fixed_coef", "mu_group_coef")
df <- data.frame()
for(param in params){
  print(param)
  result <- s[grep(glue::glue("^{param}"), rownames(s)),c(1,4,8,9,10), drop=FALSE] %>% data.frame()

  df <- df %>% bind_rows(result)
}
df

result <- s[grep(glue::glue("Corr\\["), rownames(s)),c(1), drop=FALSE][,1] %>% matrix(nrow=4) %>% as.data.frame()
colnames(result) <- c("sender:bully", "receiver:bully", "sender:victim", "receiver:victim")
rownames(result) <- c("sender:bully", "receiver:bully", "sender:victim", "receiver:victim")
result

lower <- s[grep(glue::glue("Corr\\["), rownames(s)),c(4), drop=FALSE][,1] %>% matrix(nrow=4) %>% as.data.frame()
upper <- s[grep(glue::glue("Corr\\["), rownames(s)),c(8), drop=FALSE][,1] %>% matrix(nrow=4) %>% as.data.frame()
result <- lower*upper > 0
colnames(result) <- c("sender:bully", "receiver:bully", "sender:victim", "receiver:victim")
rownames(result) <- c("sender:bully", "receiver:bully", "sender:victim", "receiver:victim")
result

## Intercepts
# param = "mu"

for (param in c("mu", "rho")) {
  df1 <- s[grep(glue::glue("^{param}_random_PS.*1]$"), rownames(s)),c(1,4,8,9,10)] %>% data.frame()
  df1$class <- class_size_id
  df1$layer <- "Bully1"
  df2 <- s[grep(glue::glue("^{param}_random_PS.*2]$"), rownames(s)),c(1,4,8,9,10)] %>% data.frame()
  df2$class <- class_size_id
  df2$layer <- "Bully1_rev"
  df <- rbind(df1, df2)
  df$type = "random_intercept"

  overall_df <- s[grep(glue::glue("^{param}_PS"), rownames(s)),c(1,4,8,9,10)] %>% data.frame() %>% mutate(
    class = "Overall",
    layer = c("Bully1", "Bully1_rev"),
    type = "overall_intercept"
  )

  df <-  rbind(df, overall_df)
  colnames(df) <- c("mean", "lower", "upper", "n_eff", "Rhat", "class", "layer", "type")



  # Intercepts
  int_plt <- df %>%
    mutate(class  = as.factor(class)) %>%
    ggplot(aes(x = mean, y = reorder_within(class, mean, layer), fill = type, color = type)) +
    geom_errorbarh(aes(xmin = lower,
                      xmax = upper, 
                  color = type),
                  height = 0.2) +
    geom_point(size = 3,
              shape = 21,
              color = "black") +
    geom_vline(xintercept = 0,
              color = "red",
              size = 0.8,
              linetype = "dashed") +
  #   scale_x_continuous(breaks = seq(-60, 60, 20), limits = c(-20, 10)) +
    labs(x = "Intercept Estimate",
        y = "Classroom",
        title = param) +
    custom_theme +
    facet_wrap(~ layer, ncol = 1, scales = "free_y") +  # Allow free y-axis scales
    scale_y_reordered()
  int_plt

  ggsave(glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/plot/{param}_intercept_plot_{file_name}.png"), plot = int_plt, width = 8, height = 10, dpi = 300)
}

####### slope ########
for (param in c("alpha", "beta")){ 
   s[grep(glue::glue("^{param}"), rownames(s)),c(1,4,8,9,10)]

  df1 <- s[grep(glue::glue("^{param}_coef\\[1,"), rownames(s)),c(1,4,8,9,10)] %>% data.frame()
  df1$class <- class_size_id
  df1$layer <- "bully1"
  df2 <- s[grep(glue::glue("^{param}_coef\\[3,"), rownames(s)),c(1,4,8,9,10)] %>% data.frame()
  df2$class <- class_size_id
  df2$layer <- "bully1_rev"
  df <- rbind(df1, df2)
  df$type = "random_slope"

  overall_df <- s[grep(glue::glue("^{param}_fixed_coef\\[1|^{param}_fixed_coef\\[3"), rownames(s)),c(1,4,8,9,10)] %>% data.frame() %>% mutate(
    class = "Overall",
    layer = c("bully1", "bully1_rev"),
    type = "overall_intercept"
  )

  df <-  rbind(df, overall_df)
  colnames(df) <- c("mean", "lower", "upper", "n_eff", "Rhat", "class", "layer", "type")



  # Intercepts
  slp_plt <- df %>%
    mutate(class  = as.factor(class)) %>%
    ggplot(aes(x = mean, y = reorder_within(class, mean, layer), fill = type, color = type)) +
    geom_errorbarh(aes(xmin = lower,
                      xmax = upper, 
                  color = type),
                  height = 0.2) +
    geom_point(size = 3,
              shape = 21,
              color = "black") +
    geom_vline(xintercept = 0,
              color = "red",
              size = 0.8,
              linetype = "dashed") +
  #   scale_x_continuous(breaks = seq(-60, 60, 20), limits = c(-20, 10)) +
    labs(x = "Slope Estimate",
        y = "Classroom",
        title = glue::glue("Is_female: {param}")) +
    custom_theme +
    facet_wrap(~ layer, ncol = 1, scales = "free_y") +  # Allow free y-axis scales
    scale_y_reordered()
  slp_plt

  ggsave(glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/plot/{param}_slope_plot_female{file_name}.png"), plot = slp_plt, width = 8, height = 10, dpi = 300)
}

####### slope ########
for (param in c("alpha", "beta")){ 
   s[grep(glue::glue("^{param}"), rownames(s)),c(1,4,8,9,10)]

  df1 <- s[grep(glue::glue("^{param}_coef\\[2,"), rownames(s)),c(1,4,8,9,10)] %>% data.frame()
  df1$class <- class_size_id
  df1$layer <- "bully1"
  df2 <- s[grep(glue::glue("^{param}_coef\\[4,"), rownames(s)),c(1,4,8,9,10)] %>% data.frame()
  df2$class <- class_size_id
  df2$layer <- "bully1_rev"
  df <- rbind(df1, df2)
  df$type = "random_slope"

  overall_df <- s[grep(glue::glue("^{param}_fixed_coef\\[2|^{param}_fixed_coef\\[4"), rownames(s)),c(1,4,8,9,10)] %>% data.frame() %>% mutate(
    class = "Overall",
    layer = c("bully1", "bully1_rev"),
    type = "overall_intercept"
  )

  df <-  rbind(df, overall_df)
  colnames(df) <- c("mean", "lower", "upper", "n_eff", "Rhat", "class", "layer", "type")



  # Intercepts
  slp_plt <- df %>%
    mutate(class  = as.factor(class)) %>%
    ggplot(aes(x = mean, y = reorder_within(class, mean, layer), fill = type, color = type)) +
    geom_errorbarh(aes(xmin = lower,
                      xmax = upper, 
                  color = type),
                  height = 0.2) +
    geom_point(size = 3,
              shape = 21,
              color = "black") +
    geom_vline(xintercept = 0,
              color = "red",
              size = 0.8,
              linetype = "dashed") +
  #   scale_x_continuous(breaks = seq(-60, 60, 20), limits = c(-20, 10)) +
    labs(x = "Slope Estimate",
        y = "Classroom",
        title = glue::glue("Low_ses: {param}")) +
    custom_theme +
    facet_wrap(~ layer, ncol = 1, scales = "free_y") +  # Allow free y-axis scales
    scale_y_reordered()
  slp_plt

  ggsave(glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/plot/{param}_slope_plot_ses{file_name}.png"), plot = slp_plt, width = 8, height = 10, dpi = 300)
}




## cross_layer effects

for (param in c("cross_mu", "cross_rho")) {  
  df1 <- s[grep(glue::glue("^{param}_random_PS.*1]$"), rownames(s)),c(1,4,8,9,10)] %>% data.frame()
  df1$class <- class_size_id
  df1$layer <- "bully1 x bully1_rev"

  df <- df1
  df$type = "random_intercept"
  colnames(df) <- c("mean", "lower", "upper", "n_eff", "Rhat", "class", "layer", "type")
  overall_df <- s[grep(glue::glue("^{param}_PS"), rownames(s)),c(1,4,8,9,10)]

  # s[grep(glue::glue("^{param}\\[1\\]"), rownames(s)),c(1,4,8,9,10)]

  # s[grep(glue::glue("^{param}_random\\["), rownames(s)),c(1,4,8,9,10)] %>% data.frame()

  overall_df <- as.data.frame(t(overall_df)) %>% mutate(
    class = "Overall",
    layer =  "bully1 x bully1_rev",
    type = "overall_intercept"
  )
  colnames(overall_df) <- c("mean", "lower", "upper", "n_eff", "Rhat", "class", "layer", "type")
  df <-  rbind(df, overall_df)

  # Intercepts
  int_plt_cross <- df %>%
    mutate(class  = as.factor(class)) %>%
    ggplot(aes(x = mean, y = reorder_within(class, mean, layer), fill = type, color = type)) +
    geom_errorbarh(aes(xmin = lower,
                      xmax = upper, 
                  color = type),
                  height = 0.2) +
    geom_point(size = 3,
              shape = 21,
              color = "black") +
    geom_vline(xintercept = 0,
              color = "red",
              size = 0.8,
              linetype = "dashed") +
  #   scale_x_continuous(breaks = seq(-60, 60, 20), limits = c(-20, 10)) +
    labs(x = "Intercept Estimate",
        y = "Classroom",
        title = param) +
    custom_theme +
    facet_wrap(~ layer, ncol = 1, scales = "free_y") +  # Allow free y-axis scales
    scale_y_reordered()
  int_plt_cross

  ggsave(glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/plot/{param}_intercept_plot_{file_name}.png"), plot = int_plt_cross, width = 8, height = 8, dpi = 300)
}
####### cross layer slope  ########
for (param in c("cross_rho")){ 
   s[grep(glue::glue("^{param}"), rownames(s)),c(1,4,8,9,10)]
  df1 <- s[grep(glue::glue("^{param}_coef"), rownames(s)),c(1,4,8,9,10)] %>% data.frame()
  df1$class <- class_size_id
  df1$layer <- "bully1 x bully1_rev"
  df <- df1
  df$type = "random_slope"
  colnames(df) <- c("mean", "lower", "upper", "n_eff", "Rhat", "class", "layer", "type")

  overall_df <- s[grep(glue::glue("^{param}_fixed"), rownames(s)),c(1,4,8,9,10)]

  # s[grep(glue::glue("^{param}\\[1\\]"), rownames(s)),c(1,4,8,9,10)]

  # s[grep(glue::glue("^{param}_random\\["), rownames(s)),c(1,4,8,9,10)] %>% data.frame()

  overall_df <- as.data.frame(t(overall_df)) %>% mutate(
    class = "Overall",
    layer = "bully1 x bully1_rev",
    type = "overall_slope"
  )
  colnames(overall_df) <- c("mean", "lower", "upper", "n_eff", "Rhat", "class", "layer", "type")
  df <-  rbind(df, overall_df)

  # Intercepts
  slp_plt_cross <- df %>%
    mutate(class  = as.factor(class)) %>%
    ggplot(aes(x = mean, y = reorder_within(class, mean, layer), fill = type, color = type)) +
    geom_errorbarh(aes(xmin = lower,
                      xmax = upper, 
                  color = type),
                  height = 0.2) +
    geom_point(size = 3,
              shape = 21,
              color = "black") +
    geom_vline(xintercept = 0,
              color = "red",
              size = 0.8,
              linetype = "dashed") +
  #   scale_x_continuous(breaks = seq(-60, 60, 20), limits = c(-20, 10)) +
    labs(x = "Slope Estimate",
        y = "Classroom",
        title = glue::glue("same ses on: {param}")) +
    custom_theme +
    facet_wrap(~ layer, ncol = 1, scales = "free_y") +  # Allow free y-axis scales
    scale_y_reordered()
  slp_plt_cross

  ggsave(glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/plot/{param}_slope_plot_{file_name}.png"), plot = slp_plt_cross, width = 8, height = 8, dpi = 300)
}

##### group covariates ######
s[grep(glue::glue("group_coef"), rownames(s)),c(1,4,8,9,10)] 
s[grep(glue::glue("mu_fixed"), rownames(s)),c(1,4,8,9,10)]
