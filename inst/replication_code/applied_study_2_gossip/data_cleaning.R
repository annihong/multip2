library(dplyr)
library(igraph)

#######helper functions#######
# Replace missing data codes with NA
recode_missing <- function(dat, old_code = 10){
  dat[dat == old_code] <- NA
  return(dat)
}

# Check if both rows and columns for a student are entirely missing
check_student_missing <- function(network){
  row <- row_missing(network) 
  col <- col_missing(network) 
  return(row & col)
}

# Extract the last n characters from a string
extract_last_n_char <- function(x, n){
  substr(x, nchar(x) - n + 1, nchar(x))
}

# Compare two networks and calculate similarity (Jaccard index)
graph_similarity <- function(matrix1, matrix2){
  matrix1 <- binary_network(matrix1)
  matrix2 <- binary_network(matrix2)
  both <- sum((matrix1 == matrix2) & (matrix1 == 1), na.rm = T) 
  total_1 <- sum(matrix1 == 1, na.rm = T)
  only_1 <- sum((matrix1 != matrix2) & (matrix1 == 1), na.rm = T)
  total_2 <- sum(matrix2 == 1, na.rm = T)
  only_2 <- sum((matrix1 != matrix2) & (matrix2 == 1), na.rm = T)
  total <- (both + only_1 + only_2)
  jaccard <- both / total
  
  return(list(both=both, total_1=total_1, only_1=only_1, total_2=total_2, only_2=only_2, total=total, jaccard=jaccard))
}

# Calculate the density of ties in the network
get_density <- function(dat, tie = 1){
  n = nrow(dat)
  ties <- sum(dat == tie, na.rm=T)
  d <- n * (n -1)
  return(ties/d)
}
#####################################

# Define constants and read in data
PATH <- "./Data/Raw/GossipBully" #path to the raw data
classes <- list.files(PATH)
class_ids <- regmatches(classes, regexpr("[0-9]+", classes))
SUMMARY_PATH <- "./Data/Raw/documentation/wellbeing_gossip20200120.csv" #path to the summary data
summary <- read.csv(SUMMARY_PATH)

# Define a function to extract network summary statistics for a given class and wave
extract_network_summary <- function(class_id, wave, summary) {
  # Load network data for the given class and wave
  load(file.path(PATH, paste0(class_id, ".RData")))
  # Keep only individuals from the given class
  keep <- grep(paste0("^", substr(class_id, start = 1, stop = 3)), rownames(gossip1[[wave]]))
  sts_raw <- recode_missing(gossip1[[wave]])#[keep, keep]
  tst_raw <- recode_missing(gossip1_rev[[wave]])#[keep, keep]
  # Count the number of individuals not in the given class
  num_not_in_class <- length(rownames(gossip1[[wave]])[-keep])
  # Filter the summary data to only include individuals from the given classv (just to see, actual data does not exclude students not in the class)
  filtered_summary <- summary %>%
    #filter(class == class_id) %>%
    select(idcode, gender) %>%
    distinct()
  # Merge the network data with the summary data to get gender information
  merged <- left_join(data.frame(idcode = as.numeric(rownames(sts_raw))), filtered_summary, by = "idcode")
  is_female <- merged$gender == 2
  names(is_female) <- as.character(merged$idcode)
  # Remove absent students from the networks
  networks <- list(sts = sts_raw, tst = tst_raw)
  all_missing <- data.frame(lapply(networks, check_student_missing))
  both_missing <- apply(all_missing, 1, function(x) all(x))
  absent_students <- rownames(all_missing)[both_missing]
  absent_students <- sapply(absent_students, extract_last_n_char, n = 2)
  num_all_missing <- length(absent_students)
  sts <- sts_raw[!both_missing, !both_missing]
  tst <- tst_raw[!both_missing, !both_missing]
  cleaned_data <- list(sts=sts, tst=tst, is_female=is_female)
  save(cleaned_data, file = file.path(PATH, paste0(class_id,"_wave_",wave,"_cleaned.RData")))
  # Compute summary statistics
  sts_igraph <- igraph::graph_from_adjacency_matrix(sts)
  tst_igraph <- igraph::graph_from_adjacency_matrix(tst)
  no_missing <- !sum(is.na(sts)) & !sum(is.na(sts))
  sender_psyche_sim <- graph_similarity(sts, tst)
  target_percep_sim <- graph_similarity(sts, t(tst))
  prop_female <- mean(is_female[rownames(sts)], na.rm = T)
  sts_density <- get_density(sts)
  tst_density <- get_density(tst)
  # Create a list with the computed summary statistics
  res <- list(
    class_id = class_id,
    wave = wave,
    num_removed = num_not_in_class,
    num_absent = num_all_missing,
    sts_density = sts_density,
    sts_reciprocity = igraph::reciprocity(sts_igraph),
    tst_density = tst_density,
    tst_reciprocity = igraph::reciprocity(tst_igraph),
    sender_psyche_jaccard = sender_psyche_sim$jaccard,
    target_percep_jaccard = target_percep_sim$jaccard,
    num_present = nrow(sts),
    prop_female = prop_female,
    missing_gender = mean(is.na(is_female[rownames(sts)]))
  )
  return(res)
}

# loop through all the classes and waves
results_df <- data.frame()
for (class_id in class_ids) {
  for (wave in 4) {
    # call the extract_network_summary function for each class and wave
    info <- extract_network_summary(class_id, wave, summary)
    # add the resulting information to the results data frame
    results_df <- bind_rows(results_df, info)
  }
}


filtered_df <- results_df %>% 
  filter(wave == 4, missing_gender <= 0, tst_density >= 0, sts_density >= 0) %>%
  unique()

save(filtered_df, file = "./filtered_df_no_missing_wave_4.RData")


