library(dplyr)
library(igraph)

#######helper functions#######
# Replace missing data codes with NA
recode_missing <- function(dat, old_code = 10){
  dat[dat == old_code] <- NA
  return(dat)
}

binary_network <- function(network){
  return(apply(network, 1:2, function(x) as.numeric(x > 0)))
}

# Check if both rows and columns for a student are entirely missing
check_student_missing <- function(network){
  row <- row_missing(network) 
  col <- col_missing(network) 
  return(row & col)
}

#given a network, returns a logical vector on whether the entire row is missing
row_missing <- function(network){
  apply(network, 1, function(x) all(is.na(x)))
}

col_missing <- function(network){
  apply(network, 2, function(x) all(is.na(x)))
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
PATH <- "/home/annihong/projects/multiplex-social-networks/stan/coding/gossip_perception/Data/Raw/GossipBully" #path to the raw data
classes <- list.files(PATH)
class_ids <- regmatches(classes, regexpr("[0-9]+", classes))
SUMMARY_PATH <- "/home/annihong/projects/multip2/data-raw/gossip/summary.csv" #path to the summary data
summary <- read.csv(SUMMARY_PATH)

# Define a function to extract network summary statistics for a given class and wave
extract_network_summary <- function(class_id, wave, summary) {
  # Load network data for the given class and wave
  load(file.path(PATH, paste0(class_id, "Nynke.RData")))
  # Keep only individuals from the given class
  keep <- grep(paste0("^", substr(class_id, start = 1, stop = 3)), rownames(bully1[[wave]]))
  bully1_raw <- recode_missing(bully1[[wave]])#[keep, keep]
  bully1_rev_raw <- recode_missing(bully1_rev[[wave]])#[keep, keep]
  # Count the number of individuals not in the given class
  num_not_in_class <- length(rownames(bully1[[wave]])[-keep])
  # # Filter the summary data to only include individuals from the given classv (just to see, actual data does not exclude students not in the class)
  filtered_summary <- summary %>%
    #filter(class == class_id) %>%
    select(idcode, gender) %>%
    distinct()
  # Merge the network data with the summary data to get gender information

  # Remove absent students from the networks
  networks <- list(bully1 = bully1_raw, bully1_rev = bully1_rev_raw)
  all_missing <- data.frame(lapply(networks, check_student_missing))
  both_missing <- apply(all_missing, 1, function(x) all(x))
  absent_students <- rownames(all_missing)[both_missing]
  absent_students <- sapply(absent_students, extract_last_n_char, n = 2)
  num_all_missing <- length(absent_students)
  bully1 <- bully1_raw[!both_missing, !both_missing]
  bully1_rev <- bully1_rev_raw[!both_missing, !both_missing]

  merged <- left_join(data.frame(idcode = as.numeric(rownames(bully1))), filtered_summary, by = "idcode")
  is_female <- merged$gender == 2
  names(is_female) <- as.character(merged$idcode)

  cleaned_data <- list(bully1=bully1, bully1_rev=bully1_rev, is_female=is_female)
  save(cleaned_data, file = file.path(PATH, paste0(class_id,"_bully1_wave_",wave,"_cleaned.RData")))
  # Compute summary statistics
  bully1_igraph <- igraph::graph_from_adjacency_matrix(bully1)
  bully1_rev_igraph <- igraph::graph_from_adjacency_matrix(bully1_rev)
  no_missing <- !sum(is.na(bully1)) & !sum(is.na(bully1))
  sender_psyche_sim <- graph_similarity(bully1, bully1_rev)
  target_percep_sim <- graph_similarity(bully1, t(bully1_rev))
  prop_female <- mean(is_female[rownames(bully1)], na.rm = TRUE)
  bully1_density <- get_density(bully1)
  bully1_rev_density <- get_density(bully1_rev)
  # Create a list with the computed summary statistics
  res <- list(
    class_id = class_id,
    wave = wave,
    num_removed = num_not_in_class,
    num_absent = num_all_missing,
    bully1_density = bully1_density,
    bully1_reciprocity = igraph::reciprocity(bully1_igraph),
    bully1_rev_density = bully1_rev_density,
    bully1_rev_reciprocity = igraph::reciprocity(bully1_rev_igraph),
    sender_psyche_jaccard = sender_psyche_sim$jaccard,
    target_percep_jaccard = target_percep_sim$jaccard,
    num_present = nrow(bully1),
    prop_female = prop_female,
    missing_gender = mean(is.na(is_female[rownames(bully1)]))
  )
  return(res)
}

# loop through all the classes and waves
results_df <- data.frame()
for (class_id in class_ids) {
  for (wave in 1:6) {
    # call the extract_network_summary function for each class and wave
    info <- extract_network_summary(class_id, wave, summary)
    # add the resulting information to the results data frame
    results_df <- bind_rows(results_df, info)
  }
}

results_df %>% filter(wave == 4) %>% unique() %>% dim() # 37 classrooms 
results_df %>% filter(wave == 4, missing_gender == 0) %>% unique() %>% dim() # 34 classrooms, 3 classrooms excluded due to missing gender.  


filtered_df <- results_df %>% 
  filter(wave == 4, missing_gender == 0, bully1_rev_density > 0, bully1_density > 0) %>%
  unique()


write.csv(filtered_df, file = "/home/annihong/projects/multip2/data-raw/bully/filtered_df_no_missing_wave_4_bully1.csv")


