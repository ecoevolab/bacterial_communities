# Author: Mayra Beatriz Mendoza Velazquez
# Title: Exploration of reconstruction methods

library(readr)
library(mlBioNets)

# Load csv archives 
metadata <- as.data.frame(read_tsv(file = "01_RawData/metadata_clean.tsv"))
f_clean <- as.data.frame(read_tsv("01_RawData/f_clean.tsv"))
rzcompositiondata <- read.csv("01_RawData/rzcomposition.csv")

#### Function ####
# Function for selecting specific columns based on the community sample 
community_isolation <- function(rcommunities, sample, sample_col, rcommunities_col, df, arrangev, interest_column, dfwvals){

   library(dplyr)
# create an empty list   
community_list <- list()
k <- 1

# subset the commuinity values based on the day and community name 

for (i in 1:length(rcommunities)){
  for (x in 1:length(sample)){
    
    community_list[[k]] <- subset(df, sample_col %in% c(0, sample[x])  & rcommunities_col == rcommunities[i]) %>% 
      arrange(arrangev) %>% 
      pull(interest_column) %>% 
      as.character()
    
    k <- k + 1
  }
}

# create an empty list for the abundances 
abundances_tables <- list()

# Sbased on the column names, select the values from f_clean 
#
for (id in seq_along(community_list)) {
  abundances_tables[[id]] <- dfwvals %>%
    dplyr::select(1 | all_of(community_list[[id]])) %>%
    column_to_rownames(dfwvas[,1])
  
}

return(abundances_tables)

}

## Function testing 
# Change the NA values in the data frame (for a correct functioning of the subsetting function)
metadata[is.na(metadata)] <- 0

comms_rhiz <- unique(metadata$community) # to use all the community values for subsetting the df 
temps_rhiz <- c(28, 32) # without the 0 temperature, because it is already being used in the function 

# Using function 
rz_communities <- community_isolation(rcommunities = comms_rhiz , sample = temps_rhiz , sample_col = "temp", 
                                      rcommunities_col = "community", df = metadata, arrangev = "day", interest_column = "label", dfwvals = f_clean)


## FUNCTION DOESN'T WORK ## 

# isolated "for loops" to subset the data.frames #

community_list <- list()
k <- 1
# subset the commuinity values based on the day and community name 

for (i in 1:length(comms_rhiz)){
  for (x in 1:length(temps_rhiz)){
    
    community_list[[k]] <- subset(metadata, temp %in% c(0, temps_rhiz[x])  & community == comms_rhiz[i]) %>% 
      arrange(day) %>% 
      pull(label) %>% 
      as.character()
    
    k <- k + 1
  }
}

# create an empty list for the abundances 
abundances_tables <- list()

# Based on the column names, select the values from f_clean 
for (id in seq_along(community_list)) {
  abundances_tables[[id]] <- f_clean %>%
    dplyr::select(1 | all_of(community_list[[id]]) )  %>%
    column_to_rownames(var = "row.names")
}

# filter the list of data.frames based on the community composition 
abnds_filtered <- list()

for (g in seq_len(ncol(rzcompositiondata))) {
  
  # select the community name / for arranging the list 
  comm_name <- colnames(rzcompositiondata)[g]
  
  # select the specific bacterial id for the selected community 
  spp <- rzcompositiondata[, g]
  spln <- 1:length(spp)
  # indexes [to select in the list of data.frames]
  idx <- (2*g - 1):(2*g)
  
  for (j in seq_along(idx)) {
    
    k <- idx[j]
    
    abnds_filtered[[paste0(comm_name, "_", temps_rhiz[j])]] <- abundances_tables[[k]] %>% 
      filter(row.names(abundances_tables[[k]]) %in% spp)

    
  }
}

abnds_filtered

##### Testing network inference algortihms ####

# aracne - exploring the method and specifically the network generated for R11 (we only have one sample)
R1_0_inference <- net_inference(taxa_abs = t(abundances_tables[[1]]), method = "aracne")
plot(R1_0_inference)
R11_3_inference <- net_inference(taxa_abs = t(abundances_tables[[44]]), method = "aracne")               
plot(R11_3_inference)
