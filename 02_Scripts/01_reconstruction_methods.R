# Author: Mayra Beatriz Mendoza Velazquez
# Title: Exploration of reconstruction methods

library(readr)
# 1. Subset communities 

metadata <- as.data.frame(read_tsv(file = "01_RawData/metadata_clean.tsv"))
f_clean <- as.data.frame(read_tsv("01_RawData/f_clean.tsv"))


rzcomunidades <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12")
dias <- 0:3

community_isolation <- function(rcommunities, sample, sample_col, rcommunities_col, df, abnd_tbl, interest_column, dfwvals, ){

   library(dplyr)
# create an empty list   
community_list <- list()
k <- 1

# subset the commuinity values based on the day and community name 

for (i in 1:length(rcommunities)){
  for (x in 1:length(sample)){
    
    community_list[[k]] <- subset(df, sample_col == sample[x] & rcommunities_col == rcommunities[i]) %>% 
    pull(label) %>% 
    as.character()
    
    k <- k + 1
  }
}

# create an empty list for the abundances 
abundances_tables <- list()

# select, based on the column names, the values from f_clean 
#
for (id in seq_along(community_list)) {
  abundances_tables[[id]] <- dfwvals %>%
    dplyr::select(all_of(community_list[[id]]))
}

return(abundances_tables)
# Note: The data.frame (df), must be pre-processed by one of the variable of interest 
}
