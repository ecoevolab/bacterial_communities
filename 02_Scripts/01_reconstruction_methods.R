# Author: Mayra Beatriz Mendoza Velazquez
# Title: Exploration of reconstruction methods

library(readr)
library(mlBioNets)

# 1. Subset communities 

metadata <- as.data.frame(read_tsv(file = "01_RawData/metadata_clean.tsv"))
f_clean <- as.data.frame(read_tsv("01_RawData/f_clean.tsv"))

# Function for selecting specific columns based on the community sample 
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

# Sbased on the column names, select the values from f_clean 
#
for (id in seq_along(community_list)) {
  abundances_tables[[id]] <- dfwvals %>%
    dplyr::select(all_of(community_list[[id]]))
}

return(abundances_tables)

}

## Function testing 
rz_communities <- community_isolation(rcommunities =  , sample, sample_col, rcommunities_col, df, abnd_tbl, interest_column, dfwvals)

abundances_tables
##### Testing network inference algortihms ####

# aracne - exploring the method and specifically the network generated for R11 (we only have one sample)
R1_0_inference <- net_inference(taxa_abs = t(abundances_tables[[1]]), method = "aracne")
plot(R1_0_inference)
R11_3_inference <- net_inference(taxa_abs = t(abundances_tables[[44]]), method = "aracne")               
plot(R11_3_inference)

rz_comunidades <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12")
temperatura <- c(28, 32)

community_list <- list()
k <- 1


for (i in 1:length(rz_comunidades)){
  for (x in 1:length(temperatura)){
    
    community_list[[k]] <- subset(metadata, temp == temperatura[x] & community == rz_comunidades[i]) %>% 
      pull(label) %>% 
      as.character()
    
    k <- k + 1
  }
}
