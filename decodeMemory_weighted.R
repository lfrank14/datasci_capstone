########################################
## Decoding Memory for Visual Stimuli ##
########################################

## Setup 
#=========
library(tidyverse)
# install.packages("e1071")
# install_github("lfrank14/lurr")
source("mvpaFunctions.R")


## Define Conditions
#====================
# specify subject ids
ssids <- str_pad(seq(1,35), width = 2, 
                 side = "left", pad = "0")

# specify trial timing conditions
timings <- c("quick4n")

# specify regions of interest (ROIs)
rois <- c("pfus")


## Create Data Frame for Conditions
#===================================
conds <- expand.grid(ssid = ssids,
                     timing = timings,
                     roi = rois) %>% 
  # not sure why creating factors
  mutate_if(is.factor, as.character)


## Add Data 
#===========
conds <- conds %>% 
  rowwise() %>% 
  mutate(df = list(getData(ssid, timing, roi)))


## Create Folds 
#===============
conds <- conds %>% 
  rowwise() %>% 
  mutate(folds = list(cvfold(df,
                             fold_by = "run",
                             dv = "remembered")))


## Run Cross Validation 
#=======================
tictoc::tic()
conds <- conds %>% 
  rowwise() %>% 
  mutate(output = list(runCV(folds = folds,
                             weighted = TRUE)),
         nremembered = sum(df$remembered==1)/2,
         nforgotten = sum(df$remembered==0)/2) %>% 
  suppressWarnings()
tictoc::toc()
# 553.2 s


## Save All Results 
#===================
saveRDS(conds,
        "models/decodeWMem.rds")


## Format and Save Metrics Only 
#===============================
results <- conds %>% 
  select(ssid, timing, roi, nremembered, nforgotten, output) %>%
  unnest(output)

saveRDS(results,
        "models/decodeWMem_results.rds")

