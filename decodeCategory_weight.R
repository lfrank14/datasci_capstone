#########################################
## Decoding Category of Visual Stimuli ##
#########################################

## Setup ##
library(tidyverse)
# install.packages("e1071")
# install_github("lfrank14/lurr")
source("mvpaFunctions.R")


## Define Conditions ##
# specify subject ids
ssids <- str_pad(seq(1,35), width = 2, 
                 side = "left", pad = "0")

# specify trial timing conditions
timings <- c("quick4n")

# specify regions of interest (ROIs)
rois <- c("pfus")

# specify ratios for different class sizes
ratios <- c("c(1,1)",
            "c(.5,.5)",
            "c(2/3,1)",
            "c(1/2,1)",
            "c(1/3,1)")


## Create Data Frame for Conditions ##
conds <- expand.grid(ssid = ssids,
                     timing = timings,
                     roi = rois,
                     ratio = ratios) %>% 
  # not sure why creating factors
  mutate_if(is.factor, as.character)


## Add Data ##
conds <- conds %>% 
  rowwise() %>% 
  mutate(df = list(getData(ssid, timing, roi)))


## Create Folds ##
conds <- conds %>% 
  rowwise() %>% 
  mutate(folds = list(cvfold(df,
                             fold_by = "run",
                             dv = "category")))


## Run Cross Validation ##
tictoc::tic()
conds <- conds %>% 
  rowwise() %>% 
  mutate(output = list(runCV(folds,
                             eval(parse(text = ratio)),
                             nsamples = 5,
                             nfeats = 100,
                             weighted = TRUE)))
tictoc::toc()
# 2564.987 s


## Save All Results ##
saveRDS(conds,
        "models/decodeWCat.rds")


## Format and Save Metrics Only ##
results <- conds %>% 
  ungroup() %>% 
  select(ssid, timing, roi, ratio, output) %>% 
  mutate(ratio = fct_recode(ratio,
                            "1:1" = "c(1,1)",
                            "2:3" = "c(2/3,1)",
                            ".5:.5" = "c(.5,.5)",
                            "1:2" = "c(1/2,1)",
                            "1:3" = "c(1/3,1)")) %>% 
  unnest(output)

saveRDS(results,
        "models/decodeWCat_results.rds")


