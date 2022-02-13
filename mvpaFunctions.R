########################
## Functions for MVPA ##
########################

###################################################################

## Loads and combines behavioral and fMRI data 

getData <- function(ssid, cond, roi) {
  
  # load behavioral data
  fname <- sprintf("beh/sub-%s_task-study_items.tsv", ssid)
  beh <- rio::import(fname) %>% 
    mutate(design = factor(design, labels = c("med8","quick4n","quick6",
                                              "slow10","slow12")),
           category = factor(category, labels = c("animal","tool")))
  
  # load betas for each voxel
  ## exclude first 3 rows, just coordinates of voxel
  fname <- sprintf("betas/item_betas_sub_%s_%s_b_%s.txt", ssid, cond, roi)
  betas <- rio::import(fname) %>%
    slice(-1:-3)
  
  # Combine data for classification
  df <- cbind(beh, betas) %>%
    filter(design == cond)
  
  return(df)
}

###################################################################

## Split the data into folds
# outputs list with train/test splits
# each element in list contains 1 train set and 1 test set
cvfold <- function(df, fold_by, dv) {
  
  # # make these executable in the code below
  #fold_by <- eval(substitute(fold_by), df, parent.frame())
  
  # number of folds should equal number of levels
  num_folds <- length(unique(df[,fold_by]))
  
  # create train/test set for each fold
  folds <- list()
  for (i in c(1:num_folds)) {
    train <- df %>%
      filter_at(vars(contains(fold_by)), ~ .x != i) %>% 
      select(dv, contains("V"))
    test <- df %>%
      filter_at(vars(contains(fold_by)), ~ .x == i) %>% 
      select(dv, contains("V"))
    folds[[i]] <- list(train = train,
                       test = test)
  }
  
  return(folds)
}

###################################################################

## Subset the data using given ratios
# iterate over X number of times for more stable estimates
# will really only work with 2 groups
# outputs all possible training sets

subSample <- function(train, ratio, times, ngroups = 2) {
  
  # get number of observations per class
  class_size <- nrow(train)/ngroups
  
  # if the class sizes are balanced
  if (ratio[1] == ratio[2]) {
    
    # determine size of each class
    class_size <- ratio[1] * (nrow(train)/ngroups)
    
    # down sample each class
    trainsub <- train %>% 
      group_by_at(1) %>% 
      slice_sample(n = class_size) %>% 
      ungroup()
    
    trainsub <- list(trainsub)
    
    # if the class sizes are not balanced  
  } else {
    
    # determine size of the smaller class
    subsize <- ratio[1] * (nrow(train)/ngroups)
    
    # downsample the smaller class
    trainsub <- list()
    for (i in 1:times) {
      
      # split training data
      trainsplit <- train %>% 
        group_by_at(1) %>% 
        group_split()
      
      # downsample group 1
      sampGroup1 <- trainsplit[[1]] %>% 
        slice_sample(n = subsize) 
      
      # downsample group 2
      sampGroup2 <- trainsplit[[2]] %>% 
        slice_sample(n = subsize) 
      
      trainsub[[i]] <- list(rbind(sampGroup1, trainsplit[[2]]),
                            rbind(trainsplit[[1]], sampGroup2))
    }
    
    trainsub <- flatten(trainsub)
  }
  
  return(trainsub)
}

###################################################################

## Run ANOVA feature selection
# input train set where column 1 = outcome
# outputs index for selected voxels

featSelect <- function(train, nfeatures) {
  
  train_split <- train %>% 
    group_by_at(1) %>% 
    group_split()
  
  fvals <- vector()
  for (fs in 2:ncol(train)) {
    tmp <- t.test(train_split[[1]][,fs],
                  train_split[[2]][,fs])
    fvals[fs-1] <- tmp$statistic[[1]]^2
  }
  
  top_features <- sort(fvals, decreasing = TRUE, index.return = TRUE)
  idx <- top_features$ix[1:nfeatures] + 1 # need to shift indices up 1 to match columns
}

###################################################################

## Run SVM classification

runSVM <- function(train, test) {
  
  # # troubleshoot
  # train <- folds[[1]][["train"]]
  # test <- folds[[1]][["test"]]
  
  # train model
  x <- train[,-1]
  y <- train[,1]
  train_mdl <- e1071::svm(x, y,
                          scale = TRUE,
                          type = "C-classification",
                          kernel = "linear")
  
  # predict test set 
  pred <- predict(train_mdl, test[,-1], decision.values = FALSE)
  #decvals <- attr(pred, "decision.values")
  y <- test[,1]
  
  # evaluate performance
  results <- caret::confusionMatrix(pred, y)
  
  test_metrics <- tibble(acc = results$overall[["Accuracy"]],
                         hit_class1 = results$byClass[["Sensitivity"]],
                         hit_class2 = results$byClass[["Specificity"]],
                         fa_class1 = results$table[1,2]/sum(results$table[,2]),
                         fa_class2 = results$table[2,1]/sum(results$table[,2]))
  
  return(test_metrics)
}
