#############################################################
####                   Load libraries                    ####
#############################################################
library(data.table)
library(dplyr)
library(keras)
library(caret)
library(pROC)

#############################################################
####                     Load data                       ####
#############################################################
setwd("./result")
dat <- readRDS("CD4dat.rds")
metadata <- dat@meta.data
metadata$samp <- rownames(metadata)
dat <- as.data.frame(t(as.matrix(dat@assays$RNA@scale.data)))

#############################################################
####                     DE genes                        ####
#############################################################
output <- fread("CD4_SeverevsHC.txt", sep="\t")
output1 <- na.omit(output)
dat <- dat[, which(colnames(dat) %in% output1$gene)]

#############################################################
####                    Severe vs Healthy                ####
#############################################################
group_index <- which(metadata$group == "Yes" | metadata$group == "Healthy")
#### test ####
healthy_index <- which(metadata$group == "Healthy")
covid_index <- which(metadata$group == "Yes")
set.seed(123)
healthy_index_test <- sample(healthy_index, round(length(healthy_index) / 5, 0), replace = F)
covid_index_test <- sample(covid_index, round(length(covid_index) / 5, 0), replace = F)
dat_test <- dat[c(healthy_index_test, covid_index_test), ]
metadata_test <- metadata[c(healthy_index_test, covid_index_test), ]
metadata_test$group_id <- c(rep(0, length(healthy_index_test)), rep(1, length(covid_index_test)))

healthy_index_training <- healthy_index[!(healthy_index %in% healthy_index_test)]
covid_index_training <- covid_index[!(covid_index %in% covid_index_test)]
dat_training <- dat[c(healthy_index_training, covid_index_training), ]
metadata_training <- metadata[c(healthy_index_training, covid_index_training), ]
metadata_training$group_id <- c(rep(0, length(healthy_index_training)), rep(1, length(covid_index_training)))

rm(dat, output, output1)


metadata <- rbind(metadata_training, metadata_test)
dat <- rbind(dat_training, dat_test)
id_group <- unique(metadata[,c("cname", "group_id")])

#### Final results ####
results_noweight_adam <- results_noweight_rmsprop <- results_weight_adam <- results_weight_rmsprop <- as.data.frame(matrix(NA, nrow = nrow(id_group), ncol = 6))
colnames(results_noweight_adam) <- colnames(results_noweight_rmsprop) <- colnames(results_weight_adam) <- colnames(results_weight_rmsprop) <- 
    c("ID", "group", "n_training", "n_test", "n_predict", "accuracy")

results_CV <- as.data.frame(matrix(NA, nrow = nrow(id_group), ncol = 6))
colnames(results_CV) <- c("ID", "group", "noweight_adam", "noweight_rmsprop", "weight_adam", "weight_rmsprop")

#### MLP ####
for (i in 1:nrow(id_group)) {
    ## data pre-process
    test_index <- which(metadata$cname == id_group[i, "cname"])
    dat_training <- dat[-test_index,]
    dat_test <- dat[test_index,]
    metadata_training <- metadata[-test_index,]
    metadata_test <- metadata[test_index,]
    
    ## data shuffle
    temp_train <- dat_training
    temp_meta_train <- metadata_training
    set.seed(i)
    shuffle_index <- sample(1:nrow(temp_train), nrow(temp_train), replace = F)
    dat_training <- temp_train[shuffle_index, ]
    metadata_training <- temp_meta_train[shuffle_index, ]
    
    results_noweight_adam[i, c("ID", "group")] <- results_noweight_rmsprop[i, c("ID", "group")] <- 
        results_weight_adam[i, c("ID", "group")] <- results_weight_rmsprop[i, c("ID", "group")] <- id_group[i, ]
    results_noweight_adam[i, "n_training"] <- results_noweight_rmsprop[i, "n_training"] <- results_weight_adam[i, "n_training"] <- results_weight_rmsprop[i, "n_training"] <- nrow(dat_training)
    results_noweight_adam[i, "n_test"] <- results_noweight_rmsprop[i, "n_test"] <- results_weight_adam[i, "n_test"] <- results_weight_rmsprop[i, "n_test"] <- nrow(dat_test)
    
    results_CV[i, c("ID", "group")] <- id_group[i, ]
    
    
    ## 1: rmsprop, with weight
    model <- keras_model_sequential()
    model %>%
        layer_dense(units = 128, activation = 'relu', input_shape = c(ncol(dat_training))) %>%
        layer_dense(units = 128, activation = 'relu') %>%
        layer_dropout(rate=0.5) %>%
        layer_dense(units = 128, activation = 'relu') %>%
        layer_dropout(rate=0.5) %>%
        layer_dense(units = 2, activation = 'sigmoid')
    
    # Compile the model
    model %>% compile(
        loss = 'binary_crossentropy',
        optimizer = 'rmsprop',
        metrics = c('accuracy')
    )
    
    # Store the fitting history in `history`
    history <- model %>% fit(
        as.matrix(dat_training),
        to_categorical(metadata_training$group_id),
        epochs = 100,
        batch_size = 32,
        verbose = 0,
        class_weight = list("0" = table(metadata_training$group_id)[2]/table(metadata_training$group_id)[1], "1" = 1),
        validation_split = 0.2
    )
    results_CV[i, "weight_rmsprop"] <- history$metrics$val_accuracy[100]
    classes <- model %>% predict(as.matrix(dat_test)) %>% `>`(0.5) %>% k_cast("int32")
    results_weight_rmsprop[i, "n_predict"] <- length(which(as.matrix(classes)[, 2] == id_group[i, "group_id"]))
    results_weight_rmsprop[i, "accuracy"] <- results_weight_rmsprop[i, "n_predict"]/results_weight_rmsprop[i, "n_test"]
    
    
    ## 2: rmsprop, without weight
    model <- keras_model_sequential()
    model %>%
        layer_dense(units = 128, activation = 'relu', input_shape = c(ncol(dat_training))) %>%
        layer_dense(units = 128, activation = 'relu') %>%
        layer_dropout(rate=0.5) %>%
        layer_dense(units = 128, activation = 'relu') %>%
        layer_dropout(rate=0.5) %>%
        layer_dense(units = 2, activation = 'sigmoid')
    
    # Compile the model
    model %>% compile(
        loss = 'binary_crossentropy',
        optimizer = 'rmsprop',
        metrics = c('accuracy')
    )
    
    # Store the fitting history in `history`
    history <- model %>% fit(
        as.matrix(dat_training),
        to_categorical(metadata_training$group_id),
        epochs = 100,
        batch_size = 32,
        verbose = 0,
        #class_weight = list("0" = table(metadata_training$group_id)[2]/table(metadata_training$group_id)[1], "1" = 1),
        validation_split = 0.2
    )
    results_CV[i, "noweight_rmsprop"] <- history$metrics$val_accuracy[100]
    classes <- model %>% predict(as.matrix(dat_test)) %>% `>`(0.5) %>% k_cast("int32")
    results_noweight_rmsprop[i, "n_predict"] <- length(which(as.matrix(classes)[, 2] == id_group[i, "group_id"]))
    results_noweight_rmsprop[i, "accuracy"] <- results_noweight_rmsprop[i, "n_predict"]/results_noweight_rmsprop[i, "n_test"]
    
    
    ## 3: adam, weight
    model <- keras_model_sequential()
    model %>%
        layer_dense(units = 128, activation = 'relu', input_shape = c(ncol(dat_training))) %>%
        layer_dense(units = 128, activation = 'relu') %>%
        layer_dropout(rate=0.5) %>%
        layer_dense(units = 128, activation = 'relu') %>%
        layer_dropout(rate=0.5) %>%
        layer_dense(units = 2, activation = 'sigmoid')
    
    # Compile the model
    model %>% compile(
        loss = 'binary_crossentropy',
        optimizer = 'adam',
        metrics = c('accuracy')
    )
    
    # Store the fitting history in `history`
    history <- model %>% fit(
        as.matrix(dat_training),
        to_categorical(metadata_training$group_id),
        epochs = 100,
        batch_size = 32,
        verbose = 0,
        class_weight = list("0" = table(metadata_training$group_id)[2]/table(metadata_training$group_id)[1], "1" = 1),
        validation_split = 0.2
    )
    results_CV[i, "weight_adam"] <- history$metrics$val_accuracy[100]
    classes <- model %>% predict(as.matrix(dat_test)) %>% `>`(0.5) %>% k_cast("int32")
    results_weight_adam[i, "n_predict"] <- length(which(as.matrix(classes)[, 2] == id_group[i, "group_id"]))
    results_weight_adam[i, "accuracy"] <- results_weight_adam[i, "n_predict"]/results_weight_adam[i, "n_test"]
    
    
    ## 4: adam, without weight
    model <- keras_model_sequential()
    model %>%
        layer_dense(units = 128, activation = 'relu', input_shape = c(ncol(dat_training))) %>%
        layer_dense(units = 128, activation = 'relu') %>%
        layer_dropout(rate=0.5) %>%
        layer_dense(units = 128, activation = 'relu') %>%
        layer_dropout(rate=0.5) %>%
        layer_dense(units = 2, activation = 'sigmoid')
    
    # Compile the model
    model %>% compile(
        loss = 'binary_crossentropy',
        optimizer = 'adam',
        metrics = c('accuracy')
    )
    
    # Store the fitting history in `history`
    history <- model %>% fit(
        as.matrix(dat_training),
        to_categorical(metadata_training$group_id),
        epochs = 100,
        batch_size = 32,
        verbose = 0,
        #class_weight = list("0" = table(metadata_training$group_id)[2]/table(metadata_training$group_id)[1], "1" = 1),
        validation_split = 0.2
    )
    results_CV[i, "noweight_adam"] <- history$metrics$val_accuracy[100]
    classes <- model %>% predict(as.matrix(dat_test)) %>% `>`(0.5) %>% k_cast("int32")
    results_noweight_adam[i, "n_predict"] <- length(which(as.matrix(classes)[, 2] == id_group[i, "group_id"]))
    results_noweight_adam[i, "accuracy"] <- results_noweight_adam[i, "n_predict"]/results_noweight_adam[i, "n_test"]
}


save(list = c('results_CV',
              'results_noweight_adam',
              'results_noweight_rmsprop',
              'results_weight_adam',
              "results_weight_rmsprop"
),
file = "CD4_MLP_severe_results_id.RData")