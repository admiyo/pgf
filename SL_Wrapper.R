# GENERALIZED ADDITIVE MODELS
SL.gam.1 <- function(..., deg.gam = 1) { 
  SL.gam(..., deg.gam = deg.gam)
}
SL.gam.2 <- function(..., deg.gam = 2) { 
  SL.gam(..., deg.gam = deg.gam)
}
SL.gam.3 <- function(..., deg.gam = 3) { 
  SL.gam(..., deg.gam = deg.gam)
}
SL.gam.4 <- function(..., deg.gam = 4) { 
  SL.gam(..., deg.gam = deg.gam)
}
# GENERALIZED BOOSTED MODELS
SL.gbm.1 <- function(...) { 
  SL.gbm(..., interaction.depth = 1)
}
SL.gbm.2 <- function(...) { 
  SL.gbm(..., interaction.depth = 2)
}
SL.gbm.3 <- function(...) { 
  SL.gbm(..., interaction.depth = 3)
}
SL.gbm.4 <- function(...) { 
  SL.gbm(..., interaction.depth = 4)
}
# RANDOM FOREST
tuneGrid <- expand.grid(mtry = c(500, 1000, 2200), nodesize = c(1, 5, 10))
for(mm in seq(nrow(tuneGrid))) { eval(parse(file = "", text =
                                              paste("SL.randomForest.", mm,
                                                    "<- function(..., mtry = ", tuneGrid[mm, 1], ", nodesize = ", tuneGrid[mm, 2], ") { SL.randomForest(..., mtry = mtry,
                                                    nodesize = nodesize) }", sep = "")))
  }


## LASSO
SL.glmnet.0.25 <- function(...) { 
  SL.glmnet(..., alpha = .25)
}

SL.glmnet.0.5 <- function(...) { 
  SL.glmnet(..., alpha = .5)
}

SL.glmnet.0.75 <- function(...) { 
  SL.glmnet(..., alpha = .75)
}

SL2.svm <- function(Y, X, newX, family, type.reg = "C-regression", type.class = "nu-classification",
        nu = 0.5,  n_controls_total = 1132116, ...) {
                   tsvm <- try(tune.svm (y=factor(Y), x = X, cost = c(1, 2), coef0 =c(-1, 0, 1), sampling = "cross", cross = 5))
                   if(class(tsvm) == "try-error"){
                   tsvm <- NULL
                   tsvm$best.parameters$cost <- 1
                   tsvm$best.parameters$coef0 <- 0
                   }
                   fit.svm <- svm(y=factor(Y), x = X, nu = mean(Y), type = type.reg, fitted = FALSE, probability = TRUE, kernel = "sigmoid", cost = tsvm$best.parameters$cost, coef0 = tsvm$best.parameters$coef0)
                   pred <- attr(predict(fit.svm, newdata = newX, probability = TRUE),
                   "prob")[, "1"]
                   fit <- list(object = fit.svm)
                   out <- list(pred = pred, fit = fit)
                   class(out$fit) <- c("SL.svm")
                   return(out)
                   }

# create.SL.glmnet <- function(alpha = c(0.25, 0.50, 0.75)) {
#   for(mm in seq(length(alpha))){
#     eval(parse(text = paste('SL.glmnet.', alpha[mm], '<- function(..., alpha = ', alpha[mm], ') SL.glmnet(..., alpha = alpha)', sep = '')), envir = .GlobalEnv)
#   }
#   invisible(TRUE)
# }

## LASSO
# SL.glmnet <- function(Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, nlambda = 100, useMin = TRUE, ...) {
#   .SL.require('glmnet')
#   # X must be a matrix, should we use model.matrix or as.matrix
#   if(!is.matrix(X)) {
#     X <- model.matrix(~ -1 + ., X)
#     newX <- model.matrix(~ -1 + ., newX)
#   }
#   # now use CV to find lambda
#   fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, lambda = NULL, type.measure = 'deviance', nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda)
#   # two options for lambda, fitCV$lambda.min and fitCV$lambda.1se
#   pred <- predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = 'response')
#   fit <- list(object = fitCV, useMin = useMin)
#   class(fit) <- 'SL.glmnet'
#   out <- list(pred = pred, fit = fit)
#   return(out)
# }
# 
# SL.glmnet.0.25 <- function(Y, X, newX, family, obsWeights, id, alpha = .25, nfolds = 10, nlambda = 100, useMin = TRUE, ...) {
#   .SL.require('glmnet')
#   # X must be a matrix, should we use model.matrix or as.matrix
#   if(!is.matrix(X)) {
#     X <- model.matrix(~ -1 + ., X)
#     newX <- model.matrix(~ -1 + ., newX)
#   }
#   # now use CV to find lambda
#   fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, lambda = NULL, type.measure = 'deviance', nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda)
#   # two options for lambda, fitCV$lambda.min and fitCV$lambda.1se
#   pred <- predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = 'response')
#   fit <- list(object = fitCV, useMin = useMin)
#   class(fit) <- 'SL.glmnet'
#   out <- list(pred = pred, fit = fit)
#   return(out)
# }
# 
# SL.glmnet.0.5 <- function(Y, X, newX, family, obsWeights, id, alpha = .5, nfolds = 10, nlambda = 100, useMin = TRUE, ...) {
#   .SL.require('glmnet')
#   # X must be a matrix, should we use model.matrix or as.matrix
#   if(!is.matrix(X)) {
#     X <- model.matrix(~ -1 + ., X)
#     newX <- model.matrix(~ -1 + ., newX)
#   }
#   # now use CV to find lambda
#   fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, lambda = NULL, type.measure = 'deviance', nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda)
#   # two options for lambda, fitCV$lambda.min and fitCV$lambda.1se
#   pred <- predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = 'response')
#   fit <- list(object = fitCV, useMin = useMin)
#   class(fit) <- 'SL.glmnet'
#   out <- list(pred = pred, fit = fit)
#   return(out)
# }
# 
# SL.glmnet.0.75 <- function(Y, X, newX, family, obsWeights, id, alpha = .75, nfolds = 10, nlambda = 100, useMin = TRUE, ...) {
#   .SL.require('glmnet')
#   # X must be a matrix, should we use model.matrix or as.matrix
#   if(!is.matrix(X)) {
#     X <- model.matrix(~ -1 + ., X)
#     newX <- model.matrix(~ -1 + ., newX)
#   }
#   # now use CV to find lambda
#   fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, lambda = NULL, type.measure = 'deviance', nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda)
#   # two options for lambda, fitCV$lambda.min and fitCV$lambda.1se
#   pred <- predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = 'response')
#   fit <- list(object = fitCV, useMin = useMin)
#   class(fit) <- 'SL.glmnet'
#   out <- list(pred = pred, fit = fit)
#   return(out)
# }