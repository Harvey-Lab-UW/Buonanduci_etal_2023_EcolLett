# Requires modified prediction function
source("Scripts/predict.gcrq.inputform.R")

# Quantile loss function for evaluating model fit
loss_fun <- function(tau, obs, pred){
  mean( pmax(tau * (obs - pred), (tau - 1) * (obs - pred)) )
}


# form = model formula
# taus = vector of taus to evaluate
# data = dataframe, which must include 'subsample' variable
# n_s = number of subsamples, or folds

kfold_fun <- function(form, taus, data, n_s = 10){

  results <- list()
  
  for (k in 1:n_s){
    
    for(i in 1:length(taus)){
      # Fit model
      mod <- gcrq(form, tau = taus[i], data = filter(data, subsample != k))
      # Predict from model
      pred <- predict.gcrq.inputform(mod, form = form,
                                     newdata = as.data.frame(filter(data, subsample == k)))
      # Calculate loss
      loss <- loss_fun(tau = taus[i], 
                       obs = filter(data, subsample == k)[[ all.vars(form)[1] ]], 
                       pred = pred)
      
      
      # Store tau, FRG, model, loss
      results[[(k-1)*length(taus) + i]] <- tibble(resp = all.vars(form)[1],
                                                  tau = taus[i],
                                                  loss = loss)
    }
  }
  return( bind_rows(results) )
}


# Alternate function for fitting and predicting from model
# fit to two separate datasets
# form = model formula
# taus = vector of taus to evaluate
# data = dataframe, which must include 'subsample' variable
# n_s = number of subsamples, or folds

kfold_2d_fun <- function(form, taus, data1, data2, n_s = 10){
  
  results <- list()
  
  for (k in 1:n_s){

    for(i in 1:length(taus)){
      # Fit models
      mod1 <- gcrq(form, tau = taus[i], data = filter(data1, subsample != k))
      mod2 <- gcrq(form, tau = taus[i], data = filter(data2, subsample != k))
      # Predict from models
      pred1 <- predict.gcrq.inputform(mod1, form = form,
                                      newdata = as.data.frame(filter(data1, subsample == k)))
      pred2 <- predict.gcrq.inputform(mod2, form = form,
                                      newdata = as.data.frame(filter(data2, subsample == k)))
      # Calculate loss
      loss <- loss_fun(tau = taus[i], 
                       obs = c(filter(data1, subsample == k)[[ all.vars(form)[1] ]],
                               filter(data2, subsample == k)[[ all.vars(form)[1] ]]), 
                       pred = c(pred1, pred2))
      
      # Store tau, FRG, model, loss
      results[[(k-1)*length(taus) + i]] <- tibble(resp = all.vars(form)[1],
                                                  tau = taus[i],
                                                  loss = loss)
    }
  }
  return( bind_rows(results) )
}
