#' @title Automated Backward Stepwise GLM
#'
#' @description Takes in a dataframe and the dependent variable (in quotes) as arguments, splits the data into testing and training,
#' and uses automated backward stepwise selection to build a series of multiple regression models on the training data.
#' Each model is then evaluated on the test data and model evaluation metrics are computed for each model. These metrics
#' are provided as plots. Additionally, the model metrics are ranked and average rank is taken. The model with the best
#' average ranking among the metrics is displayed (along with its formula). By default, metrics are all given the same
#' relative importance (i.e., weights) when calculating average model metric rank, but if the user desires to give more
#' weight to one or more metrics than the others they can specify these weights as arguments (default for weights is 1).
#' As of v 0.2.0, only the family = gauissian(link = 'identity') argument is provided within the glm function.

#' @param data A dataframe with one column as the dependent variable and the others as independent variables
#' @param dv The column name of the (continuous) dependent variable (must be in quotes, i.e., 'Dependent_Variable')
#' @param aic_wt Weight given to the rank value of the AIC of the model fitted on the training data (used when calculating mean model performance, default = 1)
#' @param r_wt Weight given to the rank value of the Pearson Correlation between the predicted and actual values on the test data (used when calculating mean model performance, default = 1)
#' @param mae_wt Weight given to the rank value of Mean Absolute Error on the test data (used when calculating mean model performance, default = 1)
#' @param r_squ_wt Weight given to the rank value of R-Squared on the test data (used when calculating mean model performance, default = 1)
#' @param train_prop Proportion of the data used for the training data set
#' @param random_seed Random seed to use when splitting into training and testing data
#'
#' @return This function returns a plot for each metric by model and the best overall model with the formula used when fitting that model
#'
#' @examples
#' dt <- mtcars
#' stepwise_model <- backwd_stepwise_glm(data = dt,
#'                                       dv = 'mpg',
#'                                       aic_wt = 1,
#'                                       r_wt = 0.8,
#'                                       mae_wt = 1,
#'                                       r_squ_wt = 0.8,
#'                                       train_prop = 0.6,
#'                                       random_seed = 5)
#' stepwise_model
#' @export


# dv argument must be a string (i.e., it must have quotes)
backwd_stepwise_glm = function(data, dv, aic_wt=1, r_wt=1, mae_wt=1, r_squ_wt=1, train_prop = 0.7, random_seed = 7){

  # save data as dt
  dt = data

  # rename the given dv to 'DV' in dt
  names(dt)[names(dt) == dv] = 'DV'

  # Split into testing/training
  nrows = nrow(dt)
  train_size = floor(train_prop*nrows)
  set.seed(random_seed)
  train_idx = sample(nrows, train_size, replace = F)
  dt_train = dt[train_idx, ]
  dt_test = dt[-train_idx, ]

  # instantiate an empty vector to store model formulas
  model_formulas_vector = vector()
  # Instantiate an empty vector of AIC for which to append AIC values
  model_AIC_vector = vector()
  # instantiate empty vector to store predictions
  predictions_vector = vector()
  # instantiate an empty vector to store residuals
  residuals_vector = vector()
  # Instantiate an empty vector of pearson correlations between predicted and actual for which to append pearson r values
  model_r_vector = vector()
  # Instantiate an empty vector of R-squared for which to append R-Squared values
  model_R_Squared_vector = vector()
  # Instantiate an empty vector of MAE for which to append MAE values
  model_MAE_vector = vector()

  # store the first formula with every variable except the DV as a predictor
  formula = as.formula('DV ~ .')

  # instantiate a counter
  count = 0

  for (i in 1:(length(names(dt_test))-2)) { # -2 because we don't want a simple linear regression at the end
    # add 1 to the counter
    count = count + 1

    # fit model
    model = glm(formula, data = dt_train, family = gaussian(link = 'identity'))

    # save the model formula to a vector
    model_formula = model$formula
    # add model_formula to model_formulas_vector
    model_formulas_vector = c(model_formulas_vector, model_formula)

    # get aic
    model_aic = model$aic
    # append this to the vector model_first_aic
    model_AIC_vector = c(model_AIC_vector, model_aic)

    # get predictions
    predictions = predict(model, dt_test)
    # append to predictions_vector
    predictions_vector = c(predictions_vector, predictions)

    # get the pearson correlation between predictions and actual values
    correlation = cor(predictions, dt_test$DV)
    # append correlation to model_r_vector
    model_r_vector = c(model_r_vector, correlation)

    # get the r-squared value
    r_squared = 1 - (sum((dt_test$DV-predictions)^2)/sum((dt_test$DV-mean(dt_test$DV))^2))
    # append to model_R_Squared_vector
    model_R_Squared_vector = c(model_R_Squared_vector, r_squared)

    # get residuals
    resids = dt_test$DV - predictions
    # append to residuals_vector
    residuals_vector = c(residuals_vector, resids)

    # get the absolute values of residuals
    abs_residuals = abs(resids)
    # get mean of the absolute values of residuals (MAE)
    MAE = mean(abs_residuals)
    # Append to the model_MAE_vector
    model_MAE_vector = c(model_MAE_vector, MAE)

    # Get variable importance
    #library(caret)
    var_imp = as.data.frame(varImp(model))
    # Turn into a data frame
    imp = data.frame(overall = var_imp$Overall, names = rownames(var_imp))
    # Order it by importance
    imp_ordered = imp[order(imp$overall, decreasing = TRUE),]
    # Get only names from this vector with the last name removed
    names_vector = as.vector(imp_ordered[1:nrow(imp_ordered)-1,2])
    # remove quotes and add plus sig for formula
    names_no_quotes = noquote(paste(names_vector, collapse='+'))

    # save the new formula
    formula = as.formula(paste(paste('DV ~'), names_no_quotes))
  }

  # Print a message if the loop was completed successfuly
  if (count == length(dt_test) - 2) {
    print(noquote(paste("Success: There were", count, 'models build and there were', length(dt_test) - 1,
                        'independent variables in the origninal model')))
  } else
    print(noquote('Failure: There was an error, review code'))

  # rank the values in each of the vectors (1=worst, largest number = best)
  # AIC: Lower is better
  model_AIC_vector_rank = rank(model_AIC_vector) * aic_wt
  # Pearson r: Higher is better
  model_r_vector_rank = rank(-model_r_vector) * r_wt
  # R-Squared: Higher is better
  model_R_Squared_vector_rank = rank(-model_R_Squared_vector) * r_squ_wt
  # MAE: Lower is better
  model_MAE_vector_rank = rank(model_MAE_vector) * mae_wt

  # create a data frame with these ranked vectors
  dt_rank = data.frame(model_AIC_vector_rank, model_r_vector_rank, model_R_Squared_vector_rank, model_MAE_vector_rank)

  # get an average ranking for each model
  dt_rank$Avg_Rank = rowMeans(dt_rank, na.rm = FALSE, dims = 1)

  # which row (i.e., model has lowest ranking)
  best_model = which.min(dt_rank$Avg_Rank)

  # insert plots
  # AIC plot
  # get the index of the minimum AIC value to programmatically plot it
  best_model_AIC = which.min(model_AIC_vector)
  # find the value of the lowest AIC
  lowest_aic = min(model_AIC_vector)
  # round this value
  rounded_lowest_aic = round(lowest_aic, 3)
  # Plot the AIC (y) by model number (x)
  ymin = min(model_AIC_vector, na.rm = TRUE)
  ymax = max(model_AIC_vector, na.rm = TRUE)
  plot(model_AIC_vector, type="o", col="blue", ylim=c(ymin, ymax), ylab='Model AIC', xlab='Model Number',
       main= print(paste('Model', best_model_AIC, 'has the lowest AIC:', rounded_lowest_aic)))

  # Pearson correlation
  # get the index of the minimum AIC value to programmatically plot it
  best_model_r = which.max(model_r_vector)
  # find the value of the highest r
  highest_r = max(model_r_vector)
  # round this value
  rounded_highest_r = round(highest_r, 4)
  # Plot the r (y) by model number (x)
  ymin = min(model_r_vector, na.rm = TRUE)
  ymax = max(model_r_vector, na.rm = TRUE)
  plot(model_r_vector, type="o", col="blue", ylim=c(ymin, ymax), ylab='Model r', xlab='Model Number',
       main= print(paste('Model', best_model_r, 'has the highest r:', rounded_highest_r)))

  # MAE
  # get the index of the minimum MAE value to programmatically plot it
  best_model_MAE = which.min(model_MAE_vector)
  # find the value of the minimum MAE
  lowest_MAE = min(model_MAE_vector)
  # round this value
  rounded_lowest_MAE = round(lowest_MAE, 4)
  # Plot the MAE (y) by model number (x)
  ymin = min(model_MAE_vector, na.rm = TRUE)
  ymax = max(model_MAE_vector, na.rm = TRUE)
  plot(model_MAE_vector, type="o", col="blue", ylim=c(ymin, ymax), ylab='Model MAE', xlab='Model Number',
       main= print(paste('Model', best_model_MAE, 'has the lowest MAE:', rounded_lowest_MAE)))

  # R squared
  # get the index of the maximum r squared value to programmatically plot it
  best_model_R_squared = which.max(model_R_Squared_vector)
  # find the value of the maximum R-squared
  highest_r_squared = max(model_R_Squared_vector)
  # round this value
  rounded_highest_r_squared = round(highest_r_squared, 4)
  # Plot the R-Squared (y) by model number (x)
  ymin = min(model_R_Squared_vector, na.rm = TRUE)
  ymax = max(model_R_Squared_vector, na.rm = TRUE)
  plot(model_R_Squared_vector, type="o", col="blue", ylim=c(ymin, ymax), ylab='Model R-Squared', xlab='Model Number',
       main= print(paste('Model', best_model_R_squared, 'has the highest R-Squared:', rounded_highest_r_squared)))

  # create vector from Avg_Rank
  Avg_Rank_vector = as.vector(dt_rank$Avg_Rank)
  # Plot mean rank of each model
  ymin = min(Avg_Rank_vector, na.rm = TRUE)
  ymax = max(Avg_Rank_vector, na.rm = TRUE)
  plot(Avg_Rank_vector, type="o", col="blue", ylim=c(ymin, ymax), ylab='Avg. Weighted Rank', xlab='Model Number',
       main = print(paste('Model', best_model, 'has the lowest weighted average rank (lower is better)\n amongst model evaluation metrics')))

  # Build the model with the best score and return it, so the user can call summary
  model = glm(paste(model_formulas_vector[best_model]), data = dt_train, family = gaussian(link = 'identity'))
  return(model)
}






