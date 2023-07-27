#' Peptide Docking Prediction Function
#'
#' This function performs model training and prediction on a given set of peptide sequences.
#' It can deal with both binary classification and regression tasks, with support for both QSAR and OHE property types.
#' The function makes use of LightGBM for model training and prediction.
#'
#' @param x A dataframe containing the sequences and docking information
#' @param docking_col A string or integer indicating the column name or number in data frame 'x' for docking scores.
#' @param seq_col A string or integer indicating the column name or number in data frame 'x' for sequences.
#' @param train_pc A numeric value indicating the percentage of the data to be used for training in the LightGBM algorithm.
#' @param properties A string indicating the type of properties to be used (either 'QSAR' or 'OHE')
#' @param method A string indicating the prediction method (either 'All' for predicting a full library - all tetramers, all pentamers ... determined by n - or 'Custom' if 'to_pred' is given)
#' @param n An integer indicating the length of the sequences (tetramers = 4, pentamers = 5, and so on). All sequences must be the same length.
#' @param to_pred A dataframe or vector of sequences to be predicted. If a vector is provided the features will be extracted. If a data frame es provided it must contain the same features as the training data.
#' @param objective A string indicating the objective for LightGBM model ('binary', 'regression_l1', 'regression_l2' are accepted for now. Other ones are available if the raw lightgbm package is used)
#' @param metric A string indicating the evaluation metric for the LightGBM model (eg. 'auc', 'l1', 'l2')
#' @param best_pc A numeric value indicating the percentage of 'best-performers' groups for binary objective.
#' @param nround An integer indicating the maximum number of iterations for the LightGBM model
#' @param early_stop An integer indicating the number of rounds without improvement after which training will stop
#' @param seed An integer to set the seed for random processes for reproducibility
#' @param monte_carlo An integer indicating the number of times the prediction process is run.
#' @return A list containing average prediction, all predictions, all models and predicted sequences
#' @examples
#' # To use the function
#' # Generate random sequences of 4 letters
#' sequences <- unique(replicate(200, paste(sample(c("A", "C", "T", "G"), 4,
#' replace = TRUE), collapse = "")))
#' # Generate random docking scores
#' dockings <- rnorm(length(sequences), 0, 1)
#' df <- data.frame(
#' Sequence = sequences,
#' Docking = dockings)
#' # Predict a vector without features extracted
#' to_pred <- c("GTAC", "ACTG", "HFGL", "REST", "DYLN")
#' results <- run_prediction(x = df, docking_col = "Docking",
#' seq_col = "Sequence", train_pc = 0.75, properties = "QSAR", method = "Custom",
#' n = 4, to_pred = to_pred, objective = "regression_l1", metric = "l1",
#' best_pc = 0.2, nround = 100, early_stop = 5, seed = 1234, monte_carlo = 2)
#' @export
#'

run_prediction <- function(x, docking_col = "X", seq_col = "Sequence", train_pc = .75, properties = "OHE",
                                method = "Custom", n = NULL, to_pred = NULL, objective = "binary", metric = "auc", best_pc = .2, nround = 1000, early_stop = 50, seed = 1234, monte_carlo = 10) {

  # Load the required library
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is needed for this function to work. Please install it.")
  }

  # Validating method and properties arguments
  if(!method %in% c("All", "Custom")){
    stop("Method should be 'All' to predict the full n-library; or 'Custom' to predict a custom vector of peptides supplied in 'to_pred' argument." )
  }
  if(!properties %in% c("QSAR", "OHE")){
    stop("The properties values supported are 'QSAR' or 'OHE'.")
  }

  # Setting the seed for reproducibility
  set.seed(seed)

  # Checking for NA or empty values
  if (any(is.na(x[,seq_col]) | x[,seq_col] == "")) {
    stop("Sequence column contains NA or empty values.")
  }

  # Select only the columns we need
  x <- x[,c(seq_col, docking_col)]

  # Sort by docking_col
  x <- x[order(x[,docking_col]),]
  rownames(x) <- x[,seq_col]

  # If objective is binary, create binary groups
  if(objective == "binary"){
    # Validate best_pc argument. Should be between 0 and 1.
    if(best_pc <= 0 | best_pc >= 1){
      stop("Percentage for 'best-performing' group must be a valid integer between 0 and 1 for binary objective.")
    }
    # Creating the binary groups
    x1 <- rep(nrow(x)*best_pc,1)
    x2 <- rep(nrow(x)*(1-best_pc),0)
    y <- data.frame(y = c(x1, x2))
    y <- rbind(y, data.frame(y = rep(0, nrow(x) - nrow(y))))

    # Updating the data frame x
    x <- data.frame(TARGET = y$y, Sequence = x[,seq_col])

    # Printing the group sizes
    cat("Worse-performers group '1', size = ", summary(as.factor(x$TARGET))[1],
        "\nBest-performers group '0', size = ", summary(as.factor(x$TARGET))[2], "\n")
  }

  # If the objective is regression, we prepare the data accordingly
  if(objective %in% c("regression_l1", "regression_l2")){
    x <- data.frame(TARGET = x[,docking_col], Sequence = x[,seq_col])
  }

  cat("Extracting '", properties, "' properties. Refer to the documentation to see other options.\n")

  # Extracting features based on property type
  if(properties == "QSAR") {
    extract_prop <- extract_features_QSAR(n = n, custom.list = T, PeList = x$Sequence, norm = T)
  } else if(properties == "OHE") {
    extract_prop <- extract_features_OHE(df = x)
  }else {
    stop("Properties value is invalid. It should be either 'QSAR' or 'OHE'.")
  }

  # Add the target column and rearrange the columns
  extract_prop <- cbind(extract_prop, TARGET = x$TARGET)
  extract_prop <- extract_prop[,c(1, ncol(extract_prop), 2:(ncol(extract_prop)-1))]

  # Check 'n' validity
  if(!is.null(n)){
    if(n < 2 | n > 7) {
      stop("This function only support n between 2 and 7.")
    } else{
      if (n > 4 & method == "All") {
        cat("Peptide library to predict = ", 20^n , ". This might take a lot of time.")
      }
    }
  }

  # Extract the features for the peptides to be predicted
  if (method == "All") {
    if (properties == "QSAR") {
      to_pred <- extract_features_QSAR(n=n, norm = T)
    } else if (properties == "OHE") {
      to_pred <- extract_features_OHE(n = n) # Pass the computed peptide combinations to OHE function
    }
  }

  if (method == "Custom") {
    if(is.null(to_pred)){
      if (properties == "QSAR") {
        to_pred <- extract_features_QSAR(n=n, custom.list = T, PeList = to_pred, norm = T)
      } else if (properties == "OHE") {
        to_pred <- extract_features_OHE(df = to_pred) # Pass the computed peptide combinations to OHE function
      }else {
        stop("Method should be 'All' or 'Custom'. If 'Custom' is set, you need to provide a character vector in 'to_pred'.")
      }
    }else {
      if (is.vector(to_pred)) {
        if (properties == "QSAR") {
          to_pred <- extract_features_QSAR(n=n, custom.list = T, PeList = to_pred, norm = T)
        } else if (properties == "OHE") {
          to_pred <- extract_features_OHE(df = to_pred) # Pass the computed peptide combinations to OHE function
        }
      }
    }
  }


  # Preparing the prediction matrix and sequence list
  X_topred <- to_pred$Sequence
  Y_topred <- as.matrix(to_pred[,-1])

  # Initialize empty list to store predictions
  pred_list <- list()
  models_list <- list()

  # Run the prediction process multiple times according to monte_carlo
  for(mc in 1:monte_carlo) {
    tries <- 0
    while(tries < 200) {
      # Sample training indices and split the data
      train_idx <- sample(seq_along(extract_prop$TARGET), train_pc * nrow(extract_prop))

      # Splitting the data into training and testing set
      train_data <- extract_prop[train_idx, ]
      test_data <- extract_prop[-train_idx, ]

      # Preparing the training and testing matrices
      X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Sequence"))])
      y_train <- train_data$TARGET

      X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Sequence"))]
      y_test <- test_data$TARGET

      # Checks for y_train and y_test to ensure they both have more than one unique value
      if (length(unique(y_train)) > 1 & length(unique(y_test)) > 1) {
        break
      }

      tries <- tries + 1
    }

    if(tries == 20) {
      stop("Failed to generate distinct y_train and y_test after 200 attempts.")
    }

    # Define the parameters for the lightgbm model
    params <- list(
      objective = objective,
      metric = metric,
      boosting = "gbdt")

    # Prepare the training and validation datasets
    dtrain = lightgbm::lgb.Dataset(data = as.matrix(X_train), label = y_train)
    valids = list(validation = lightgbm::lgb.Dataset(data = as.matrix(X_test), label = y_test))

    # Training the lightgbm model
    model <- lightgbm::lgb.train(
      data = dtrain,
      valids = valids,
      params = params,
      nrounds = nround,
      early_stopping_rounds = early_stop)

    pred <- stats::predict(model, as.matrix(Y_topred), reshape = T)

    prediction <- data.frame(Sequence = X_topred, pred)

    # Store the prediction for this iteration
    pred_list[[mc]] <- prediction
    models_list[[mc]] <- model
  }

  # Check if all predictions have the same sequences in the same order
  if(!all(sapply(pred_list, function(df) all(df$Sequence == pred_list[[1]]$Sequence)))) {
    stop("Sequences in predictions do not match.")
  }

  # Calculate the average prediction
  avg_pred <- rowMeans(do.call(cbind, lapply(pred_list, `[[`, "pred")))

  # Combine the average prediction with the sequences
  avg_prediction <- data.frame(Sequence = pred_list[[1]]$Sequence, Pred = avg_pred)
  avg_prediction <- avg_prediction[order(avg_prediction[,2], decreasing = T),]

  # Create an empty data frame to store the scores
  models_scores <- data.frame(
    model = character(),
    score = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop over the models and add their scores to the data frame
  for (i in 1:mc) {
    models_scores <- rbind(models_scores,
                           data.frame(model = i,
                                      score = models_list[[i]][["best_score"]],
                                      stringsAsFactors = FALSE)
    )
  }

  # Find the model with the highest score
  highest_score_index <- which.max(models_scores$score)
  highest_score_model <- models_scores$model[highest_score_index]
  highest_score <- models_scores$score[highest_score_index]

  # Find the model with the lowest score
  lowest_score_index <- which.min(models_scores$score)
  lowest_score_model <- models_scores$model[lowest_score_index]
  lowest_score <- models_scores$score[lowest_score_index]

  # Print the result
  cat("Model", highest_score_model, "has the highest", metric, "=", highest_score,
      ", model", lowest_score_model, "has the lowest", metric, "=", lowest_score, ".")

  rm(models_scores)

  if(method == "Custom"){
    return(list(avg_prediction = avg_prediction, mc_predictions = pred_list, models = models_list, to_pred = "Custom"))

  }else{
    return(list(avg_prediction = avg_prediction, mc_predictions = pred_list, models = models_list, predicted_matrix = to_pred))
  }

}

