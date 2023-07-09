#' Extract One-Hot Encoded (OHE) Features from Peptide Sequences
#'
#' This function takes a data frame or a vector of peptide sequences and generates
#' a one-hot encoded data frame representing each amino acid in the sequences.
#' It can also include additional data (such as docking information), if provided.
#'
#' @param df A data frame or a vector of peptide sequences.
#' @param sequence_col A string representing the name of the column containing the peptide sequences.
#' @param docking_col A string representing the name of the column containing the docking information.
#' @return A data frame containing one-hot encoded peptide sequences and, if provided, docking information.
#' @examples
#' # Load required library caret
#' library(caret)
#' # Generate a mock data frame of peptide sequences
#' df <- data.frame(Sequence = c("AVILG", "VILGA", "ILGAV", "LGAVI"), X = c(1,1,2,3))
#' # Apply the function to the mock data
#' extract_features_OHE(df)
#' @export
#'
extract_features_OHE <- function(df, sequence_col = "Sequence", docking_col = NULL) {

    # Load the required library
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is needed for this function to work. Please install it.")
  }

  # Check if data is a vector. If it is, convert it to a data frame.
  if(is.vector(df)) {
    data <- data.frame(Sequence = df)
    sequence_col <- "Sequence"
  } else if (!is.data.frame(df)) {
    stop("'data' must be a data frame or a vector.")
  }

  # Check if the sequence column exists
  if (!(sequence_col %in% names(df))) {
    stop(paste("Column", sequence_col, "not found in the data frame."))
  }

  # Check if the docking column exists if provided
  if (!is.null(docking_col) && !(docking_col %in% names(df))) {
    stop(paste("Column", docking_col, "not found in the data frame."))
  }

  # Split sequences into individual amino acids and convert list of sequences to a data frame
  sequences <- df[[sequence_col]]
  split_sequences <- strsplit(sequences, "")
  sequence_df <- do.call("rbind", lapply(split_sequences, function(sequence) {
    data.frame(t(sequence))
  }))

  # One-hot encode the sequences
  dummies <- caret::dummyVars(~ ., data = sequence_df)
  one_hot_encoded_sequences <- data.frame(stats::predict(dummies, newdata = sequence_df))

  # Move the specified sequence column to the start of the data frame
  one_hot_encoded_sequences[[sequence_col]] <- df[[sequence_col]]

  if(!is.null(docking_col)) {
    # Include docking data if it is provided
    one_hot_encoded_sequences[[docking_col]] <- df[[docking_col]]
    new_col_order <- c(sequence_col, docking_col, setdiff(names(one_hot_encoded_sequences), c(sequence_col, docking_col)))
  } else {
    # If no docking data, only include the sequence data
    new_col_order <- c(sequence_col, setdiff(names(one_hot_encoded_sequences), sequence_col))
  }

  one_hot_encoded_sequences <- one_hot_encoded_sequences[,new_col_order]

  return(one_hot_encoded_sequences)
}
