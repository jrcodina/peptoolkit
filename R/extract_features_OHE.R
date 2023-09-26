#' Extract One-Hot Encoded (OHE) Features from Peptide Sequences
#'
#' This function takes a data frame or a vector of peptide sequences and generates
#' a one-hot encoded data frame representing each amino acid in the sequences.
#' It can also include additional data (such as docking information), if provided.
#' Furthermore, it can generate a peptide library of specified length n.
#'
#' @param df A data frame or a vector of peptide sequences. If 'df' is provided, 'n' will be ignored.
#' @param sequence_col A string representing the name of the column containing the peptide sequences.
#' @param docking_col A string representing the name of the column containing the docking information.
#' @param n An integer representing the length of the peptide library to be generated. If 'df' is provided, 'n' will be ignored.
#' @return A data frame containing one-hot encoded peptide sequences and, if provided, docking information.
#' @export
#' @import caret
#'
#' @examples
#' # Load required library caret
#' library(caret)
#' extract_features_OHE(df = c('ACA', 'EDE'))
extract_features_OHE <- function(df = NULL, sequence_col = "Sequence", docking_col = NULL, n = NULL) {

  # Load the required library
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is needed for this function to work. Please install it.")
  }
  if (!requireNamespace("Peptides", quietly = TRUE)) {
    stop("Package 'Peptides' is needed for this function to work. Please install it.")
  }

  # Check if df is provided and is valid
  if (!is.null(df)) {
    if (is.vector(df)) {
      # Convert vector to data frame.
      df <- data.frame(Sequence = df)
    } else if (!is.data.frame(df)) {
      stop("'df' must be a data frame or a vector.")
    }
  } else if (!is.null(n)) {
    # If df is not provided, check for n

    # Check if n is a positive integer
    if(!is.numeric(n) || n <= 2 || (n %% 1 != 0)){
      stop("n must be a positive integer greater than 2.")
    }

    # Generate the peptide library
    PeList <- expand.grid(rep(list(Peptides::aaList()), n))
    sequences <- do.call(paste0, PeList)
    df <- data.frame(Sequence = sequences)
  } else {
    stop("Either 'df' or 'n' must be provided.")
  }

  # Store column names
  col_names <- names(df)

  # Check if the sequence column exists
  if (!(sequence_col %in% col_names)) {
    stop(paste("Column", sequence_col, "not found in the data frame."))
  }

  # Check if the docking column exists if provided
  if (!is.null(docking_col) && !(docking_col %in% col_names)) {
    stop(paste("Column", docking_col, "not found in the data frame."))
  }

  # Determine the length of the longest sequence
  sequences <- toupper(df[[sequence_col]])
  max_length <- max(nchar(sequences))

  # Warn once if sequences need padding
  if(any(nchar(sequences) < max_length)) {
    warning("Not all sequences are of the same length. Padding shorter sequences with 'Z'.")
  }

  # Pad shorter sequences with "Z"
  padded_sequences <- sapply(sequences, function(sequence) {
    if (nchar(sequence) < max_length) {
      return(paste(c(strsplit(sequence, "")[[1]], rep("Z", max_length - nchar(sequence))), collapse = ""))
    } else {
      return(sequence)
    }
  })

  # Split sequences into individual amino acids and convert list of sequences to a data frame
  split_sequences <- strsplit(padded_sequences, "")
  sequence_df <- do.call("rbind", lapply(split_sequences, function(sequence) {
    data.frame(t(sequence))
  }))

  # Check for columns with only one level
  check_single_levels <- function(df){
    single_levels <- sapply(df, function(x) length(unique(x)) < 2)
    return(names(df)[single_levels])
  }

  single_level_cols <- check_single_levels(sequence_df)
  if(length(single_level_cols) > 0) {
    warning(paste("Removing columns with only one level:", paste(single_level_cols, collapse = ", ")))
    sequence_df <- sequence_df[, !(names(sequence_df) %in% single_level_cols)]
  }

  # Create dummy variable for each position and amino acid (including 'Z')
  dummies <- caret::dummyVars(~ ., data = sequence_df)
  one_hot_encoded_sequences <- data.frame(stats::predict(dummies, newdata = sequence_df))

  # Create a "template" matrix with all possible positions and residues
  all_residues <- c(Peptides::aaList(), if(any(nchar(sequences) < max_length)) "Z")
  template_names <- unlist(sapply(1:7, function(x) paste0("X", x, all_residues)))

  # Ensure all template columns exist in the output, filled with 0 if not present in the data
  missing_cols <- setdiff(template_names, names(one_hot_encoded_sequences))
  if (length(missing_cols) > 0) {
    one_hot_encoded_sequences[, missing_cols] <- 0
  }

  # Move the specified sequence column to the start of the data frame
  one_hot_encoded_sequences[[sequence_col]] <- df[[sequence_col]]

  # Define new column order
  new_col_order <- c(sequence_col, setdiff(names(one_hot_encoded_sequences), sequence_col))

  if(!is.null(docking_col)) {
    # Include docking data if it is provided
    one_hot_encoded_sequences[[docking_col]] <- df[[docking_col]]
    new_col_order <- c(sequence_col, docking_col, setdiff(names(one_hot_encoded_sequences), c(sequence_col, docking_col)))
  }

  one_hot_encoded_sequences <- one_hot_encoded_sequences[,new_col_order]

  return(one_hot_encoded_sequences)
}
