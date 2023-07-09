#' Reduce Peptide Sequences by One Residue
#'
#' This function takes a vector of peptide sequences and generates all possible
#' sequences by removing one amino acid residue at a time. It can also associate each
#' sequence with an ID, if provided.
#'
#' @param peptides A character vector of peptide sequences.
#' @param id A character vector of IDs that correspond to the peptides.
#' @return A list of data frames, each containing all possible sequences resulting
#'         from removing one amino acid from the original sequence.
#' @examples
#' # Generate a mock vector of peptide sequences
#' peptides <- c("AVILG", "VILGA", "ILGAV", "LGAVI")
#' # Apply the function to the mock data
#' reduce_sequences(peptides)
#' @export
#'
reduce_sequences <- function(peptides, id = NULL) {
  # Check if peptides is a vector
  if(!is.vector(peptides)) {
    stop("Input 'peptides' should be a vector of peptide sequences.")
  }

  # Check if id is either NULL or a vector of the same length as peptides
  if(!is.null(id) && (length(id) != length(peptides))) {
    stop("If provided, 'id' should be a vector of the same length as 'peptides'.")
  }

  # Initialize list to store results
  reduced_peptides <- list()

  # Iterate over each peptide
  for(i in 1:length(peptides)) {
    # Check if the sequence is long enough to be reduced
    if(nchar(peptides[i]) < 2) {
      stop("All sequences should be at least two residues long.")
    }

    # Generate all possible reduced peptides by removing each amino acid sequentially
    # If an id vector is provided, include it in the data frame
    if(is.null(id)) {
      reduced_peptides[[i]] <- data.frame(Sequence = sapply(1:nchar(peptides[i]), function(j) paste0(substr(peptides[i], 1, j-1), substr(peptides[i], j+1, nchar(peptides[i])))),
                                          stringsAsFactors = FALSE)
    } else {
      reduced_peptides[[i]] <- data.frame(ID = id[i],
                                          Sequence = sapply(1:nchar(peptides[i]), function(j) paste0(substr(peptides[i], 1, j-1), substr(peptides[i], j+1, nchar(peptides[i])))),
                                          stringsAsFactors = FALSE)
    }
  }

  return(reduced_peptides)
}
