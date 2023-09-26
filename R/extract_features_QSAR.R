#' Extract QSAR Features from Peptide Sequences
#'
#' This function extracts various Quantitative Structure-Activity Relationship (QSAR) features
#' from peptide sequences. The extraction is based on a variety of amino acid properties
#' and functions from the "Peptides" package (https://github.com/dosorio/Peptides/).
#'
#' @param df A data frame or a vector of peptide sequences. If 'df' is provided, 'n' will be ignored.
#' @param sequence_col A string representing the name of the column containing the peptide sequences.
#' @param docking_col A string representing the name of the column containing the docking information.
#' @param n An integer representing the length of the peptide library to be generated. If 'df' is provided, 'n' will be ignored.
#' @param pH The pH used for calculating charge (default is 7.4).
#' @param normalize A boolean indicating if the data should be normalized (default is FALSE).
#' @return A dataframe with the calculated peptide properties.
#' @export
#' @import Peptides
#'
#' @examples
#' extract_features_QSAR(df = c('ACA', 'EDE'))
extract_features_QSAR <- function(df = NULL, n = NULL, sequence_col = "Sequence", docking_col = NULL, pH = 7.4, normalize = FALSE){

  if (!requireNamespace("Peptides", quietly = TRUE)) {
    stop("Package 'Peptides' is needed for this function to work. Please install it.")
  }

  # Check if df is provided and is valid
  if (!is.null(df)) {
    if (is.data.frame(df)){
      if (length(df) >=2 & is.null(docking_col)) {
        warning(paste0("No 'docking_col' argument provided, but 'df' contains more than 1 column. Assuming '", names(df)[2], "' as 'docking_col'."))
        docking_col <- names(df)[2]
      }
    }

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

  sequences <- toupper(df[[sequence_col]])

  peptides <- cbind(df,
                    as.data.frame(do.call(rbind, Peptides::crucianiProperties(sequences))),
                    as.data.frame(do.call(rbind, Peptides::fasgaiVectors(sequences))),
                    as.data.frame(do.call(rbind, Peptides::kideraFactors(sequences))),
                    as.data.frame(do.call(rbind, Peptides::protFP(sequences))),
                    as.data.frame(do.call(rbind, Peptides::tScales(sequences))),
                    as.data.frame(do.call(rbind, Peptides::vhseScales(sequences))),
                    as.data.frame(do.call(rbind, Peptides::zScales(sequences))))

  pI <- cbind(pI_EMBOS = Peptides::pI(sequences, pKscale = "EMBOSS"),
              pI_Bjellqvist = Peptides::pI(sequences, pKscale = "Bjellqvist"),
              pI_Lehninger = Peptides::pI(sequences, pKscale = "Lehninger"),
              pI_Murray = Peptides::pI(sequences, pKscale = "Murray"),
              pI_Rodwell = Peptides::pI(sequences, pKscale = "Rodwell"),
              pI_Sillero = Peptides::pI(sequences, pKscale = "Sillero"))

  mShft_15n <- Peptides::massShift(sequences, label = "15n", aaShift = NULL, monoisotopic = TRUE)

  charges <- cbind(Peptides::charge(sequences, pH = pH, pKscale = "EMBOSS"))

  hydrophobicity_indexes <- c("Aboderin", "AbrahamLeo", "BlackMould",
                              "BullBreese", "Casari", "Chothia", "Cid",
                              "Cowan3.4", "Cowan7.5", "Eisenberg", "Fasman",
                              "Fauchere", "Goldsack", "Guy", "HoppWoods",
                              "interfaceScale_pH2", "interfaceScale_pH8", "Janin",
                              "Jones", "Juretic", "Kuhn", "KyteDoolittle",
                              "Levitt", "Manavalan", "Miyazawa", "octanolScale_pH2",
                              "oiScale_pH2", "oiScale_pH8", "Parker", "Ponnuswamy",
                              "Prabhakaran", "Rao", "Rose", "Roseman", "Sweet",
                              "Tanford", "Welling", "Wilson", "Wolfenden",
                              "Zimmerman")

  hydrophobicity <- list()

  for (i in 1:length(hydrophobicity_indexes)) {

    hydrophobicity[[i]] <- Peptides::hydrophobicity (sequences, scale = hydrophobicity_indexes[i])

    names(hydrophobicity)[[i]] <- paste0("hb_",hydrophobicity_indexes[i])

  }


  peptides <- cbind(peptides, data.frame(aIndex = Peptides::aIndex(sequences),
                                         Boman = Peptides::boman(sequences),
                                         chrg_EMBOSS = charges,
                                         hmoment1 = Peptides::hmoment(sequences, angle = 100, window = 11),
                                         hmoment2 = Peptides::hmoment(sequences, angle = 160, window = 11),
                                         hydrophobicity,
                                         instaIndex = Peptides::instaIndex(sequences),
                                         mShft_15n,
                                         mw1 = Peptides::mw(sequences, monoisotopic = TRUE),
                                         pI))
  if (normalize) {
    X <- peptides[,!(names(peptides) %in% c(sequence_col, docking_col))]

    max = apply(X , 2 , max)
    min = apply(X, 2 , min)
    Xn = as.data.frame(scale(X, center = min, scale = max - min))

    peptides <- Xn
  }

  # Move the specified sequence column to the start of the data frame
  QSAR_sequences <- data.frame(Sequence = sequences)

  if(!is.null(docking_col)) {
    # Include docking data if it is provided
    QSAR_sequences[[docking_col]] <- df[[docking_col]]
  }

  QSAR_sequences <- cbind(QSAR_sequences, peptides[-1])

  return(QSAR_sequences)
}
