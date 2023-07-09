#' Extract QSAR Features from Peptide Sequences
#'
#' This function extracts various Quantitative Structure-Activity Relationship (QSAR) features
#' from peptide sequences. The extraction is based on a variety of amino acid properties
#' and functions from the "Peptides" package (https://github.com/dosorio/Peptides/).
#'
#' @param n The length of the peptide sequences.Must be more than 2.
#' @param pH The pH used for calculating charge (default is 7.4).
#' @param custom.list A boolean indicating if a custom peptide list is provided (default is FALSE).
#' @param PeList The custom list of peptides (required if custom.list is TRUE).
#' @param rem.cys A boolean indicating if sequences with Cys should be removed (default is FALSE).
#' @param rem.met A boolean indicating if sequences with Met should be removed (default is FALSE).
#' @param rem.sali A boolean indicating if sequences with 2 or more small aliphatic amino acids should be removed (default is FALSE).
#' @param norm A boolean indicating if the data should be normalized (default is FALSE).
#'
#' @return A dataframe with the calculated peptide properties.
#' @export
#'
#' @examples
#' extract_features_QSAR(n = 3)
#' extract_features_QSAR(n = 3, custom.list = TRUE, PeList = c('ACA', 'ADE'))
extract_features_QSAR <- function(n, pH = 7.4, custom.list = FALSE, PeList = NULL, rem.cys = FALSE, rem.met = FALSE, rem.sali = FALSE, norm = FALSE){

  # Check if custom.list is TRUE then PeList should not be NULL
  if(custom.list && is.null(PeList)){
    stop("When custom.list is TRUE, PeList must be provided.")
  }

  # Check if n is a positive integer
  if(!is.numeric(n) || n < 1 || (n %% 1 != 0)){
    stop("n must be a positive integer.")
  }

  # Check if n more than 2
  if(n <= 2){
    stop("n must be more than 2.")
  }

  system.time({
    if (!custom.list){PeList <- expand.grid(rep(list(Peptides::aaList()), n))}

    if (inherits(PeList, "data.frame") && !is.null(ncol(PeList))){
      PeList <- do.call(paste0, PeList[, 1:n])
    }

    if (rem.cys){PeList <- PeList[!grepl("C", PeList)]} # Remove sequences with Cys

    if (rem.met){PeList <- PeList[!grepl("M", PeList)]} # Remove sequences with Met

    if (rem.sali) {

      remove <- which(nchar(gsub("[AGILV]", "", PeList)) <= 2)  # Remove sequences with 2 or more small aliphatic amino acids

      if (length(remove) != 0) {
        PeList <- PeList[-remove]
      }
    }

    # Make sure PeList is not empty after all the modifications
    if(length(PeList) == 0){
      stop("After applying all filters, PeList is empty.")
    }

    peptides <- cbind(PeList,
                      as.data.frame(do.call(rbind, Peptides::crucianiProperties(PeList))),
                      as.data.frame(do.call(rbind, Peptides::fasgaiVectors(PeList))),
                      as.data.frame(do.call(rbind, Peptides::kideraFactors(PeList))),
                      as.data.frame(do.call(rbind, Peptides::protFP(PeList))),
                      as.data.frame(do.call(rbind, Peptides::tScales(PeList))),
                      as.data.frame(do.call(rbind, Peptides::vhseScales(PeList))),
                      as.data.frame(do.call(rbind, Peptides::zScales(PeList))))

    pI <- cbind(pI_EMBOS = Peptides::pI(PeList, pKscale = "EMBOSS"),
                pI_Bjellqvist = Peptides::pI(PeList, pKscale = "Bjellqvist"),
                pI_Lehninger = Peptides::pI(PeList, pKscale = "Lehninger"),
                pI_Murray = Peptides::pI(PeList, pKscale = "Murray"),
                pI_Rodwell = Peptides::pI(PeList, pKscale = "Rodwell"),
                pI_Sillero = Peptides::pI(PeList, pKscale = "Sillero"))

    mShft_15n <- Peptides::massShift(PeList, label = "15n", aaShift = NULL, monoisotopic = TRUE)

    charges <- cbind(Peptides::charge(PeList, pH = pH, pKscale = "EMBOSS"))

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
      hydrophobicity[[i]] <- Peptides::hydrophobicity (PeList, scale = hydrophobicity_indexes[i])
      names(hydrophobicity)[[i]] <- paste0("hb_",hydrophobicity_indexes[i])
    }


    peptides <- cbind(peptides, data.frame(aIndex = Peptides::aIndex(PeList),
                                           Boman = Peptides::boman(PeList),
                                           chrg_EMBOSS = charges,
                                           hmoment1 = Peptides::hmoment(PeList, angle = 100, window = 11),
                                           hmoment2 = Peptides::hmoment(PeList, angle = 160, window = 11),
                                           hydrophobicity,
                                           instaIndex = Peptides::instaIndex(PeList),
                                           mShft_15n,
                                           mw1 = Peptides::mw(PeList, monoisotopic = TRUE),
                                           pI))
    if (norm == T) {
      X <- peptides[,-1]
      Sequence <- peptides[,1]

      max = apply(X , 2 , max)
      min = apply(X, 2 , min)
      Xn = as.data.frame(scale(X, center = min, scale = max - min))

      peptides <- cbind(Sequence, Xn)
    }

    names(peptides)[1] <- "Sequence"

    return(peptides)
  })
}
