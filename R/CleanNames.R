#' CleanNames
#'
#' Cleans FASTA headers by matching species names to a reference list and renaming them in a standardized format.
#'
#' @param species_list csv file containing a "species" column
#' @param fasta_file FASTA file to be cleaned
#' @param output_file output file path / name for the cleaned FASTA
#' @param remove_unmatched Logical. If TRUE, removes sequences that don't match any species in the list.
#'
#' @return FASTA file
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom stringr str_to_title
#' @importFrom readr read_csv
#'
#' @return FASTA
#'
#' @export
#'
#' @examples
#' \dontrun{
#' CleanNames(
#' species_list = "./path.csv",
#' fasta_file = "./path.fasta",
#' output_file = "./path.fasta",
#' remove_unmatched = FALSE
#' )
#' }

CleanNames <- function(species_list, fasta_file, output_file = NULL, remove_unmatched = FALSE) {
  # Read the species list using readr
  species_df <- suppressMessages(read_csv(species_list))
  species_names <- species_df$species

  # Read the multi-FASTA file
  fasta_data <- readDNAStringSet(fasta_file)
  fasta_headers <- names(fasta_data)

  # Function to format header if species name is matched
  get_species_name <- function(header) {
    matched_species <- species_names[sapply(species_names, function(species) grepl(species, header, ignore.case = TRUE))]
    if (length(matched_species) > 0) {
      species <- matched_species[1]
      parts <- unlist(strsplit(species, " "))
      if (length(parts) == 2) {
        genus <- str_to_title(parts[1])
        sp <- tolower(parts[2])
        return(paste0(genus, "_", sp))
      }
    }
    return(NA)
  }

  # Update headers
  updated_headers <- sapply(fasta_headers, get_species_name)

  if (!remove_unmatched) {
    updated_headers[is.na(updated_headers)] <- fasta_headers[is.na(updated_headers)]
    names(fasta_data) <- updated_headers
  } else {
    matched_indices <- which(!is.na(updated_headers))
    fasta_data <- fasta_data[matched_indices]
    names(fasta_data) <- updated_headers[matched_indices]
  }

  if (is.null(output_file)) {
    output_file <- sub("\\.fasta$", "_cleaned.fasta", fasta_file)
  }

  writeXStringSet(fasta_data, filepath = output_file, format = "fasta")
  cat("Cleaned FASTA saved to:", output_file, "\n")
}
