library(Biostrings)
library(dplyr)
library(stringr)

CleanNames <- function(species_list, fasta_file, output_file = NULL) {
  # Read the species list
  species_list <- read.csv(species_list, stringsAsFactors = FALSE)
  species_names <- species_list$species

  # Read the multi-FASTA file
  fasta_data <- readDNAStringSet(fasta_file)
  fasta_headers <- names(fasta_data)

  # Function to find and format species name
  get_species_name <- function(header) {
    matched_species <- species_names[sapply(species_names, function(species) grepl(species, header, ignore.case = TRUE))]

    if (length(matched_species) > 0) {
      # Get first matched species
      species <- matched_species[1]

      # Split into genus and species
      species_parts <- unlist(strsplit(species, " "))

      if (length(species_parts) == 2) {
        genus <- str_to_title(species_parts[1])
        species <- tolower(species_parts[2])
        formatted_species <- paste0(genus, "_", species)
        return(formatted_species)
      }
    }

    # If no match found, return the original header unchanged
    return(header)
  }

  # Apply function to all headers
  updated_headers <- sapply(fasta_headers, get_species_name)

  # Keep all sequences, updating names where matches exist
  names(fasta_data) <- updated_headers

  # Set output file name if not provided
  if (is.null(output_file)) {
    output_file <- sub("\\.fasta$", "_cleaned.fasta", fasta_file)
  }

  writeXStringSet(fasta_data, filepath = output_file, format = "fasta")

  cat("Updated FASTA file saved to:", output_file, "\n")
}
