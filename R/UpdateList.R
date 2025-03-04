library(Biostrings)
library(dplyr)

UpdateList <- function(species_list, fasta_file, output_file) {
  # Read the species list
  species_list <- read.csv(species_list, stringsAsFactors = FALSE)

  # Ensure the correct column name exists
  if (!"species" %in% colnames(species_list)) {
    stop("The species list file does not have a 'species' column.")
  }

  species_names <- species_list$species

  # Read the multi-FASTA file
  fasta_data <- readDNAStringSet(fasta_file)
  fasta_headers <- names(fasta_data)

  # Extract species names from FASTA headers by matching the species list
  species_found <- sapply(species_names, function(species) {
    any(grepl(species, fasta_headers, ignore.case = TRUE))
  })

  # Filter species not found in the FASTA file
  species_remaining <- species_list[!species_found, , drop = FALSE]

  # Ensure column names are preserved when writing output
  if (nrow(species_remaining) == 0) {
    message("No species were removed. Writing an empty file with headers.")
    write.csv(species_list[0, , drop = FALSE], file = output_file, row.names = FALSE, quote = FALSE)
  } else {
    write.csv(species_remaining, file = output_file, row.names = FALSE, quote = FALSE)
  }

  cat("Updated species list saved to:", output_file, "\n")
}
