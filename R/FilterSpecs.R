#' FilterSpecs
#'
#' Filters a species list based on sequence matches in FASTA file. Optionally, returns species with no matches ("unmatched")
#' or those meeting a minimum hit count ("filtered"). Can also write a filtered FASTA containing only matched species.
#'
#' @param fasta_file FASTA file to be filtered
#' @param species_list csv file containing a "species" column
#' @param min_hits minimum number of matches required for a species to be retained
#' @param output_type Either "unmatched" (to return species with 0 hits) or "filtered" (to return species meeting hit threshold)
#' @param output_file output file path / name for resulting species csv
#' @param output_fasta output FASTA path / name for filtered FASTA
#'
#' @return CSV and FASTA
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr filter select left_join %>%
#' @importFrom stats setNames
#'
#' @export
#'
#' @examples
#' \dontrun{
#' FilterSpecs(
#' fasta_file = "./path.fasta",
#' species_list = "./path.csv",
#' min_hits = 3,
#' output_type = "filtered",
#' output_file = "./path.csv",
#' output_fasta = "./path.fasta"
#' )
#' }

FilterSpecs <- function(fasta_file, species_list, min_hits = 1, output_type = c("unmatched", "filtered"), output_file = "species_output.csv", output_fasta = NULL
) {
  output_type <- match.arg(output_type)

  # Read full species list with extra metadata
  species_df <- read_csv(species_list, show_col_types = FALSE)
  if (!"species" %in% colnames(species_df)) stop("CSV must contain a 'species' column.")

  species_names <- species_df$species

  # Build regex patterns
  species_patterns <- setNames(paste0("\\b", species_names, "\\b"), species_names)
  species_counts <- integer(length(species_names))
  names(species_counts) <- species_names

  for (file_path in fasta_file) {
    cat("Processing:", file_path, "\n")
    seqs <- readDNAStringSet(file_path)
    headers <- names(seqs)

    for (species in names(species_patterns)) {
      pattern <- species_patterns[[species]]
      matches <- grepl(pattern, headers, ignore.case = TRUE)
      species_counts[species] <- species_counts[species] + sum(matches)
    }
  }

  # Create summary table and join to original CSV
  count_df <- data.frame(species = names(species_counts), count = species_counts)
  joined_df <- left_join(species_df, count_df, by = "species")

  # Output based on desired filter
  if (output_type == "unmatched") {
    result_df <- joined_df %>% filter(count == 0 | is.na(count)) %>% select(-count)
    write_csv(result_df, output_file)
    cat("Wrote unmatched species (including metadata) to:", output_file, "\n")

  } else if (output_type == "filtered") {
    result_df <- joined_df %>% filter(count >= min_hits) %>% select(-count)
    write_csv(result_df, output_file)
    cat("Wrote filtered species list (including metadata) to:", output_file, "\n")

    if (!is.null(output_fasta)) {
      keep_species <- result_df$species
      keep_pattern <- paste0("\\b(", paste(keep_species, collapse = "|"), ")\\b")

      combined_seqs <- DNAStringSet()
      for (file_path in fasta_file) {
        seqs <- readDNAStringSet(file_path)
        headers <- names(seqs)
        matches <- grepl(keep_pattern, headers, ignore.case = TRUE)
        combined_seqs <- c(combined_seqs, seqs[matches])
      }

      writeXStringSet(combined_seqs, output_fasta)
      cat("Wrote filtered FASTA to:", output_fasta, "\n")
    }
  }
}
