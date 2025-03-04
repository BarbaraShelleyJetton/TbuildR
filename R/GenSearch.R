library(readr)
library(Biostrings)
library(DECIPHER)
library(rentrez)
library(microseq)

GenSearch <- function(species_list, gene_list, output_file, retmax_value = 10, desired_length = 18000, search_field = "GENE") {

  # Load the species list from CSV
  species_list <- suppressMessages(read_csv(species_list))

  # Checks if gene_list is a file path or a character vector
  if (is.character(gene_list) && length(gene_list) == 1 && file.exists(gene_list)) {
    # If it's a single string and a valid file path, read the gene list from the CSV file
    gene_list_df <- suppressMessages(read_csv(gene_list))
    if (!"gene" %in% names(gene_list_df)) {
      stop("The gene list CSV file must contain a 'gene' column.")
    }
    gene_list <- gene_list_df$gene
  }

  # Ensure gene_list is a character vector
  if (!is.character(gene_list)) {
    stop("gene_list must be either a character vector or a CSV file containing a 'gene' column.")
  }

  # Generate and perform search on GenBank
  generate_and_perform_search <- function(species_list, gene_list, retmax_value) {
    search_results <- list()
    for (species in species_list$species) {
      for (gene in gene_list) {
        search_term <- paste0(species, "[Organism] AND ", gene, "[", search_field, "]")
        search_result <- entrez_search(db = "nucleotide", term = search_term, retmax = retmax_value, use_history = TRUE)

        if (search_result$count > 0) {
          search_results[[paste(species, gene)]] <- search_result
        } else {
          cat("No hits found for: ", search_term, "\n")
        }
      }
    }
    return(search_results)
  }

  # Get Summaries
  get_summaries <- function(search_results) {
    cat("Pulling summaries...\n")
    summaries <- list()
    for (key in names(search_results)) {
      ids <- search_results[[key]]$ids
      if (length(ids) > 0) {
        sum_result <- tryCatch(
          entrez_summary(db = "nucleotide", id = ids),
          error = function(e) {
            cat("Error fetching summary for", key, ": ", e$message, "\n")
            return(NULL)
          }
        )
        if (!is.null(sum_result)) {
          summaries[[key]] <- sum_result
        } else {
          cat("No esummary records found for:", key, "\n")
        }
      } else {
        cat("No IDs found for:", key, "\n")
      }
    }
    return(summaries)
  }

  # Filter accessions
  get_filtered_accessions <- function(species_summaries, gene_list, desired_length) {
    filtered_accessions <- list()
    for (key in names(species_summaries)) {
      summaries <- species_summaries[[key]]
      filtered_accession <- NULL
      if (is.list(summaries) && length(summaries) == 33) {
        if (!is.null(summaries$slen) && !is.null(summaries$accession) && !is.null(summaries$title)) {
          contains_gene <- any(sapply(gene_list, function(gene) grepl(gene, summaries$title, ignore.case = TRUE)))
          if (contains_gene && (is.null(desired_length) || summaries$slen <= desired_length)) {
            filtered_accession <- summaries$accession
          }
        }
      } else if (is.list(summaries)) {
        for (i in seq_along(summaries)) {
          current_hit <- summaries[[i]]
          current_length <- current_hit$slen
          current_title <- current_hit$title
          contains_gene <- any(sapply(gene_list, function(gene) grepl(gene, current_title, ignore.case = TRUE)))
          if (contains_gene && (is.null(desired_length) || current_length <= desired_length)) {
            filtered_accession <- current_hit$accession
            break
          }
        }
      }
      filtered_accessions[[key]] <- filtered_accession
    }
    return(filtered_accessions)
  }

  # Fetch and write Fasta file
  get_fasta <- function(accession_numbers, output_file) {
    cat("Writing to file...\n")
    all_recs <- character()
    for (accession_list in accession_numbers) {
      for (accession in accession_list) {
        rec <- entrez_fetch(db = "nuccore", id = accession, rettype = "fasta")
        all_recs <- c(all_recs, rec)
      }
    }
    writeLines(all_recs, con = output_file)
  }

  # Execute
  search_results <- generate_and_perform_search(species_list, gene_list, retmax_value)
  species_summaries <- get_summaries(search_results)
  filtered_accessions <- get_filtered_accessions(species_summaries, gene_list, desired_length)
  get_fasta(filtered_accessions, output_file)
}
