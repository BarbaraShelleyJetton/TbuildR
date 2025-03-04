library(readr)
library(Biostrings)
library(DECIPHER)
library(rentrez)
library(microseq)
library(parallel)
library(stringr)


SeqSmash <- function(species_list, fasta_file, reference_file, output_file, replace_gaps = FALSE, gap_symbol = "?") {

  # Validate gap_symbol input
  if (!gap_symbol %in% c("?", "N")) {
    stop("Invalid gap_symbol. Choose '?' or 'N'.")
  }

  # Read species list
  species_list <- suppressMessages(read_csv(species_list))
  species_names <- species_list$species

  # Process the multi-FASTA file
  dat <- readDNAStringSet(fasta_file)
  all_seqs <- list()
  print("processing fasta files")

  for (i in seq_along(dat)) {
    header <- names(dat)[i]
    matches <- character(length(species_names))

    for (j in seq_along(species_names)) {
      pattern <- paste0("\\b", species_names[j], "\\b")
      match <- regmatches(header, regexpr(pattern, header, ignore.case = TRUE))
      if (length(match) > 0) {
        matches[j] <- match
      }
    }

    matches <- matches[matches != ""]

    if (length(matches) > 0) {
      species <- matches[1]

      # Format Names
      species_parts <- unlist(strsplit(species, " "))
      if (length(species_parts) == 2) {
        genus <- str_to_title(species_parts[1])
        species <- tolower(species_parts[2])
        species <- paste0(genus, "_", species)
      }

      gene <- sub(".*_(.*?)_.*", "\\1", header)

      if (!is.na(species) && !is.na(gene)) {
        if (!(species %in% names(all_seqs))) {
          all_seqs[[species]] <- list()
        }
        all_seqs[[species]][[gene]] <- dat[[i]]
      }
    }
  }

  # Remove empty elements
  all_seqs_filtered <- all_seqs[!sapply(all_seqs, function(x) all(sapply(x, is.null)))]

  # Read the reference genome
  ref_mito <- readDNAStringSet(ref_mito_file)

  # Temporary file for consensus sequences
  temp_fasta <- tempfile(fileext = ".fasta")
  file_conn <- file(temp_fasta, "w")

  # Generate consensus sequences and write to the temporary FASTA file
  for (species_name in names(all_seqs_filtered)) {
    species_data <- all_seqs_filtered[[species_name]]
    gene_seqs <- list()

    for (gene_name in names(species_data)) {
      if (!is.null(species_data[[gene_name]])) {
        gene_seqs[[gene_name]] <- species_data[[gene_name]]
      }
    }

    # Combine extracted gene sequences with the reference genome
    seqs <- c(DNAStringSet(unlist(gene_seqs)), ref_mito)
    aln <- AlignSeqs(seqs, processors = (detectCores() - 1))

    # Generate the consensus sequence
    cons_seq <- ConsensusSequence(aln[-length(aln)])

    bstring <- BStringSet(cons_seq)

    # Gap replacement
    if (replace_gaps) {
      cons_seq_mod <- gsub("-", gap_symbol, as.character(bstring))
    } else {
      cons_seq_mod <- as.character(bstring)  # Keep original sequence
    }

    line <- paste0(">", species_name, "\n", cons_seq_mod, "\n")
    writeLines(line, file_conn)
  }

  close(file_conn)
  file.copy(temp_fasta, output_file, overwrite = TRUE)

  cat("Cleaned FASTA file saved to:", output_file, "\n")
}
