#!/usr/bin/env Rscript
#
# ------------------------------------------------
# Author: Michael Jahn, PhD
# Affiliation:
#   Science for Life Lab (KTH), Stockholm, Sweden, and
#   Max Planck Unit for the Science of Pathogens
#   (MPUSP - MPG), Berlin, Germany
# Email: jahn@mpusp.mpg.de
# Website: https://github.com/m-jahn/TnSeq-pipe
# Date: last change 2022-0928
# Description: This R script downloads and modifies
# reference genome files from NCBI (gff, fasta, genbank)
# ------------------------------------------------
#
# NOTE:
# List of NCBI databases and corresponding IDs to query and fetch data:
# https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/
#
#
# Load libraries
req_packages <- c("tidyverse", "stringi", "rentrez", "RCurl")
if (any(!req_packages %in% rownames(installed.packages()))) {
    stop(paste0("Not all required R packages (",
        paste(req_packages, collapse = ", "),
        ") are available, please install them."))
} else {
    message(paste0("Required R packages ",
        paste(req_packages, collapse = ", "), " are installed"))
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(stringi)
    library(rentrez)
    library(RCurl)
})

# get command line argument (genome assembly ID)
assembly <- commandArgs(trailingOnly = TRUE)

message(paste("Searching NCBI for supplied assembly ID:", assembly))
search_result <- entrez_search(db = "assembly", term = assembly)
search_summary <- entrez_summary(db = "assembly", id = search_result$ids)
search_link <- entrez_link(
    db = "nuccore", dbfrom = "assembly",
    id = search_result$ids
)
message(paste0(
    "Found organism '", search_summary$organism,
    "' matching assembly ID"
))

result_filename <- paste(
    search_summary$assemblyaccession,
    search_summary$assemblyname, "genomic",
    sep = "_"
)

message("Fetching Genbank file from NCBI")
result_gbk <- entrez_fetch(
    db = "nuccore",
    id = search_link$links[["assembly_nuccore_refseq"]], rettype = "gb"
)

message("Fetching FASTA file from NCBI")
result_fna <- entrez_fetch(
    db = "nuccore",
    id = search_link$links[["assembly_nuccore_refseq"]], rettype = "fasta"
)

message("Fetching GFF file from NCBI")
result_ftp_path <- paste0(
    gsub("^ftp:", "https:", search_summary$ftppath_refseq), "/",
    result_filename, ".gff.gz"
)

if (url.exists(result_ftp_path)) {
    result_gff <- read_tsv(result_ftp_path,
        comment = "#",
        col_names = FALSE,
        show_col_types = FALSE
    )
} else {
    stop(paste0(
        "The FTP server at '", result_ftp_path,
        "' is not available. Check your internet connection"
    ))
}

result_gff <- result_gff %>%
    filter(X3 != "CDS", X3 != "region") %>%
    select(X1, X4, X5, X7, X3, X9) %>%
    rename_with(~ c(
        "scaffold", "begin", "end", "strand", "desc",
        "gene_ID"
    ))

tag_terms <- str_extract_all(
    result_gff[1, "gene_ID"], ";[a-z_]*locus_tag="
)[[1]]
for (term in tag_terms) {
    new_col <- str_remove_all(term, "[;=]")
    result_gff <- result_gff %>%
        mutate({{ new_col }} := str_extract(
            gene_ID, paste0(term, "[a-zA-Z0-9_]*")
        ) %>%
            str_remove(term))
}
result_gff <- select(result_gff, -gene_ID)

message(paste0(
    "Writing processed files ", result_filename,
    ".gbk (.gff, .fna) to ref/"
))
write_file(result_gbk, paste0("./ref/", result_filename, ".gbk"))
write_file(result_fna, paste0("./ref/", result_filename, ".fna"))
write_tsv(result_gff, paste0("./ref/", result_filename, ".gff"))
