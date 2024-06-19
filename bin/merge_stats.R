#!/usr/bin/env Rscript

# setup ========================================================================
library(tidyverse)

# get arguments ================================================================
args <-  commandArgs(trailingOnly=TRUE)
sample_id <- args[1]

# define input files ===========================================================
stat_files <-  list.files(
  "input/", 
  pattern = ".stat.tsv",
  full.names = TRUE
)

# merge files ==================================================================
stat_files |> 
  map(read_tsv, col_names = c("Barcode", "n_reads")) |> 
  bind_rows() |> 
  group_by(Barcode) |> 
  summarise(n_reads = sum(n_reads)) |> 
  ungroup() |> 
  arrange(desc(n_reads)) |> 
  write_tsv(paste0(sample_id, "_merged.stat.tsv"))