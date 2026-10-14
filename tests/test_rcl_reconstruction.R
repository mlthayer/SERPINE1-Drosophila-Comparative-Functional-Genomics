# Run from the repository root: Rscript tests/test_rcl_reconstruction.R
library(Biostrings)
script <- 'Scripts/03_RCL_Alignment.R'
reference <- 'Sequence/RCL/Reference_Uploads'
work <- tempfile('rcl-reconstruction-')
dir.create(work)
old_output <- Sys.getenv('SERPINE1_RCL_OUTPUT', unset = NA_character_)
old_reference <- Sys.getenv('SERPINE1_RCL_REFERENCE_DIR', unset = NA_character_)
run <- function(args = character(), expected_failure = FALSE) {
  log <- file.path(work, 'run.log')
  status <- system2(file.path(R.home('bin'), 'Rscript'), c(shQuote(script), args),
                    stdout = log, stderr = log)
  if ((status != 0L) != expected_failure) stop(paste(readLines(log), collapse = '\n'))
  readLines(log)
}
compare <- function(output) {
  for (name in c('03_full_length_alignment_ClustalOmega.fasta',
                '03_SERPINE1_defined_homologous_RCL_alignment.fasta')) {
    a <- readAAStringSet(file.path(reference, name))
    b <- readAAStringSet(file.path(output, 'Processed', name))
    stopifnot(identical(names(a), names(b)), identical(as.character(a), as.character(b)))
  }
  for (name in c('03_full_length_alignment_table.csv',
                '03_SERPINE1_defined_homologous_RCL_regions.csv',
                '03_SERPINE1_referenced_RCL_region_scores.csv',
                '03_FASTA_inventory.csv', '03_full_length_sequence_inventory.csv',
                '03_Spn42DaB_reactive_site_reference.csv', '03_homologous_RCL_alignment_QC.csv',
                '03_human_RCL_sequence_validation.csv', '03_known_reactive_site_RCL_region_comparison.csv',
                '03_known_reactive_sites_full_alignment_QC.csv', '03_reactive_site_validation.csv',
                '03_SERPINE1_homologous_RCL_alignment_reference.csv')) {
    subdir <- if (name %in% c('03_full_length_alignment_table.csv',
                             '03_SERPINE1_defined_homologous_RCL_regions.csv',
                             '03_SERPINE1_referenced_RCL_region_scores.csv',
                             '03_known_reactive_site_RCL_region_comparison.csv')) 'Processed' else 'QC'
    a <- read.csv(file.path(reference, name), stringsAsFactors = FALSE)
    b <- read.csv(file.path(output, subdir, name), stringsAsFactors = FALSE)
    if ('Accession' %in% names(a)) {
      a <- a[order(a$Accession), ]; b <- b[order(b$Accession), ]
      rownames(a) <- rownames(b) <- NULL
    }
    stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-12)))
  }
  qc <- read.csv(file.path(output, 'QC', '03_reactive_site_alignment_QC.csv'))
  stopifnot(nrow(qc) == 6L, all(qc$Residues_Match), all(qc$Within_One_Alignment_Position))
}
tryCatch({
  Sys.setenv(SERPINE1_RCL_REFERENCE_DIR = reference,
             SERPINE1_RCL_OUTPUT = file.path(work, 'validated'))
  run(); compare(file.path(work, 'validated'))
  run(); compare(file.path(work, 'validated')) # repeatability
  Sys.setenv(SERPINE1_RCL_OUTPUT = file.path(work, 'realigned'))
  run('--realign'); compare(file.path(work, 'realigned'))
  bad <- file.path(work, 'bad-reference'); dir.create(bad)
  inputs <- list.files(reference, full.names = TRUE)
  stopifnot(all(file.copy(inputs, bad)))
  path <- file.path(bad, '03_SERPINE1_defined_homologous_RCL_regions.csv')
  table <- read.csv(path)
  table$Alignment_Interval_Start <- 1L
  write.csv(table, path, row.names = FALSE)
  Sys.setenv(SERPINE1_RCL_REFERENCE_DIR = bad, SERPINE1_RCL_OUTPUT = file.path(work, 'invalid'))
  result <- run(expected_failure = TRUE)
  stopifnot(any(grepl('does not match SERPINE1 residues', result)),
            !dir.exists(file.path(work, 'invalid', 'Processed')))
  cat('PASS: reference reproduction, repeatability, FASTA realignment, six supported-site checks, and invalid-boundary rejection.\n')
}, finally = {
  if (is.na(old_output)) Sys.unsetenv('SERPINE1_RCL_OUTPUT') else Sys.setenv(SERPINE1_RCL_OUTPUT = old_output)
  if (is.na(old_reference)) Sys.unsetenv('SERPINE1_RCL_REFERENCE_DIR') else Sys.setenv(SERPINE1_RCL_REFERENCE_DIR = old_reference)
  unlink(work, recursive = TRUE)
})
