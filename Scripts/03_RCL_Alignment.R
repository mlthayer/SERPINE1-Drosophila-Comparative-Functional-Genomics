# SERPINE1–Drosophila Comparative Functional Genomics Workflow
# Reconstructed homologous RCL extraction from validated alignment/table inputs.
# The original motif-based workflow is retained under Scripts/Legacy/.
# Run from the repository root. --realign regenerates the full-length MSA.
library(Biostrings)
library(msa)

reference_dir <- Sys.getenv("SERPINE1_RCL_REFERENCE_DIR", unset = "Sequence/RCL/Reference_Uploads")
output_dir <- Sys.getenv("SERPINE1_RCL_OUTPUT", unset = "Sequence/RCL")
region_input <- file.path(reference_dir, "03_SERPINE1_defined_homologous_RCL_regions.csv")
alignment_input <- file.path(reference_dir, "03_full_length_alignment_ClustalOmega.fasta")
fasta_input <- Sys.getenv("SERPINE1_RCL_FASTA_INPUT", unset = "Sequence/RCL/selected_human_fly_Serpins_Fasta.txt")
metadata <- read.csv(region_input, stringsAsFactors = FALSE)
required <- c("Gene", "Accession", "Display_Name", "Species", "Region_Definition",
              "Reactive_Site_Assigned", "Alignment_Interval_Start", "Alignment_Interval_End")
if (!all(required %in% names(metadata)) || anyDuplicated(metadata$Gene) ||
    anyDuplicated(metadata$Accession) || !nrow(metadata)) {
  stop("Region table requires unique genes/accessions and the documented metadata columns.")
}
validated_alignment <- readAAStringSet(alignment_input)
if (!setequal(names(validated_alignment), metadata$Gene) ||
    length(unique(width(validated_alignment))) != 1L) stop("Alignment/table record mismatch or unequal widths.")
validated_alignment <- validated_alignment[metadata$Gene]
raw <- readAAStringSet(fasta_input)
accessions <- vapply(strsplit(names(raw), "|", fixed = TRUE), function(x) {
  if (length(x) < 2L) stop("Original FASTA headers must contain UniProt accessions.")
  x[2]
}, character(1))
if (anyDuplicated(accessions) || !setequal(accessions, metadata$Accession)) stop("Original FASTA accession mismatch.")
raw <- raw[match(metadata$Accession, accessions)]
original_headers <- names(raw)
names(raw) <- metadata$Gene
if (!identical(unname(gsub("-", "", as.character(validated_alignment), fixed = TRUE)),
               unname(as.character(raw)))) stop("Validated alignment differs from original FASTA sequences.")

# SERPINE1 P17–P8' is explicitly defined as residues 353–377.
# Audit the validated table against that biological interval.
reference_gene <- metadata$Gene[metadata$Accession == "P05121"]
if (length(reference_gene) != 1L) stop("Exactly one SERPINE1 P05121 reference is required.")
start_column <- unique(metadata$Alignment_Interval_Start)
end_column <- unique(metadata$Alignment_Interval_End)
if (length(start_column) != 1L || length(end_column) != 1L ||
    !is.finite(start_column) || !is.finite(end_column) ||
    start_column != floor(start_column) || end_column != floor(end_column) ||
    start_column < 1L || end_column < start_column || end_column > width(validated_alignment)[1]) {
  stop("The region table must define one valid common alignment interval.")
}
reference_chars <- strsplit(as.character(validated_alignment[[reference_gene]]), "", fixed = TRUE)[[1]]
reference_positions <- which(reference_chars != "-")
reference_start <- 353L
reference_end <- 377L
if (!identical(as.integer(match(start_column, reference_positions)), reference_start) ||
    !identical(as.integer(match(end_column, reference_positions)), reference_end)) {
  stop("Validated interval does not match SERPINE1 residues 353–377.")
}
if (substr(as.character(raw[[reference_gene]]), 369L, 370L) != "RM") {
  stop("SERPINE1 reactive bond must be R369–M370 in the supplied full-length FASTA.")
}

aligned <- validated_alignment
if ("--realign" %in% commandArgs(trailingOnly = TRUE)) {
  aligned <- as(msa::msa(raw, method = "ClustalOmega", order = "input"), "AAStringSet")
  aligned <- aligned[metadata$Gene]
}
reference_chars <- strsplit(as.character(aligned[[reference_gene]]), "", fixed = TRUE)[[1]]
reference_positions <- which(reference_chars != "-")
start_column <- reference_positions[reference_start]
end_column <- reference_positions[reference_end]
# Supplied site annotations are used only to check alignment consistency.
# Coordinates/residues are preserved; sequence agreement is not experimental
# validation. See Sequence/RCL/Citation_Audit for primary-source audit status.
supported_sites <- data.frame(
  Accession = c("P05121", "P29508", "P48594", "Q96P15", "Q99574", "A1Z6R1"),
  P1_Position = c(369L, 354L, 354L, 357L, 362L, 374L),
  P1prime_Position = c(370L, 355L, 355L, 358L, 363L, 375L),
  P1 = c("R", "S", "L", "K", "R", "R"),
  P1prime = c("M", "S", "S", "S", "M", "A"), stringsAsFactors = FALSE)
supported_sites$Gene <- metadata$Gene[match(supported_sites$Accession, metadata$Accession)]
if (anyNA(supported_sites$Gene)) stop("Independently supported site records are missing.")
for (i in seq_len(nrow(supported_sites))) {
  site <- supported_sites[i, ]
  chars <- strsplit(as.character(aligned[[site$Gene]]), "", fixed = TRUE)[[1]]
  positions <- which(chars != "-")
  if (site$P1prime_Position > length(positions)) stop("Supported reactive site is outside its protein.")
  supported_sites$P1_Alignment_Column[i] <- positions[site$P1_Position]
  supported_sites$P1prime_Alignment_Column[i] <- positions[site$P1prime_Position]
  supported_sites$Observed_P1[i] <- chars[positions[site$P1_Position]]
  supported_sites$Observed_P1prime[i] <- chars[positions[site$P1prime_Position]]
}
supported_sites$P1_Offset_From_SERPINE1 <- supported_sites$P1_Alignment_Column - reference_positions[369L]
supported_sites$P1prime_Offset_From_SERPINE1 <- supported_sites$P1prime_Alignment_Column - reference_positions[370L]
supported_sites$Residues_Match <- supported_sites$Observed_P1 == supported_sites$P1 &
  supported_sites$Observed_P1prime == supported_sites$P1prime
supported_sites$Within_One_Alignment_Position <-
  abs(supported_sites$P1_Offset_From_SERPINE1) <= 1L &
  abs(supported_sites$P1prime_Offset_From_SERPINE1) <= 1L
if (!all(supported_sites$Residues_Match & supported_sites$Within_One_Alignment_Position)) {
  stop("Independently supported reactive-site alignment QC failed.")
}

regions <- data.frame(
  Gene = metadata$Gene,
  Homologous_RCL_Alignment = substring(as.character(aligned), start_column, end_column),
  Alignment_Interval_Start = start_column,
  Alignment_Interval_End = end_column,
  Alignment_Interval_Length = end_column - start_column + 1L,
  stringsAsFactors = FALSE
)
regions$Ungapped_Region <- gsub("-", "", regions$Homologous_RCL_Alignment, fixed = TRUE)
regions$Ungapped_Length <- nchar(regions$Ungapped_Region)
for (column in c("Accession", "Display_Name", "Species", "Region_Definition", "Reactive_Site_Assigned")) {
  regions[[column]] <- metadata[[column]]
}
# Reactive-site annotations are retained metadata, not newly inferred sites.
data("BLOSUM62", package = "Biostrings")
reference <- strsplit(regions$Homologous_RCL_Alignment[regions$Gene == reference_gene], "", fixed = TRUE)[[1]]
used <- reference != "-"
if (!any(used)) stop("SERPINE1 region has no reference residues.")
scores <- do.call(rbind, lapply(seq_len(nrow(regions)), function(i) {
  candidate <- strsplit(regions$Homologous_RCL_Alignment[i], "", fixed = TRUE)[[1]]
  paired <- used & candidate != "-"
  if (any(!candidate[paired] %in% rownames(BLOSUM62))) stop("Unsupported residue in BLOSUM62 comparison.")
  values <- BLOSUM62[cbind(reference[paired], candidate[paired])]
  identical_residues <- sum(candidate[used] == reference[used])
  similar_residues <- sum(candidate[paired] == reference[paired] | values > 0)
  data.frame(
    Gene = regions$Gene[i], Reference_Residues = sum(used),
    Identical_Residues = identical_residues,
    RCL_Region_Identity_Percent = 100 * identical_residues / sum(used),
    Similar_Residues_BLOSUM62_Positive = similar_residues,
    RCL_Region_Similarity_Percent = 100 * similar_residues / sum(used),
    Candidate_Gaps_at_SERPINE1_Residues = sum(used & candidate == "-"),
    Candidate_Residues_at_SERPINE1_Gap_Columns = sum(!used & candidate != "-"),
    Mean_BLOSUM62_NonGap = if (length(values)) mean(values) else NA_real_,
    Display_Name = regions$Display_Name[i], Accession = regions$Accession[i],
    Species = regions$Species[i], Ungapped_Length = regions$Ungapped_Length[i],
    Reactive_Site_Assigned = regions$Reactive_Site_Assigned[i],
    Region_Definition = regions$Region_Definition[i], stringsAsFactors = FALSE
  )
}))
processed_dir <- file.path(output_dir, "Processed")
qc_dir <- file.path(output_dir, "QC")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
writeXStringSet(aligned, file.path(processed_dir, "03_full_length_alignment_ClustalOmega.fasta"))
write.csv(data.frame(Gene = names(aligned), Aligned_Full_Length = as.character(aligned),
                     Alignment_Length = width(aligned)),
          file.path(processed_dir, "03_full_length_alignment_table.csv"), row.names = FALSE)
rcl_set <- AAStringSet(setNames(regions$Homologous_RCL_Alignment, regions$Gene))
writeXStringSet(rcl_set, file.path(processed_dir, "03_SERPINE1_defined_homologous_RCL_alignment.fasta"))
write.csv(regions, file.path(processed_dir, "03_SERPINE1_defined_homologous_RCL_regions.csv"), row.names = FALSE)
write.csv(scores, file.path(processed_dir, "03_SERPINE1_referenced_RCL_region_scores.csv"), row.names = FALSE)
inventory <- data.frame(Original_Header = original_headers, Accession = metadata$Accession,
                        Protein_Length = width(raw), Gene = metadata$Gene,
                        Display_Name = metadata$Display_Name, Species = metadata$Species)
write.csv(inventory, file.path(qc_dir, "03_FASTA_inventory.csv"), row.names = FALSE)
write.csv(inventory[c("Gene", "Display_Name", "Accession", "Species", "Protein_Length")],
          file.path(qc_dir, "03_full_length_sequence_inventory.csv"), row.names = FALSE)
write.csv(data.frame(Reference_Accession = "P05121", Reference_Start = reference_start,
                     Reference_End = reference_end, Alignment_Start = start_column,
                     Alignment_End = end_column, Reference_Residues = sum(used),
                     Reconstruction = "Explicit author-defined SERPINE1 P17-P8prime residues 353-377; validated-table boundary audit"),
          file.path(qc_dir, "03_Reconstruction_Provenance.csv"), row.names = FALSE)
write.csv(supported_sites, file.path(qc_dir, "03_reactive_site_alignment_QC.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "03_sessionInfo.txt"))

# Recreate the supplied QC records independently from mapped coordinates.
supported_sites$Evidence_Class <- ifelse(supported_sites$Accession == "A1Z6R1",
  "experimentally_characterized_literature_R342_A343_full_FASTA_R374_A375", "validated")
site_validation <- supported_sites[c("Gene", "P1_Position", "P1", "P1prime_Position", "P1prime", "Evidence_Class", "Observed_P1", "Observed_P1prime")]
site_validation$P1_OK <- supported_sites$Observed_P1 == supported_sites$P1
site_validation$P1prime_OK <- supported_sites$Observed_P1prime == supported_sites$P1prime
site_validation$Reactive_Site_OK <- site_validation$P1_OK & site_validation$P1prime_OK
write.csv(site_validation, file.path(qc_dir, "03_reactive_site_validation.csv"), row.names = FALSE)
site_metadata <- metadata[match(supported_sites$Accession, metadata$Accession), ]
full_site_qc <- data.frame(
  Gene = supported_sites$Gene, Display_Name = site_metadata$Display_Name,
  Accession = supported_sites$Accession, Species = site_metadata$Species,
  P1_Position = supported_sites$P1_Position, P1 = supported_sites$P1,
  P1prime_Position = supported_sites$P1prime_Position, P1prime = supported_sites$P1prime,
  P1_Alignment_Column = supported_sites$P1_Alignment_Column,
  P1prime_Alignment_Column = supported_sites$P1prime_Alignment_Column,
  Alignment_Bond_Position = paste(supported_sites$P1_Alignment_Column, supported_sites$P1prime_Alignment_Column, sep = "|"),
  P1_Offset_From_SERPINE1 = supported_sites$P1_Offset_From_SERPINE1,
  P1prime_Offset_From_SERPINE1 = supported_sites$P1prime_Offset_From_SERPINE1,
  Evidence_Class = supported_sites$Evidence_Class, stringsAsFactors = FALSE)
write.csv(full_site_qc, file.path(qc_dir, "03_known_reactive_sites_full_alignment_QC.csv"), row.names = FALSE)
homologous_qc <- regions[c("Gene", "Display_Name", "Accession", "Species", "Alignment_Interval_Length")]
homologous_qc$Number_of_Gaps <- nchar(regions$Homologous_RCL_Alignment) - regions$Ungapped_Length
for (column in c("Ungapped_Length", "Reactive_Site_Assigned", "Region_Definition")) homologous_qc[[column]] <- regions[[column]]
write.csv(homologous_qc, file.path(qc_dir, "03_homologous_RCL_alignment_QC.csv"), row.names = FALSE)
write.csv(scores[scores$Accession %in% supported_sites$Accession, ],
          file.path(processed_dir, "03_known_reactive_site_RCL_region_comparison.csv"), row.names = FALSE)
human_sites <- supported_sites[supported_sites$Accession != "A1Z6R1", ]
human_validation <- data.frame(
  Gene = human_sites$Gene, RCL_Start = human_sites$P1_Position - 16L,
  RCL_End = human_sites$P1prime_Position + 7L, RCL_Length = 25L,
  RCL_Sequence = vapply(seq_len(nrow(human_sites)), function(i) {
    substr(as.character(raw[[human_sites$Gene[i]]]), human_sites$P1_Position[i] - 16L,
           human_sites$P1prime_Position[i] + 7L)
  }, character(1)), stringsAsFactors = FALSE)
for (column in c("P1_Position", "P1", "P1prime_Position", "P1prime", "Evidence_Class")) human_validation[[column]] <- human_sites[[column]]
human_validation$RCL_Definition <- "P17-P8prime"
human_validation$Site_Status <- "validated"
# Expected human sequences are retained from the supplied validation records.
expected_human_rcls <- c(
  SERPINE1 = "ESGTVASSSTAVIVSARMAPEEIIM", SERPINB3 = "EEGAEAAAATAVVGFGSSPTSTNEE",
  SERPINB4 = "EEGVEAAAATAVVVVELSSPSTNEE", SERPINB11 = "EEGTEAAAATGDSIAVKSLPMRAQF",
  SERPINI1 = "EEGSEAAAVSGMIAISRMAVLYPQV")
human_validation$Expected_RCL <- unname(expected_human_rcls[human_validation$Gene])
human_validation$Sequence_OK <- human_validation$RCL_Sequence == human_validation$Expected_RCL
if (!all(human_validation$Sequence_OK)) stop("Known human RCL sequence validation failed.")
write.csv(human_validation, file.path(qc_dir, "03_human_RCL_sequence_validation.csv"), row.names = FALSE)
ref_site <- supported_sites[supported_sites$Accession == "P05121", ]
reference_qc <- data.frame(
  Reference_Gene = reference_gene, Reference_Window = "P17-P8prime",
  Full_Length_RCL_Start = reference_start, Full_Length_P1 = ref_site$P1_Position,
  P1_Residue = ref_site$P1, Full_Length_P1prime = ref_site$P1prime_Position,
  P1prime_Residue = ref_site$P1prime, Full_Length_RCL_End = reference_end,
  Full_Alignment_Start_Column = start_column, Full_Alignment_P1_Column = ref_site$P1_Alignment_Column,
  Full_Alignment_P1prime_Column = ref_site$P1prime_Alignment_Column,
  Full_Alignment_End_Column = end_column,
  Local_P1_Column = ref_site$P1_Alignment_Column - start_column + 1L,
  Local_P1prime_Column = ref_site$P1prime_Alignment_Column - start_column + 1L)
write.csv(reference_qc, file.path(qc_dir, "03_SERPINE1_homologous_RCL_alignment_reference.csv"), row.names = FALSE)
fly_site <- supported_sites[supported_sites$Accession == "A1Z6R1", ]
fly_reference <- data.frame(
  Gene = fly_site$Gene, Accession = fly_site$Accession,
  Literature_P1_Position = 342L, Literature_P1prime_Position = 343L,
  FASTA_P1_Position = fly_site$P1_Position, P1_Residue = fly_site$P1,
  FASTA_P1prime_Position = fly_site$P1prime_Position, P1prime_Residue = fly_site$P1prime,
  Full_Alignment_P1_Column = fly_site$P1_Alignment_Column,
  Full_Alignment_P1prime_Column = fly_site$P1prime_Alignment_Column,
  Evidence = "experimentally_structurally_supported_R-A_reactive_bond")
write.csv(fly_reference, file.path(qc_dir, "03_Spn42DaB_reactive_site_reference.csv"), row.names = FALSE)

message("Completed reconstructed homologous RCL extraction: ", nrow(regions),
        " proteins; ", sum(used), " reference residues across ", length(reference), " alignment columns.")
