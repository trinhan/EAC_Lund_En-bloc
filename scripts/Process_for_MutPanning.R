# Function to process selected variants for mutpanning
process_variants_mutpanning <- function(data, include_hgvsp = FALSE) {
  # Select columns
  selected_cols <- c("Hugo_Symbol", "CHROM", "POS", "STRAND", "Consequence",
                     "REF", "ALT", "Tumor_Sample_Barcode")
  
  if (include_hgvsp) {
    selected_cols <- c(selected_cols, "HGVSp_Short")
  }
  
  output <- data[, selected_cols]
  
  # Calculate variant properties
  result <- calculate_variant_properties(output)
  output <- result$data
  ref_lengths <- result$ref_lengths
  alt_lengths <- result$alt_lengths
  
  # Map consequences
  output <- map_consequence_annotations(output, ref_lengths, alt_lengths)
  
  return(output)
}

# Function to map consequence annotations for mutpanning
map_consequence_annotations <- function(data, ref_lengths, alt_lengths) {
  # Define consequence mappings
  consequence_patterns <- c(
    "stop_gained", "frameshift", "missense", "3_prime", "5_prime", "splice",
    "inframe_ins", "inframe_del", "synonymous", "intergenic", "intron",
    "non_coding_transcript", "stop_retained"
  )
  
  consequence_values <- c(
    "Nonsense_Mutation", "frameshift", "Missense_Mutation", "3'UTR", "5'UTR",
    "Splice_Site", "In_Frame_Ins", "In_Frame_Del", "Silent", "IGR", "Intron",
    "RNA", "RNA"
  )
  
  # Apply consequence mappings
  for (i in seq_along(consequence_patterns)) {
    matching_rows <- grep(consequence_patterns[i], data$Consequence)
    data$Consequence[matching_rows] <- consequence_values[i]
  }
  
  # Handle frameshift variants specially (depends on indel direction)
  frameshift_indices <- which(data$Consequence == "frameshift")
  data$Consequence[frameshift_indices] <- ifelse(
    ref_lengths[frameshift_indices] > alt_lengths[frameshift_indices],
    "Frame_Shift_Del",
    "Frame_Shift_Ins"
  )
  
  return(data)
}

calculate_variant_properties <- function(data) {
  ref_lengths <- nchar(data$REF)
  alt_lengths <- nchar(data$ALT)
  max_allele_length <- rowMaxs(cbind(ref_lengths, alt_lengths))
  
  data$END <- data$POS + max_allele_length - 1
  data$Variant_Type <- ifelse(
    max_allele_length == 1, "SNP",
    ifelse(ref_lengths > alt_lengths, "DEL", "INS")
  )
  
  return(list(data = data, ref_lengths = ref_lengths, alt_lengths = alt_lengths))
}



format_maf_output <- function(data, include_hgvsp = FALSE) {
  # Reorder columns
  column_order <- c("Hugo_Symbol", "CHROM", "POS", "END", "STRAND",
                    "Consequence", "Variant_Type", "REF", "REF", "ALT",
                    "Tumor_Sample_Barcode")
  
  if (include_hgvsp) {
    column_order <- c(column_order, "HGVSp_Short")
  }
  
  data <- data[, column_order]
  
  # Rename columns to MAF standard
  standard_names <- c("Hugo_Symbol", "Chromosome", "Start_Position",
                      "End_Position", "Strand", "Variant_Classification",
                      "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1",
                      "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")
  
  if (include_hgvsp) {
    standard_names <- c(standard_names, "HGVSP_SHORT")
  }
  
  colnames(data) <- standard_names
  
  return(data)
}

# Function to create sample metadata
create_sample_metadata <- function(sample_ids) {
  sample_data <- data.frame(
    ID = substr(sample_ids$SAMPLE_ID, 2, 3),
    Sample = sample_ids$SAMPLE_ID,
    Cohort = "Oesophagus",
    SAMPLE_ID = sample_ids$SAMPLE_ID,
    ONCOTREE_CODE = "EGC"
  )
  
  sample_data$ID <- gsub("\\.", "", sample_data$ID)
  
  return(sample_data)
}
