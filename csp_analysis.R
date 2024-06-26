
library(ggplot2)
library(dplyr)
library(tidyr)
library(Biostrings) 
library(progress)
library(fs)
library(DECIPHER) 

######################################################################
#-----------------------INPUT AND PREP DATA--------------------------#
######################################################################

#######################################################
# 1.- metadata
#######################################################

db <- readxl::read_excel('HFS22_ID_run_list_csp_030524.xlsx')
colnames(db)[1]<- c("NIDA2")  #just to avoid changing variables in the script from Simone's analysis

length(unique(db$NIDA2))

# 1a.- complete run_id column with data from the actual runs

# Initialize an empty data frame to store the results
result_df <- data.frame(NIDA2 = character(), NEW_run_id = character(), stringsAsFactors = FALSE)

directory_path <- "../results_v0.1.8_RESMARKERS_FIX/"

# Create a progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(dir_ls(path = directory_path, regexp = "_RESULTS_v0.1.8$"))
)

# Iterate through the folders ending with _RESULTS_v0.1.8
for (folder_path in dir_ls(path = directory_path, regexp = "_RESULTS_v0.1.8$")) {
  pb$tick()  # Update progress bar
  
  folder_name <- path_file(folder_path)
  file_path <- file.path(folder_path, "amplicon_coverage.txt")
  
  cat("/n")
  print(folder_name)
  
  # Read the contents of /quality_report/amplicon_stats.txt
  if (file.exists(file_path)) {
    sample_coverage_content <- readLines(file_path)
    
    # Truncate each line after the first tab (/t) and return unique values and edit nida format to match the one from the db
    truncated_values <- unique(sapply(strsplit(sample_coverage_content, "/t"), function(x) x[1]))
    truncated_values <- gsub("_S.*$", "", truncated_values)
    truncated_values <- gsub("_", ".", truncated_values)
    truncated_values <- gsub("-", ".", truncated_values)
    truncated_values <- gsub("N", "", truncated_values)
    
    # Extract NIDA2 from runs using grep
    nida_values <- grep(paste(db$NIDA2, collapse = "|"), truncated_values, value = TRUE) #not all samples from runs are needed, only those in the db, hence this
    
    # Create a data frame with NIDA2 and folder_name
    if (length(nida_values) > 0) {
      temp_df <- data.frame(NIDA2 = nida_values, NEW_run_id = folder_name, stringsAsFactors = FALSE)
      
      # Append the results to the main data frame
      result_df <- rbind(result_df, temp_df)
    }
  }
}

result_df$NEW_run_id <- sub("_RESULTS_v0.1.8$", "", result_df$NEW_run_id)

#this should be TRUE if the same number of nidas is in the db and the results from the grep
length(result_df$NIDA2) == length(db$NIDA2) #it's not....
length(result_df$NIDA2) - length(db$NIDA2)

#check for repeated nidas 
repeated_nidas <- names(table(result_df$NIDA2)[table(result_df$NIDA2) > 1]) #there are duplicate nidas OOF
repeated_nidas_df<- result_df[result_df$NIDA2 %in% repeated_nidas,]
repeated_nidas_df <- repeated_nidas_df[order(repeated_nidas_df$NIDA2), ]
length(repeated_nidas_df$NIDA2)
length(unique(repeated_nidas_df$NIDA2))

#ask team about these nidas
#write.csv(repeated_nidas_df, "repeated_nidas.csv", row.names = F)


#######################################################
# 2.- genomic data from all runs
#######################################################

runs <- unique(paste0(result_df$NEW_run_id, "_RESULTS_v0.1.8_FILTERED"))

folder_path <- paste0("./", runs)
allele_data_files <- file.path(folder_path, "allele_data_global_max_0_filtered.csv")

# Import filtered allele data tables to the list
allele_data_list <- list()
for (file in allele_data_files) {
  allele_data <- read.csv(file)
  allele_data_list <- append(allele_data_list, list(allele_data))
}

#format the imported dfs
for (i in seq_along(allele_data_list)) {
  df <- allele_data_list[[i]]
  
  df$sampleID <- gsub("_S.*$", "", df$sampleID)
  df$sampleID <- gsub("_", ".", df$sampleID)
  df$sampleID <- gsub("-", ".", df$sampleID)
  df$sampleID <- gsub("N", "", df$sampleID)
  df$sampleID <- gsub("//.0", "", df$sampleID)
  
  colnames(df)[1] <- "sample_id"
  
  df <- df %>% ### TEHERE IS A BUG WITH THE MASK THAT GENERATES MORE ALLELES THAN THERE REALLY ARE. THIS SNIPPET OF CODE COLLAPSES REPETITIONSAND SUMS THE READS AND FREQS. CRITICAL!!!
    group_by(sample_id, locus, pseudo_cigar, asv) %>%
    summarize(reads = sum(reads),
              norm.reads.locus = sum(norm.reads.locus)) %>%
    mutate(allele = paste(locus, ".", row_number(), sep = ""))
  
  
  # Update the modified data frame back in the list
  df<- as.data.frame(df)
  allele_data_list[[i]] <- df
}

names(allele_data_list) <- runs


#get rid of replicate nidas keep the sample with the most reads across runs for each one of the replicates
sum_reads <- function(df) {
  aggregate(reads ~ sample_id, data = df, sum)
}


summed_reads <- lapply(allele_data_list, sum_reads)

# Change colnames of reads for the name of the respective df
summed_reads <- lapply(names(summed_reads), function(df_name) {
  df <- summed_reads[[df_name]]
  names(df)[2] <- df_name
  return(df)
})

# Merge all data frames by sample_id
merged_df_dups <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), summed_reads)

#keep rows with replicates
merged_df_dups <- merged_df_dups[rowSums(is.na(merged_df_dups)) < length(colnames(merged_df_dups)) - 2, ] #keep rows with more than 14 NAs (that is, that have reads in more than one run)

# Exclude the first column and find the column with the maximum value for each row
merged_df_dups$BEST_RUN <- colnames(merged_df_dups)[apply(merged_df_dups[-1], 1, which.max) + 1]

#remove replicates, keep the best
for (i in 1:nrow(merged_df_dups)) {
  best_run <- merged_df_dups$BEST_RUN[i]
  sample_id <- merged_df_dups$sample_id[i]
  
  # Loop through allele_data_list
  for (j in seq_along(allele_data_list)) {
    df <- allele_data_list[[j]]
    
    # Exclude the df named after BEST_RUN
    if (names(allele_data_list)[j] != best_run) {
      allele_data_list[[j]] <- df[!(df$sample_id %in% sample_id), ]
    }
  }
}

#final check:
summed_reads <- lapply(allele_data_list, sum_reads)
summed_reads <- lapply(names(summed_reads), function(df_name) {
  df <- summed_reads[[df_name]]
  names(df)[2] <- df_name
  return(df)
})
merged_df_dups <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), summed_reads)
merged_df_dups <- merged_df_dups[rowSums(is.na(merged_df_dups)) < length(colnames(merged_df_dups)) - 2, ] #keep rows with more than 14 NAs (that is, that have reads in more than one run)

if (dim(merged_df_dups)[1] == 0){
  print("NO MORE REPLICATES.")
}else{
  "grab a coffee"
}

# Define a function to filter rows based on criteria: remove controls, basically
#filter_rows <- function(df) {
  # Filter rows based on criteria
 # filtered_df <- df[!grepl("^[A-Za-z]|3D", df$sample_id), ]
 # return(filtered_df)
#}

# Apply the filter_rows function to each dataframe in allele_data_list
#allele_data_list <- lapply(allele_data_list, filter_rows)

#save allele_data_list
saveRDS(allele_data_list, "allele_data_list.RDS")

#######################################################
# 3.- GENOMIC + DB MERGING, FILTERING ETC. (DATA PREP)
#######################################################

allele_data_list <- readRDS("allele_data_list.RDS")

# concat all dataframes together
combined_df <- bind_rows(allele_data_list)

# calculate n.alleles for each locus of each sample if not done already during contaminant filtering
if (!("n.alleles" %in% colnames(combined_df))){
  combined_df <- combined_df %>%
    group_by(sample_id, locus) %>%
    mutate(n.alleles = n_distinct(pseudo_cigar))
}

# merge with metadata
colnames(combined_df)[1]<- c("NIDA2")
combined_df_merged <- merge(combined_df, db[c("NIDA2", "year", "province", "region")], by="NIDA2", all.y =T) #forcing adding all db nidas

#check for columns with NAs. THESE SAMPLES WERE BAD QUALITY AND THUS FILTERED OUT DURING CONTAMINANTS FILTERING
removed_samples <- combined_df_merged[is.na(combined_df_merged$Category),]
removed_samples$NIDA2

# Remove rows with NIDA2 matching values in removed_samples
combined_df_merged <- combined_df_merged[!combined_df_merged$NIDA2 %in% removed_samples$NIDA2, ]


#sanity check
if( sum(!(combined_df_merged$NIDA2 %in% db$NIDA2)) == 0){
  print("All nidas in combined_merged_df are also the metadata db. No weird samples ✔")
}else{
  print("grab another coffee.")
}


## FURTHER FILTERING

# 1) MAF filering (< 0.01)
combined_df_merged <- combined_df_merged[combined_df_merged$norm.reads.locus  > 0.01, ]

## Keep only csp amplicons
csp_amps <- c("Pf3D7_03_v3-221185-221419-1B", "Pf3D7_03_v3-221463-221712-1B", "Pf3D7_03_v3-221295-221532-2")
combined_df_merged_csp <- combined_df_merged[combined_df_merged$locus %in% csp_amps,]

# 2) check for coverage (>100 loci with >threshold reads)
# Define thresholds for read depth
thresholds <- c(25, 50, 100, 200)
count_list <- list()

# Loop over each threshold
for (threshold in thresholds) {
  # Calculate unique loci counts for the current threshold
  count <- combined_df_merged_csp %>%
    group_by(NIDA2, locus) %>%
    summarize(total_reads = sum(reads)) %>%
    group_by(NIDA2) %>%
    filter(total_reads > threshold) %>%
    summarize(!!paste("unique_loci_", threshold, sep = "") := n_distinct(locus))
  
  count_list[[paste("count_", threshold, sep = "")]] <- count
}

result_df <- Reduce(function(x, y) left_join(x, y, by = "NIDA2"), count_list)

#count cells above 1 for each column with the "unique" substring: here, i'm calculating the sample size for each reads count threshold using 100 loci as cutoff: samples with <100 read count per loci below the threshold should be removed
count_above_1 <- function(x) sum(x >= 1, na.rm = TRUE)
unique_columns <- grep("unique", colnames(result_df), value = TRUE)

count_results <- result_df %>%
  summarise_at(vars(unique_columns), count_above_1)

count_results

# decided going with a threshold of >= 100 read depth,
samples_to_keep <- result_df[result_df$unique_loci_100 >= 1, ]$NIDA2

combined_df_merged_csp <- combined_df_merged_csp[combined_df_merged_csp$NIDA2 %in% samples_to_keep, ]


# 3) remove bad loci: those that are not in at least 100 samples with a read depth of 100. 
# Group by NIDA2 and locus, then summarize the total reads
locus_read_depth <- combined_df_merged_csp %>%
  group_by(NIDA2, locus) %>%
  summarize(read_depth = sum(reads))

# Count the number of samples with read depth greater than 100 for each locus
locus_counts <- locus_read_depth %>%
  group_by(locus) %>%
  summarize(samples_above_100 = sum(read_depth > 100))

loci_to_keep <- locus_counts$locus


combined_df_merged_csp <- combined_df_merged_csp[combined_df_merged_csp$locus %in% loci_to_keep, ]

#recount n.alleles after filtering
combined_df_merged_csp <- combined_df_merged_csp %>%
  group_by(NIDA2, locus) %>%
  mutate(n.alleles = n_distinct(pseudo_cigar))

combined_df_merged_csp <- as.data.frame(combined_df_merged_csp)



#rename alleles
combined_df_merged_csp$allele <- paste0(combined_df_merged_csp$locus, "_", combined_df_merged_csp$pseudo_cigar)

# small edit for cabo delgado
combined_df_merged_csp$province <- gsub(" ", "_", combined_df_merged_csp$province) 


# FILTERING RESULTS
SS <- length(unique(combined_df_merged_csp$NIDA2))
cat("Final sample size is", as.character(SS))

LC <- length(unique(combined_df_merged_csp$locus))
cat("csp loci count is", as.character(LC))


#save genomic db
saveRDS(combined_df_merged_csp, "combined_df_merged_csp_only_testing.RDS")



#######################################################
# 4.- CHECK SAMPLE SIZES
#######################################################

combined_df_merged_csp <- readRDS("combined_df_merged_csp_only_testing.RDS")

sample_size_provinces <- combined_df_merged_csp %>%
  group_by(province) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_provinces

sample_size_regions <- combined_df_merged_csp %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions



######################################################################
#-----------------------     ALIGNMENT     --------------------------#
######################################################################

combined_df_merged_csp <- readRDS("combined_df_merged_csp_only_testing.RDS")

# csp reference gene (https://plasmodb.org/plasmo/app/record/gene/PF3D7_0304600#category:gene-structure):
csp_ref <- readDNAStringSet("csp_plasmodb.fasta")


#to avoid doing excessive amounts of alignment, perform only on unique alleles (goes from 3000+ to just 86!), fill table later
unique_alleles <- combined_df_merged_csp[c("locus", "asv")]
unique_alleles <- distinct(unique_alleles)
dim(unique_alleles)


# Loop through each sequence in unique_alleles$asv
amp_rev_comps<-c()

for (seq in unique_alleles$asv) {
  ok <- DNAString(seq)
  rev_comp_seq <- as.character(reverseComplement(ok)) # Reverse complement the sequence
  amp_rev_comps <- c(amp_rev_comps, rev_comp_seq)
}


unique_alleles$amp_reverse_compleemntary <- amp_rev_comps

#get dna alignments
aligned_amp_rev_comp<- c()
aligned_amp_rev_comp_aln<- c()

for (amprevcomp in unique_alleles$amp_reverse_compleemntary){
  
  # Create example DNA sequences as a DNAStringSet
  sequences <- DNAStringSet(c(csp_ref, amprevcomp))
  
  # Perform pairwise sequence alignment
  alignment <- AlignSeqs(sequences)
  
  aligned_amp_rev_comp <- c(aligned_amp_rev_comp, as.character(alignment[2]))
}

unique_alleles$aligned_amp_rev_comp <- aligned_amp_rev_comp


#output full alignment just because
full_alignment_concat <- DNAStringSet(c(csp_ref, aligned_amp_rev_comp))


# rename unique alleles
counter <- 1

modified_names <- character(length(unique_alleles$locus))

for (unique_element in unique(unique_alleles$locus)) {
  indices <- which(unique_alleles$locus == unique_element)
  modified_names[indices] <- paste0("HAP_", counter:(counter + length(indices) - 1), "_", unique_element) #puedes cambiar HAP_ por el prefijo que quieras
  counter <- counter + length(indices)
}

# Replace loci names
modified_names <- gsub("Pf3D7_03_v3-221185-221419-1B", "amp1", modified_names)
modified_names <- gsub("Pf3D7_03_v3-221295-221532-2", "amp2", modified_names)
modified_names <- gsub("Pf3D7_03_v3-221463-221712-1B", "amp3", modified_names)

unique_alleles$hap_name <- modified_names

names(full_alignment_concat)[2:length(full_alignment_concat)]<- unique_alleles$hap_name


writeXStringSet(full_alignment_concat, "full_csp_aligment_Dd2_test1K.fasta", format = "fasta")


# trim amplicons that are outside the csp gene 
csp_length <- length(csp_ref$`PF3D7_0304600.1  | Plasmodium falciparum 3D7 | circumsporozoite (CS) protein | CDS | length=1194`)

trim_string <- function(x, csp_length) {
  if (nchar(x) > csp_length) {
    return(substr(x, 1, csp_length))
  } else {
    return(x)
  }
}

# Apply the function to each element
unique_alleles$aligned_amp_rev_comp <- sapply(unique_alleles$aligned_amp_rev_comp, trim_string, csp_length)

#check if length of sequences is the same now:
if (length(unique(lapply(unique_alleles$aligned_amp_rev_comp, nchar))) == 1){
  print("Amplicons have been successfully trimmed to the reference length")
}else{
  print("grab a coffee.")
}



######################################################################
#------------------------   TRANSLATION    --------------------------#
######################################################################

# translate aligned amplicons (reverse compliment)
getCodons <- function(myAln) {
  seqs <- as.character(myAln)
  len <- width(myAln)[1]
  starts <- seq(from=1, to=len, by=3)
  ends <- starts + 2
  myViews <- lapply(myAln, function(x) { 
    Views(x, starts, ends)
  })
  myCodons <- lapply(myViews, function(x) {
    as.character(DNAStringSet(x))
  })
  myCodons
}

## translateCodons - takes a character vector of codons as input, outputs the corresponding amino acids
translateCodons <- function(myCodons, unknownCodonTranslatesTo="-") {
  ## make new genetic code
  gapCodon <- "-"
  names(gapCodon) <- "---"
  my_GENETIC_CODE <- c(GENETIC_CODE, gapCodon)
  
  ## translate the codons
  pep <- my_GENETIC_CODE[myCodons]
  
  ## check for codons that were not possible to translate, e.g. frameshift codons
  if (sum(is.na(pep))>0) {
    cat("/nwarning - there were codons I could not translate. Using this character", unknownCodonTranslatesTo, "/n/n")
    pep[ which(is.na(pep)) ] <- unknownCodonTranslatesTo
  }
  
  ## prep for output
  pep <- paste(pep, collapse="")
  return(pep)
}

## wrap the getCodons and translateCodons functions together into one:
translateGappedAln <- function(myAln, unknownCodonTranslatesTo="-") {
  myCodons <- getCodons(myAln)
  myAAaln <- AAStringSet(unlist(lapply(myCodons, translateCodons, unknownCodonTranslatesTo=unknownCodonTranslatesTo)))
  return(myAAaln)
}

#translate
aas <- translateGappedAln(unique_alleles$aligned_amp_rev_comp, unknownCodonTranslatesTo="-") #change to x to see imcomplete codons if needed
unique_alleles$translated_aligned_amp_rev_comp <- as.character(aas)

# identify all non-synonymous mutations 
csp_ref_prot <- translate(csp_ref)
aa_alignment <- c(csp_ref_prot, aas)

#rename seqs
names(aa_alignment)[2:length(aa_alignment)]<- unique_alleles$hap_name

#export aminoacid alignment WITH REFERENCE! just because
writeXStringSet(aa_alignment, "full_csp_aligment_Dd2_test1K.faa", format = "fasta")

aa_matrix <- matrix("", nrow = length(aa_alignment), ncol = width(aa_alignment[1]))
dim(aa_matrix) #check: good

# Fill the matrix character by character
for (i in 1:length(aa_alignment)) {
  current_seq <- as.character(aa_alignment[i])
  
  for (j in 1:nchar(current_seq)) {
    aa_matrix[i, j] <- substr(current_seq, j, j)
  }
}

aa_df <- as.data.frame(aa_matrix)
colnames(aa_df) <- 1:ncol(aa_df)

# Compare each sequence with the reference
reference_seq_df <- aa_df[1,]
aa_df <- aa_df[-1,]
rownames(aa_df)<- c()

nsym <- data.frame(rowname = character(), non_synonymous_codon = character(), stringsAsFactors = FALSE)

for (i in 1:nrow(aa_df)) {
  current_seq <- as.character(aa_df[i, ])
  differences <- which(current_seq != reference_seq_df & (current_seq != "-" & current_seq != "X"))
  
  if (length(differences) > 0) {
    nsym <- rbind(nsym, data.frame(rowname = rownames(aa_df)[i], non_synonymous_codon = differences))
  } else {
    nsym <- rbind(nsym, data.frame(rowname = rownames(aa_df)[i], non_synonymous_codon = NA))
  }
}

nsym$ALT <- character(nrow(nsym))
nsym$REF<- character(nrow(nsym))

# Populate 'ALT' column based on the coordinates in nsym
for (i in 1:nrow(nsym)) {
  if (is.na(nsym$non_synonymous_codon[i])) {
    nsym$ALT[i] <- NA
  } else {
    row_num <- as.numeric(nsym$rowname[i])
    col_num <- as.numeric(nsym$non_synonymous_codon[i])
    nsym$ALT[i] <- aa_df[row_num, col_num]
  }
}

# Populate 'REF' column based on the coordinates in nsym
for (i in 1:nrow(nsym)) {
  if (!is.na(nsym$non_synonymous_codon[i])) {
    row_num <- 1
    col_num <- as.numeric(nsym$non_synonymous_codon[i])
    nsym$REF[i] <- reference_seq_df[row_num, col_num]
  } else {
    
  }
}

unique_alleles$rowname<- 1:length(rownames(unique_alleles))
unique_alleles_complete <- merge(nsym, unique_alleles, by = "rowname", all.x = TRUE)

#unique_alleles_complete <- unique_alleles_complete[complete.cases(unique_alleles_complete$ALT), ] #removed NA rows (don't have a nsym mutation so not needed)

#keep reference sequence for calcualting allele freqs later.
# Replace NAs and blank values with "WT"
unique_alleles_complete <- unique_alleles_complete %>%
  mutate_all(~replace(., is.na(.) | . == "", "WT"))

# Identify rows where both `asv` and `hap_name` are "WT"
wt_rows <- unique_alleles_complete %>%
  filter(REF == "WT" & ALT == "WT")

# Keep only one "WT" row per unique `locus`
wt_rows_filtered <- wt_rows %>%
  group_by(locus) %>%
  filter(row_number()==1)
  ungroup()

# Non-"WT" rows
non_wt_rows <- unique_alleles_complete %>%
  filter(!(REF == "WT" & ALT == "WT"))

# Combine the filtered "WT" rows with the non-"WT" rows
unique_alleles_complete <- bind_rows(wt_rows_filtered, non_wt_rows)


# Merge based on "locus" and "asv" 
merged_data <- merge(combined_df_merged_csp, unique_alleles_complete, by = c("locus", "asv"), all = T) #all = T for including samples with only synonymous mutations. remove if needed. NA means no nsym mutations found.


######################################################################
#---------------------   OUTPUT FORMATTING    -----------------------#
######################################################################

#sort and subset final table
FINAL_TABLE_sorted <- merged_data[order(merged_data$NIDA2, merged_data$locus, merged_data$pseudo_cigar), ]
FINAL_TABLE_sorted <- FINAL_TABLE_sorted[c("NIDA2", "locus", "allele", "asv", "pseudo_cigar", "reads", "norm.reads.locus", "n.alleles", "non_synonymous_codon", "REF", "ALT", "region", "province", "hap_name")]

FINAL_TABLE_sorted <- FINAL_TABLE_sorted[!is.na(FINAL_TABLE_sorted$non_synonymous_codon),]

#output results
write.csv(FINAL_TABLE_sorted, "csp_nsym_mutations_Dd2_Test1K.csv", row.names = F)

