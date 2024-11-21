
##################################################
# amplicon_seq_analysis.R
# Written by Richard Li and Dorthy Fang 10/19/2022, updated 1/14/2023 (v3), updated 8/13/2024 (v5)
# For compiling mutation sites across Azenta
# amplicon sequencing data tables
##################################################

# SETTINGS
#install.packages("readxl")
require(readxl)
manual_ref_seq <- NA # Paste in ref seq here, or default to ref seq in source
# dataframe
filepath <- "INSERT FILE NAME HERE.xlsx"

###################################################
# FUNCTION FOR GETTING INDICES OF SEQUENCE MISMATCH
###################################################
# Function takes in single reference and target sequence and returns
# a named list including all indices of snps, insertions, and deletions, and
# their lengths
find_mismatch_indices <- function(ref_seq, targ_seq){
  
  # Create empty numeric vectors we can append values to
  snp_indices <- numeric()
  insertion_indices <- numeric()     # start index of insertion
  insertion_lengths <- numeric()
  deletion_indices <- numeric()      # start index of deletion
  deletion_indices_full <- numeric() # all indices of deletion
  deletion_lengths <- numeric()
  
  ###############################################
  # FIND INDICES AND LENGTHS OF INSERTIONS
  # Index returned is site where insertion begins
  ###############################################
  # Boolean vector classifying letters in targ_seq as lowercase
  # (insertions, TRUE) or uppercase (FALSE)
  is_lowercase <- unlist(strsplit(targ_seq,'')) %in% letters
  
  # Loop that identifies location of first insertion + length, then removes it
  # from is_lowercase, looping until no remaining insertions.
  while(TRUE %in% is_lowercase){
    
    # What is the start index of the first remaining insertion?
    current_start_index <- which(is_lowercase == TRUE)[1]
    
    # Set the end index of this insertion to the start index, and increment
    # index until we find end of insertion (no longer lowercase)
    current_end_index <- current_start_index
    while(is_lowercase[current_end_index + 1] == TRUE){
      current_end_index <- current_end_index + 1
      # If we reach the end of the whole string break and do not continue
      if(current_end_index == length(is_lowercase)) break
    }
      
    current_length <- current_end_index - current_start_index + 1
    
    # Append current index + length to vectors
    insertion_indices <- append(insertion_indices, current_start_index)
    insertion_lengths <- append(insertion_lengths, current_length)
    
    # Remove the current insertion from is_lowercase and start loop again
    is_lowercase <- is_lowercase[-c(current_start_index:current_end_index)]
  }
  
  ##############################################
  # FIND INDICES AND LENGTHS OF DELETIONS
  # Index returned is site where deletion begins
  ##############################################
  # Remove insertions from targ_seq before finding deletion indices
  targ_seq_del <- gsub("[a-z ]", "", targ_seq)
  
  # Classifies letters in targ_seq to dashes (deletions, TRUE) or uppercase
  # letters (FALSE)
  is_dash <- unlist(strsplit(targ_seq_del,'')) == "-"
  
  # If deletion found
  if(TRUE %in% is_dash){
    
    # Iterate over letters in target sequence
    i <- 1
    while(i <= length(is_dash)){
      
      # If not dash, increment i and continue
      if(!is_dash[i]){
        i <- i + 1
      }
      
      # If dash, increment i until you hit the next uppercase letter
      # (end of deletions). Track how many increments are done.
      else{
        start_index_targ <- i
        while(is_dash[i] & i <= length(is_dash)){
          i <- i + 1
        }
        end_index_targ <- i
        
        # target indices equivalent to reference indices for deletions
        start_index_ref <- start_index_targ
        end_index_ref <- end_index_targ
        
        # Add converted indices and insertion lengths to vectors
        deletion_indices <- append(deletion_indices, start_index_ref)
        deletion_indices_full <- append(deletion_indices_full,
                                        start_index_ref:end_index_ref-1)
        deletion_lengths <- append(deletion_lengths,
                                   end_index_ref - start_index_ref)
      }
    }
  }
  
  ####################################
  # FIND INDICES OF MUTATIONS
  # Index returned is site of mutation
  ####################################
  # Use version of targ_seq with insertions removed
  # Split target and reference sequences into vectors of single characters
  targ_seq_del_split <- unlist(strsplit(targ_seq_del, ""))
  ref_seq_split <- unlist(strsplit(ref_seq, ""))
  
  # Find locations in targ_seq_del where there are dashes and copy over to
  # same locations in ref_seq
  dash_indices <- which(unlist(strsplit(targ_seq_del, "")) == "-")
  ref_seq_split[dash_indices] <- "-"
  
  # Compare target and reference sequences
  snp_indices <- which(ref_seq_split != targ_seq_del_split)
  
  # Compile to named list to return
  all_indices <- list(
    snp_indices,
    insertion_indices,
    insertion_lengths,
    deletion_indices,
    deletion_indices_full,
    deletion_lengths
  )
  names(all_indices) <- c("snp_indices","insertion_indices","insertion_lengths",
                          "deletion_indices", "deletion_indices_full",
                          "deletion_lengths")
  
  return(all_indices)
}

###############################################################################
###############################################################################

# MAIN CODE

file <- read_xlsx(filepath)
# Set reference sequence to manual input or take from source data
ref_seq <- ifelse(is.na(manual_ref_seq), file$WideReference[1], manual_ref_seq)


###################################################
# FIND INSERTIONS, DELETIONS, AND SNPS FOR ALL RUNS
###################################################
# Create table of outputs from find_mismatch_indices()
mismatch_runs <- data.frame(targ_seq = file$TargetSequence,
                            ref_seq = ref_seq,
                            reads = file$Reads,
                            snp_indices = character(nrow(file)),
                            insertion_indices = character(nrow(file)),
                            insertion_lengths = character(nrow(file)),
                            deletion_indices = character(nrow(file)),
                            deletion_indices_full = character(nrow(file)),
                            deletion_lengths = character(nrow(file)))

# Find indices of mismatch locations and store into dataframe. When vectors
# are multiple elements longs separate elements with | symbol
for(i in 1:nrow(mismatch_runs)){
  if(i == 1) print("Finding insertions, deletions, and mutations")
  if(i %% 1000 == 0) print(paste0("Finished row: ", i))
  row_mismatches <- find_mismatch_indices(mismatch_runs$ref_seq[i],
                                          mismatch_runs$targ_seq[i])
  mismatch_runs$snp_indices[i] <-
    ifelse(length(row_mismatches$snp_indices) == 0, NA,
           paste(row_mismatches$snp_indices, collapse = "|"))
  mismatch_runs$insertion_indices[i] <-
    ifelse(length(row_mismatches$insertion_indices) == 0, NA,
           paste(row_mismatches$insertion_indices, collapse = "|"))
  mismatch_runs$insertion_lengths[i] <-
    ifelse(length(row_mismatches$insertion_lengths) == 0, NA,
           paste(row_mismatches$insertion_lengths, collapse = "|"))
  mismatch_runs$deletion_indices[i] <-
    ifelse(length(row_mismatches$deletion_indices) == 0, NA,
           paste(row_mismatches$deletion_indices, collapse = "|"))
  mismatch_runs$deletion_indices_full[i] <-
    ifelse(length(row_mismatches$deletion_indices_full) == 0, NA,
           paste(row_mismatches$deletion_indices_full, collapse = "|"))  
  mismatch_runs$deletion_lengths[i] <-
    ifelse(length(row_mismatches$deletion_lengths) == 0, NA,
           paste(row_mismatches$deletion_lengths, collapse = "|"))
}

############################################
# AGGREGATE SEQUENCE MISMATCHES ACROSS SITES
############################################
# Create table compiling number of reads across sites
site_analysis <- data.frame(index = 1:nchar(ref_seq),
                            nucleotide = unlist(strsplit(ref_seq, "")),
                            snp_reads = numeric(nchar(ref_seq)),
                            snp_percent = numeric(nchar(ref_seq)),
                            insertion_reads = numeric(nchar(ref_seq)),
                            insertion_percent = numeric(nchar(ref_seq)),
                            deletion_reads = numeric(nchar(ref_seq)),
                            deletion_percent = numeric(nchar(ref_seq)))
for(i in 1:nrow(mismatch_runs)){
  if(i == 1) print("Aggregating mismatches across sites")
  if(i %% 1000 == 0) print(paste0("Finished row: ", i))
  # Get SNP indices from that run, split by character "|". Note that we need
  # to double escape the symbol "|" to have it recognized
  snp_sites <- as.numeric(unlist(strsplit(
    mismatch_runs$snp_indices[i], split = "\\|")))
  # At each site in the snp_sites, increment the reads according to the run
  site_analysis$snp_reads[snp_sites] <-
    site_analysis$snp_reads[snp_sites] + mismatch_runs$reads[i]
  
  insertion_sites <- as.numeric(unlist(strsplit(
    mismatch_runs$insertion_indices[i], split = "\\|")))
  site_analysis$insertion_reads[insertion_sites] <-
    site_analysis$insertion_reads[insertion_sites] + mismatch_runs$reads[i]
  
  deletion_sites <- as.numeric(unlist(strsplit(
    mismatch_runs$deletion_indices_full[i], split = "\\|")))
  site_analysis$deletion_reads[deletion_sites] <-
    site_analysis$deletion_reads[deletion_sites] + mismatch_runs$reads[i]
}

###############
# OVERALL STATS
###############
# Get percentages, format to floating point number
total_reads <- sum(mismatch_runs$reads)
site_analysis$snp_percent <-
  formatC(site_analysis$snp_reads/total_reads*100, digits = 3, format = 'f')
site_analysis$insertion_percent <-
  formatC(site_analysis$insertion_reads/total_reads*100, digits = 3,
          format = 'f')
site_analysis$deletion_percent <-
  formatC(site_analysis$deletion_reads/total_reads*100, digits = 3,
          format = 'f')

write.csv(mismatch_runs,"INSERT FILE NAME HERE.csv", row.names=FALSE)
write.csv(site_analysis,"INSERT FILE NAME HERE.csv", row.names=FALSE)
