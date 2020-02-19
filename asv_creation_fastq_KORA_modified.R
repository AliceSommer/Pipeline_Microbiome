
# set working directory
# setwd("/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
# 
# biocLite(suppressUpdates = FALSE)
# biocLite("ShortRead", suppressUpdates = FALSE)
BiocManager::install("ShortRead", suppressUpdates = FALSE)
library(ShortRead)

# biocLite("devtools")
BiocManager::install("devtools")
library("devtools")

devtools::install_github("benjjneb/dada2")
library('dada2'); packageVersion("dada2")

library(Biostrings)

#######################################################################################################

##################
# Getting ready  #
##################

# not all KORA subjects gave permission to use their data
# use the same sample ids than in the OTU table (created at the TUM)
samples_to_keep <- readRDS('dada2code/samples_to_keep.rds')
length(samples_to_keep) # 2034

# path to the sequence (fastq_files) folder 
# path <- "dada2code/KORA_fasq_files"

list.files(path)[1:10]
tail(list.files(path))

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
# remove "DNAstab", "Wasser" with regular expression
fnFs <- sort(list.files(path, pattern="^[0-9.]+_R1.fastq.gz$", full.names = TRUE))
fnFs[1:5]
fnRs <- sort(list.files(path, pattern="^[0-9.]+_R2.fastq.gz$", full.names = TRUE))
length(fnFs) # 2186

# Forward and reverse fastq filenames for the samples pipetted with the robot
path_robot <- paste(path, "/sample-sequences", sep = '')
fnFs_Robot <- sort(list.files(path_robot, pattern="[0-9]{5}.roboter@F.fastq", full.names = TRUE))
fnRs_Robot <- sort(list.files(path_robot, pattern="[0-9]{5}.roboter@R.fastq", full.names = TRUE))

# Remove the file names of the fastq list that are in the roboter lists
identical(sapply(strsplit(fnFs_Robot, "[.]"), "[[", 2), sapply(strsplit(fnRs_Robot, "[.]"), "[[", 2))
robot_sample_ids <- sapply(strsplit(fnFs_Robot, "[.]"), "[[", 2)

fnFs_small <- fnFs[!grepl(paste0(robot_sample_ids, collapse = "|"), fnFs)]
fnRs_small <- fnRs[!grepl(paste0(robot_sample_ids, collapse = "|"), fnRs)]

# Create two lists of paths for the R and F with manunal and robot sequencing combined
# with only the "approved" sample ids from the KORA OTU-table
fnFs_combined <- c(fnFs_small,fnFs_Robot)
fnFs_final <- fnFs_combined[grepl(paste0(samples_to_keep, collapse = "|"), fnFs_combined)]

fnRs_combined <- c(fnRs_small,fnRs_Robot)
fnRs_final <- fnRs_combined[grepl(paste0(samples_to_keep, collapse = "|"), fnRs_combined)]
  
length(fnFs_final)

# Extract sample names from origial filenames where format: SAMPLENAME_XXX.fastq 
# Sample name function
get.sample.name <- function(fname) {
  if(grepl('roboter', fname)==TRUE) { 
    sub_str = strsplit(basename(fname), ".roboter")[[1]][1]
    substr(sub_str, 5, 9)
  } else { 
    strsplit(basename(fname), "_")[[1]][1]
  }
}

sample.names <- unname(sapply(fnFs_final, get.sample.name))
head(sample.names); length(sample.names)

#####################
# Identify primers #
####################

FWD <- "CCTACGGGNGGCWGCAG"  
REV <- "GACTACHVGGGTATCTAATCC" 

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Count the number of times the primers appear in the forward and reverse read, 
# while considering all possible primer orientations.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs_final[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs_final[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs_final[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs_final[[1]]))

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs_final[length(fnFs_final)]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs_final[length(fnRs_final)]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs_final[length(fnFs_final)]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs_final[length(fnRs_final)]))

#################################
# Inspect read quality profiles #
#################################

l <- length(fnFs_final)
plotQualityProfile(fnFs_final[(l-1):l])
plotQualityProfile(fnFs_final[1:2])
plotQualityProfile(fnRs_final[(l-1):l])
plotQualityProfile(fnRs_final[c(200,201)])

###################
# Filter and Trim #
###################

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# trim the nucleotides at the end of the sequence (often bad quality)
out <- filterAndTrim(fnFs_final, filtFs,
                     fnRs_final, filtRs,
                     trimLeft = c(17,21), # remove the primers
                     truncLen=c(245,245), 
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
# compress: TRUE, the output fastq file(s) are gzipped.
# multithread: TRUE, input files are filtered in parallel via mclapply
head(out[,2],20)
sum(out[,2] < 1)
min(out[,2])

#####################
# Learn error rates #
#####################

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#####################################################
# Dereplication, sample inference, merge pair reads #
#####################################################

### Sample inference
ddF <- dada(filtFs, err=errF, multithread=TRUE)
ddR <- dada(filtRs, err=errR, multithread=TRUE)

### merge pair reads
mergers <- mergePairs(ddF, filtFs, ddR, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

############################
# Construct sequence table #
############################

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

###################
# Remove chimeras #
###################

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

table(nchar(getSequences(seqtab.nochim)))

seq_depth <- apply(seqtab.nochim, 1, function(x) sum(x))
hist(seq_depth)
min(seq_depth)
max(seq_depth)
median(seq_depth)
table(seq_depth)

getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(ddF, getN), sapply(ddR, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track,10)

length(unique(row.names(seqtab.nochim)))

# rename one sample (86110) for later merge with KORA covariates
rownames(seqtab.nochim) <- sample.names
rownames(seqtab.nochim)[1666]
rownames(seqtab.nochim)[1666] <- "86110"
rownames(seqtab.nochim)[1666]

###################
# Assign Taxonomy #
###################

# download version 128 of Silva's reference database 
# (https://zenodo.org/record/824551#.Xk1b_y2ZNUM)

# assignTaxPath <- "/Users/ajs654/Desktop/KORA_fasq_files/silva_download/silva_nr_v128_train_set.fa"
# assignSpeciesPath <- "/Users/ajs654/Desktop/KORA_fasq_files/silva_download/silva_species_assignment_v128.fa"

tax <- assignTaxonomy(seqtab.nochim, assignTaxPath, multithread=TRUE)
taxa <- addSpecies(tax, assignSpeciesPath)
tail(taxa)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#########
# SAVE #
########

# Write to disk
# saveRDS(seqtab.nochim, "dada2output/dada2output2020/seqtab2020.rds")
# saveRDS(taxa, "dada2output2020/taxa2020.rds") 
