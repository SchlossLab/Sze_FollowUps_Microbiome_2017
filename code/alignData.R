#Output shared file that has the same rows as the two combined metadata files. 


# Load data to be aligned

shared <- read.table("data/process/unmatched.an.unique_list.shared", stringsAsFactors=F, header=T)

followup <- read.delim("data/raw/metadata/followUps_metadata.txt", header=T, sep='\t')

initial <- read.delim("data/raw/metadata/initials_metadata.tsv", header=T, sep='\t')

# Get the sample names that are used in initial and follow up samples
tempList <- c(initial[, 'sample'], followup[, 'followUp'])

# Re-align the shared to match these two files
shared <- shared[match(shared$Group, tempList), ]

# Remove extrac sample that is in the metadata and not in shared for check
tempList <- tempList[match(shared$Group, tempList), ]

# Check to make sure that it worked
stopifnot(shared$Group == tempList)

# Write new shared file
write.table(shared, file="data/process/final.shared", quote=F, sep='\t', row.names=F)

