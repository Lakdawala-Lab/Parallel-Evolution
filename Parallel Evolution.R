# Title: Parallel Evolution in Influenza A Virus (IAV) Sequences
# Author: Jennifer E. Jones

# This analysis requires the DECIPER package. This package and its dependencies
# can be downloaded from Bioconductor at http://bioconductor.org/packages/release/bioc/html/DECIPHER.html.

# Load libraries.

library(DECIPHER)
library(ape)
library(phangorn)
library(RColorBrewer)
library(lattice)
library(plotrix)
library(VGAM)

sessionInfo() # Need DECIPHER >= v2.18.1

# Set working directory.

setwd("<<path to IAV directory>>/") # Needs trailing slash.

# FASTA FILE IMPORT AND CLEANUP

# Retrieve all fasta sequences for analysis from the Influenza Research Database. 
# Note that we are beginning with H3N2, but H3N2 and H1N1 
# sequences are interchangeable at this stage.

files <- c("./Human influenza H3N2 segment 1-PB2 sequences.fasta", 
           "./Human influenza H3N2 segment 2-PB1 sequences.fasta",
           "./Human influenza H3N2 segment 3-PA sequences.fasta",
           "./Human influenza H3N2 segment 4-HA sequences.fasta",
           "./Human influenza H3N2 segment 5-NP sequences.fasta",
           "./Human influenza H3N2 segment 6-NA sequences.fasta",
           "./Human influenza H3N2 segment 7-M sequences.fasta",
           "./Human influenza H3N2 segment 8-NS sequences.fasta")

# Fasta file QC.

length <- 8
mRNA <- vector(mode = "list", length = length)
vRNA <- vector(mode = "list", length = length)

for (i in seq_along(vRNA)) {
      mRNA[[i]] <- readDNAStringSet(files[i]) # IAV fasta files are downloaded as positive-sense (mRNA) sequences. 
      vRNA[[i]] <- reverseComplement(mRNA[[i]]) # Convert sequence files to the genomic (vRNA) sequences.
      vRNA[[i]] <- vRNA[[i]][grep("|Subtype:H3N2|", names(vRNA[[i]]), fixed=TRUE)] # Ensure that all fasta files contain the desired subtype (H1N1 or H3N2).
      vRNA[[i]] <- vRNA[[i]][grep("Host:Human", names(vRNA[[i]]), fixed=TRUE)] # Ensure that all fasta files contain human origin IAV.
      
      # Extract the strain names from the fasta files. This step is required to ensure accurate concatenation of full-length genomes.
      names(vRNA[[i]]) <- gsub(".+\\|Strain Name:(.+?)\\|Segment.+", "\\1", names(vRNA[[i]])) 
      
      # Remove duplicate entries from fasta files.
      d <- which(!duplicated(names(vRNA[[i]])))
      vRNA[[i]] <- vRNA[[i]][d]
}
  
# Use matching to remove partial sequences from analysis.

o <- order(sapply(vRNA, length)) # Return the order of the fasta files from least number of sequences to greatest. 

# Return the positions of the first argument that are found in the second argument. NAs indicate sequences
# that are not shared and must be removed.
strainnames <- match(names(vRNA[[o[2]]]), names(vRNA[[o[1]]])) 
names(vRNA[[o[1]]][which(is.na(strainnames))]) # List the sequences in vRNA[[o[2]]] that are not found in vRNA[[o[1]]].
commonnames <- which(!is.na(strainnames))
vRNA[[o[1]]] <- vRNA[[o[1]]][strainnames[commonnames]]
vRNA[[o[2]]] <- vRNA[[o[2]]][commonnames]

# Repeat for remaining segments until matching is complete. Once all eight sets 
# of fasta files are ordered, concatenate full-length genomes and write these files to the working directory.

IAVFL <- xscat(vRNA[[1]], vRNA[[2]], vRNA[[3]], vRNA[[4]], vRNA[[5]], vRNA[[6]], vRNA[[7]], vRNA[[8]])
names(IAVFL) <- names(vRNA[[1]])

# Write all files to the designated working directory. [[Here is a potential stopping point.]]

writeXStringSet(vRNA[[1]], file =  paste(substring(files[1], 1, nchar(files[1]) - nchar("sequences.fasta")), "vRNA sequences.fasta.gz", sep=""))
writeXStringSet(vRNA[[2]], file =  paste(substring(files[2], 1, nchar(files[2]) - nchar("sequences.fasta")), "vRNA sequences.fasta.gz", sep=""))
writeXStringSet(vRNA[[3]], file =  paste(substring(files[3], 1, nchar(files[3]) - nchar("sequences.fasta")), "vRNA sequences.fasta.gz", sep=""))
writeXStringSet(vRNA[[4]], file =  paste(substring(files[4], 1, nchar(files[4]) - nchar("sequences.fasta")), "vRNA sequences.fasta.gz", sep=""))
writeXStringSet(vRNA[[5]], file =  paste(substring(files[5], 1, nchar(files[5]) - nchar("sequences.fasta")), "vRNA sequences.fasta.gz", sep=""))
writeXStringSet(vRNA[[6]], file =  paste(substring(files[6], 1, nchar(files[6]) - nchar("sequences.fasta")), "vRNA sequences.fasta.gz", sep=""))
writeXStringSet(vRNA[[7]], file =  paste(substring(files[7], 1, nchar(files[7]) - nchar("sequences.fasta")), "vRNA sequences.fasta.gz", sep=""))
writeXStringSet(vRNA[[8]], file =  paste(substring(files[8], 1, nchar(files[8]) - nchar("sequences.fasta")), "vRNA sequences.fasta.gz", sep=""))
writeXStringSet(IAVFL, file=paste(substring(files[1], 1, nchar(files[1]) - nchar("segment 1-PB2 sequences.fasta")), "concatenated full-length vRNA sequences.fasta.gz", sep=""))

# Reading files back in after stopping. 

files <- list("./Human influenza H3N2 segment 1-PB2 vRNA sequences.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA sequences.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA sequences.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA sequences.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA sequences.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA sequences.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA sequences.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA sequences.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA sequences.fasta.gz")

length <- 9
vRNA <- vector(mode = "list", length = length)
vRNA <- lapply(files, readDNAStringSet)

# Select strains from desired years to analyze. [[FOR H1N1: 2000-2008 and 2010-2018 were used.]]

year <- gsub(".+/([0-9])", "\\1", names(vRNA[[1]]))
w <- which(year > 1994 & year < 2005) # 1995 - 2004
vRNA1 <- list(vRNA[[1]][w],
              vRNA[[2]][w],
              vRNA[[3]][w],
              vRNA[[4]][w],
              vRNA[[5]][w],
              vRNA[[6]][w],
              vRNA[[7]][w],
              vRNA[[8]][w],
              vRNA[[9]][w])
w <- which(year > 2004 & year < 2015) # 2005 - 2014
vRNA2 <- list(vRNA[[1]][w],
              vRNA[[2]][w],
              vRNA[[3]][w],
              vRNA[[4]][w],
              vRNA[[5]][w],
              vRNA[[6]][w],
              vRNA[[7]][w],
              vRNA[[8]][w],
              vRNA[[9]][w])

# Align sequences and write the alignments to the designated working directory. [[Here is another potential stopping point.]]

vRNA1 <- lapply(vRNA1, AlignSeqs, processors = NULL)
vRNA2 <- lapply(vRNA2, AlignSeqs, processors = NULL)

writeXStringSet(vRNA1[[1]], file =  paste(substring(files[1], 1, nchar(files[1]) - nchar("sequences.fasta.gz")), "MSA - 1995-2004.fasta.gz", sep=""))
writeXStringSet(vRNA1[[2]], file =  paste(substring(files[2], 1, nchar(files[2]) - nchar("sequences.fasta.gz")), "MSA - 1995-2004.fasta.gz", sep=""))
writeXStringSet(vRNA1[[3]], file =  paste(substring(files[3], 1, nchar(files[3]) - nchar("sequences.fasta.gz")), "MSA - 1995-2004.fasta.gz", sep=""))
writeXStringSet(vRNA1[[4]], file =  paste(substring(files[4], 1, nchar(files[4]) - nchar("sequences.fasta.gz")), "MSA - 1995-2004.fasta.gz", sep=""))
writeXStringSet(vRNA1[[5]], file =  paste(substring(files[5], 1, nchar(files[5]) - nchar("sequences.fasta.gz")), "MSA - 1995-2004.fasta.gz", sep=""))
writeXStringSet(vRNA1[[6]], file =  paste(substring(files[6], 1, nchar(files[6]) - nchar("sequences.fasta.gz")), "MSA - 1995-2004.fasta.gz", sep=""))
writeXStringSet(vRNA1[[7]], file =  paste(substring(files[7], 1, nchar(files[7]) - nchar("sequences.fasta.gz")), "MSA - 1995-2004.fasta.gz", sep=""))
writeXStringSet(vRNA1[[8]], file =  paste(substring(files[8], 1, nchar(files[8]) - nchar("sequences.fasta.gz")), "MSA - 1995-2004.fasta.gz", sep=""))
writeXStringSet(vRNA1[[9]], file =  paste(substring(files[9], 1, nchar(files[9]) - nchar("sequences.fasta.gz")), "MSA - 1995-2004.fasta.gz", sep=""))

writeXStringSet(vRNA2[[1]], file =  paste(substring(files[1], 1, nchar(files[1]) - nchar("sequences.fasta.gz")), "MSA - 2005-2014.fasta.gz", sep=""))
writeXStringSet(vRNA2[[2]], file =  paste(substring(files[2], 1, nchar(files[2]) - nchar("sequences.fasta.gz")), "MSA - 2005-2014.fasta.gz", sep=""))
writeXStringSet(vRNA2[[3]], file =  paste(substring(files[3], 1, nchar(files[3]) - nchar("sequences.fasta.gz")), "MSA - 2005-2014.fasta.gz", sep=""))
writeXStringSet(vRNA2[[4]], file =  paste(substring(files[4], 1, nchar(files[4]) - nchar("sequences.fasta.gz")), "MSA - 2005-2014.fasta.gz", sep=""))
writeXStringSet(vRNA2[[5]], file =  paste(substring(files[5], 1, nchar(files[5]) - nchar("sequences.fasta.gz")), "MSA - 2005-2014.fasta.gz", sep=""))
writeXStringSet(vRNA2[[6]], file =  paste(substring(files[6], 1, nchar(files[6]) - nchar("sequences.fasta.gz")), "MSA - 2005-2014.fasta.gz", sep=""))
writeXStringSet(vRNA2[[7]], file =  paste(substring(files[7], 1, nchar(files[7]) - nchar("sequences.fasta.gz")), "MSA - 2005-2014.fasta.gz", sep=""))
writeXStringSet(vRNA2[[8]], file =  paste(substring(files[8], 1, nchar(files[8]) - nchar("sequences.fasta.gz")), "MSA - 2005-2014.fasta.gz", sep=""))
writeXStringSet(vRNA2[[9]], file =  paste(substring(files[9], 1, nchar(files[9]) - nchar("sequences.fasta.gz")), "MSA - 2005-2014.fasta.gz", sep=""))

# CLUSTERING INTO OPERATIONAL TAXONOMIC UNITS (OTUs) AND SEQUENCE SELECTION

# Calculate distances between full-length concatenated sequences 
# and compare clusters with different cutoffs ranging from 95-99% sequence identity.

d1 <- DistanceMatrix(vRNA1[[9]], type="dist", correction="JC", processors = NULL)
d2 <- DistanceMatrix(vRNA2[[9]], type="dist", correction="JC", processors = NULL)

otu1 <- IdClusters(d1, method = "NJ", cutoff = c(0.01, 0.02, 0.03, 0.04, 0.05), type = "clusters", myXStringSet = vRNA1[[9]], processors = NULL)
otu2 <- IdClusters(d2, method = "NJ", cutoff = c(0.01, 0.02, 0.03, 0.04, 0.05), type = "clusters", myXStringSet = vRNA2[[9]], processors = NULL)

# View the number of clusters in each species tree with each cutoff.

c(range(otu1[,1])[2], range(otu1[,2])[2], range(otu1[,3])[2], range(otu1[,4])[2], range(otu1[,5])[2])
c(range(otu2[,1])[2], range(otu2[,2])[2], range(otu2[,3])[2], range(otu2[,4])[2], range(otu2[,5])[2])

# Rerun clustering with desired cutoff (in this case, 97% sequence identity was selected) with type = "both".

otu1 <- IdClusters(d1, method = "NJ", cutoff = 0.03, type = "both", myXStringSet = vRNA1[[9]], processors = NULL)
otu2 <- IdClusters(d2, method = "NJ", cutoff = 0.03, type = "both", myXStringSet = vRNA2[[9]], processors = NULL)

# Write clustering data to working directory.

write.table(otu1[[1]], file = "H3N2 1995-2004 Clusters.csv", sep=",")
write.table(otu2[[1]], file = "H3N2 2005-2014 Clusters.csv", sep=",")

plot(otu2[[2]], leaflab = "none") # Plot species tree with strain names hidden.

# Visually inspect sequence quality and choose representative sequences from each cluster.

w <- which(otu1[[1]] == 1) 
length(w)
m <- match(names(vRNA1[[9]]), rownames(otu1[[1]])[w]) # Locate position in alignment for all sequences in cluster 1.
w <- which(!is.na(m))
BrowseSeqs(vRNA1[[9]][sample(w, 20)], colWidth = 100, highlight = 0) # View an alignment of a subset of sequences from this cluster.

# Subset the alignments of the sequences selected for analysis. otu1 - otu7 represent replicates. 
# ** The sequences analyzed here were entered into an external table, which can be found on GitHub.

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 1995-2004.fasta.gz")

length <- 9
otu <- vector(mode = "list", length = length)

for (i in seq_along(otu)) {
  otu[[i]] <- readDNAStringSet(files[i])
  
  id <- read.csv("H3N2 1995-2004 Strains Analyzed.csv")
  
  m <- match(names(otu[[i]]), id$id_1)
  w <- which(!is.na(m))
  otu1 <- otu[[i]][w]
  writeXStringSet(otu1.1, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 1.fasta.gz", sep=""))
  
  m <- match(names(otu[[i]]), id$id_2)
  w <- which(!is.na(m))
  otu2 <- otu[[i]][w]
  writeXStringSet(otu1.2, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 2.fasta.gz", sep=""))
  
  m <- match(names(otu[[i]]), id$id_3)
  w <- which(!is.na(m))
  otu3 <- otu[[i]][w]
  writeXStringSet(otu1.3, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 3.fasta.gz", sep=""))
  
  m <- match(names(otu[[i]]), id$id_4)
  w <- which(!is.na(m))
  otu4 <- otu[[i]][w]
  writeXStringSet(otu1.4, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 4.fasta.gz", sep=""))
  
  m <- match(names(otu[[i]]), id$id_5)
  w <- which(!is.na(m))
  otu5 <- otu[[i]][w]
  writeXStringSet(otu1.5, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 5.fasta.gz", sep=""))
  
  m <- match(names(otu[[i]]), id$id_6)
  w <- which(!is.na(m))
  otu6 <- otu[[i]][w]
  writeXStringSet(otu1.6, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 6.fasta.gz", sep=""))
  
  m <- match(names(otu[[i]]), id$id_7)
  w <- which(!is.na(m))
  otu7 <- otu[[i]][w]
  writeXStringSet(otu1.7, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 7.fasta.gz", sep=""))
  
}

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014.fasta.gz", 
              "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014.fasta.gz",
              "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014.fasta.gz",
              "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014.fasta.gz",
              "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014.fasta.gz",
              "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014.fasta.gz",
              "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014.fasta.gz",
              "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014.fasta.gz",
              "./Human influenza H3N2 concatenated full-length vRNA MSA - 2005-2014.fasta.gz")

length <- length(files)
otu <- vector(mode = "list", length = length)

for (i in seq_along(otu)) {
  otu[[i]] <- readDNAStringSet(files[i])
  
  id <- read.csv("./H3N2 2005-2014 Strains Analyzed.csv")
  
  m <- match(names(otu[[i]]), id$id_1)
  w <- which(!is.na(m))
  otu1 <- otu[[i]][w]
  writeXStringSet(otu1.1, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 1.fasta.gz", sep=""))
  
  m <- match(names(otu[[i]]), id$id_2)
  w <- which(!is.na(m))
  otu2 <- otu[[i]][w]
  writeXStringSet(otu1.2, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 2.fasta.gz", sep=""))
    
  m <- match(names(otu[[i]]), id$id_3)
  w <- which(!is.na(m))
  otu3 <- otu[[i]][w]
  writeXStringSet(otu1.3, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 3.fasta.gz", sep=""))
  
  m <- match(names(otu[[i]]), id$id_4)
  w <- which(!is.na(m))
  otu4 <- otu[[i]][w]
  writeXStringSet(otu1.4, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 4.fasta.gz", sep=""))
    
  m <- match(names(otu[[i]]), id$id_5)
  w <- which(!is.na(m))
  otu5 <- otu[[i]][w]
  writeXStringSet(otu1.5, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 5.fasta.gz", sep=""))
    
  m <- match(names(otu[[i]]), id$id_6)
  w <- which(!is.na(m))
  otu6 <- otu[[i]][w]
  writeXStringSet(otu1.6, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 6.fasta.gz", sep=""))
    
  m <- match(names(otu[[i]]), id$id_7)
  w <- which(!is.na(m))
  otu7 <- otu[[i]][w]
  writeXStringSet(otu1.7, file= paste(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU 7.fasta.gz", sep=""))
  
}

# RECONSTRUCTING PHYLOGENETIC TREES

# Model testing for maximum-likelihood method. 

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014.fasta.gz")

for (file in seq_along(files)) {
  seg <- read.dna(gzfile(files[file]), "fasta")
  segPD <- phyDat(seg, type = "DNA", levels = NULL)
  mt <- modelTest(segPD)
  write.table(mt, file = paste(substring(files[file], 1, nchar(files[file]) - nchar(".fasta.gz")), "_mt.csv", sep=""), sep=",", quote = FALSE, row.names = F)
}

# Perform bootstrapping.

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 1.fasta.gz", 
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 7.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 1995-2004 OTU 1.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 1995-2004 OTU 2.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 1995-2004 OTU 3.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 1995-2004 OTU 4.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 1995-2004 OTU 5.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 1995-2004 OTU 6.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 1995-2004 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 1.fasta.gz", 
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 2005-2014 OTU 7.fasta.gz")

size <- 1000 # number of bootstrap replicates

for (file in seq_along(files)) {
  vRNA <- readDNAStringSet(files[file])
  d <- DistanceMatrix(vRNA, type="dist", correction="JC")
  tree <- IdClusters(d, method="ML", type="dendrogram", model="HKY85", myXStringSet=vRNA) 
  
  f <- function(x) {
    if (is.null(attributes(x)$leaf)) {   
      x0 <- paste(sort(unlist(x)), collapse=" ")   
      x1 <- f(x[[1]])                                 
      x2 <- f(x[[2]])
      return(list(x0, x1, x2))
    } else {
      return(NULL)
    }
  }
  
  pBar <- txtProgressBar(style=3)
  bootstraps <- list()
  l <- unique(width(vRNA))
  for (i in seq_len(size)) { # size iterations
    r <- sample(l, replace=TRUE) # Sample the length of the alignment with replacement (i.e. resample sample value multiple times).
    at <- IRanges(r, width=1)            
    vRNA2 <- extractAt(vRNA, at) # Extract subsequence from alignment at the ranges of positions specified in 'at' (i.e. scramble the genome).
    vRNA2 <- lapply(vRNA2, unlist)
    vRNA2 <- DNAStringSet(vRNA2)
    
    d <- DistanceMatrix(vRNA2, type="dist", correction="JC", verbose=FALSE)
    temp <- IdClusters(d, method="ML", type="dendrogram", model="HKY85", myXStringSet=vRNA2, verbose=FALSE)   # Create OTUs from each iteration of the sampling/extraction. Now have 100 resampled dendrograms.
    bootstraps[[i]] <- unlist(f(temp))   
    setTxtProgressBar(pBar, i/size)
  }
  
  bootstraps <- table(unlist(bootstraps))     
  original <- unlist(f(tree))       
  hits <- bootstraps[original]     
  names(hits) <- original
  w <- which(is.na(hits))
  if (length(w) > 0)         
    hits[w] <- 0
  hits <- round(hits/size*100)
  
  f <- function(x) {
    if (is.null(attributes(x)$leaf)) {
      attr(x, "edgetext") <- as.character(hits[paste(sort(unlist(x)), collapse=" ")])
    }
    return(x)
  }
  d <- dendrapply(tree, f)                  
  attr(d, "edgetext") <- NULL                
  
  WriteDendrogram(d, file = paste(substring(files[file], 1,  nchar(files[file]) - nchar(".fasta.gz")), "_MLtree", sep=""))
  par(mai= c(1, 0.5, 0.25, 0.1)) # bottom, left, top, right
  plot(d, edgePar=list(t.cex=0.5), nodePar=list(lab.cex=0.7, pch=NA), edge.root=FALSE)
  
}

# CALCULATION OF ROBINSON-FOULDS DISTANCES.

# Trees built using the DECIPHER package must be read in by the ape package. 
# This converts the trees from objects of type = 'dendrogram' to type = 'phylo', 
# which introduces syntax errors that must be corrected prior to analysis.

# Working with one set of replicates from 2005-2014 at a time:

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 1_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 1_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 1_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 1_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 1_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 1_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 1_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 1_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 2005-2014 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]])) # Replace strain names with cluster IDs in tree so that replicates can be compared.
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_1[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_1[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_1[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_1[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_1[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_1[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_1[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_1[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_1[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

# Combine all eight trees from each vRNA segment into one object and perform Robinson-Foulds distance calculation.

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf1 <- dist.topo(unroot(trees), method = "PH85")

# Repeat for remaining replicates.

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 2_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 2_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 2_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 2_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 2_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 2_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 2_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 2_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 2005-2014 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_2[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_2[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_2[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_2[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_2[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_2[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_2[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_2[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_2[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf2 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 3_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 3_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 3_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 3_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 3_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 3_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 3_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 3_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 2005-2014 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_3[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_3[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_3[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_3[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_3[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_3[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_3[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_3[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_3[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf3 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 4_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 4_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 4_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 4_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 4_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 4_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 4_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 4_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 2005-2014 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]])) 
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_4[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_4[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_4[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_4[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_4[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_4[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_4[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_4[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_4[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf4 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 5_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 5_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 5_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 5_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 5_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 5_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 5_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 5_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 2005-2014 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_5[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_5[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_5[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_5[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_5[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_5[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_5[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_5[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_5[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf5 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 6_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 6_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 6_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 6_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 6_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 6_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 6_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 6_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 2005-2014 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_6[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_6[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_6[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_6[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_6[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_6[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_6[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_6[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_6[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf6 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 7_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 7_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 7_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 7_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 7_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 7_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 7_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 7_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 2005-2014 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_7[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_7[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_7[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_7[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_7[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_7[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_7[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_7[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_7[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf7 <- dist.topo(unroot(trees), method = "PH85")

l <- list(as.matrix(rf1), as.matrix(rf2), as.matrix(rf3), as.matrix(rf4), as.matrix(rf5), as.matrix(rf6), as.matrix(rf7))
rf <- do.call(cbind, l)
rf <- array(rf, dim=c(dim(l[[1]]), length(l)))

rownames(l[[1]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(l[[1]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
rownames(l[[2]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(l[[2]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
rownames(l[[3]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(l[[3]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
rownames(l[[4]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(l[[4]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
rownames(l[[5]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(l[[5]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
rownames(l[[6]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(l[[6]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
rownames(l[[7]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(l[[7]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))

write.table(l[[1]], file = "./H3N2 2005-2014 d values_replicate 1.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[2]], file = "./H3N2 2005-2014 d values_replicate 2.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[3]], file = "./H3N2 2005-2014 d values_replicate 3.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[4]], file = "./H3N2 2005-2014 d values_replicate 4.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[5]], file = "./H3N2 2005-2014 d values_replicate 5.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[6]], file = "./H3N2 2005-2014 d values_replicate 6.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[7]], file = "./H3N2 2005-2014 d values_replicate 7.csv", sep=",", quote = FALSE, row.names = T)

rfmean <- apply(rf, c(1, 2), mean) # Compute the mean Robinson-Foulds distances across replicates.
sem <- apply(rf, c(1, 2), std.error) # Compute the standard error.
write.table(rfmean, file = "./H3N2 2005-2014 d values_mean.csv", sep=",", quote = FALSE, row.names = T)
write.table(rfsd, file = "./H3N2 2005-2014 d values_sd.csv", sep=",", quote = FALSE, row.names = T)

# The data can be visualized in a number of ways. Here is how to generate the dendrograms 
# from H3N2 viruses from 2005-2014.

dend <- matrix(rfmean, nrow = 8)
row.names(dend) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
colnames(dend) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
tree <- IdClusters(dend, method="UPGMA", type="dendrogram")

par(mai = c(1, 0.75, 0.75, 0.25))  # bottom, left, top, right
plot(tree, xlab = "vRNA Segment", type = "triangle", main = "2005-2014", ylim= c(0, 5))
WriteDendrogram(tree, file = "./H3N2 2005-2014 d value dendrogram")

rownames(rfmean) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(rfmean) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))

# Here is how to generate the heatmaps presented in the manuscript.
# First, remove the duplicated portions of the data.

rfmean <- rfmean[,-8] # Remove the NS column.
rfmean <- rfmean[-1,] # Remove the PB2 row.
rfmean[1,][2:7] <- NA
rfmean[2,][3:7] <- NA
rfmean[3,][4:7] <- NA
rfmean[4,][5:7] <- NA
rfmean[5,][6:7] <- NA
rfmean[6,][7] <- NA

sem <- sem[,-8]
sem <- sem[-1,]
sem[1,][2:7] <- NA
sem[2,][3:7] <- NA
sem[3,][4:7] <- NA
sem[4,][5:7] <- NA
sem[5,][6:7] <- NA
sem[6,][7] <- NA

rownames(sem) <- as.character(c("PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(sem) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M"))


palette <- colorRampPalette(brewer.pal(7, "YlGnBu"))(20)
cols <- palette[20:1]

levelplot(rfmean, 
          at=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
          xlab = "",
          ylab = "",
          col.regions = cols)
axis(4) #  Draw the y-axis to the right of the plot area

levelplot(sem, 
          at=c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
          xlab = "",
          ylab = "",
          col.regions = cols)
axis(4)

# Now repeat for 1995-2004:

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 1_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 1_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 1_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 1_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 1_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 1_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 1_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 1_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 1995-2004 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]])) # Replace strain names with cluster IDs in tree so that replicates can be compared.
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_1[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_1[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_1[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_1[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_1[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_1[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_1[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_1[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_1[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

# Combine all eight trees from each vRNA segment into one object and perform Robinson-Foulds distance calculation.

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf1 <- dist.topo(unroot(trees), method = "PH85")

# Repeat for remaining replicates.

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 2_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 2_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 2_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 2_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 2_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 2_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 2_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 2_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 1995-2004 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_2[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_2[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_2[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_2[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_2[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_2[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_2[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_2[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_2[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf2 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 3_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 3_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 3_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 3_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 3_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 3_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 3_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 3_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 1995-2004 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_3[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_3[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_3[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_3[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_3[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_3[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_3[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_3[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_3[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf3 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 4_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 4_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 4_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 4_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 4_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 4_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 4_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 4_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 1995-2004 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]])) 
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_4[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_4[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_4[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_4[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_4[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_4[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_4[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_4[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_4[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf4 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 5_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 5_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 5_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 5_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 5_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 5_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 5_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 5_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 1995-2004 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_5[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_5[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_5[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_5[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_5[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_5[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_5[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_5[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_5[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf5 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 6_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 6_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 6_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 6_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 6_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 6_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 6_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 6_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 1995-2004 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_6[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_6[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_6[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_6[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_6[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_6[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_6[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_6[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_6[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf6 <- dist.topo(unroot(trees), method = "PH85")

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 1995-2004 OTU 7_MLtree", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 1995-2004 OTU 7_MLtree",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 1995-2004 OTU 7_MLtree",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 1995-2004 OTU 7_MLtree",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 1995-2004 OTU 7_MLtree",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 1995-2004 OTU 7_MLtree",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 1995-2004 OTU 7_MLtree",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 1995-2004 OTU 7_MLtree")

length <- 8
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  otu <- read.csv("./H3N2 1995-2004 Strains Analyzed - ape.csv")
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_7[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_7[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_7[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_7[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_7[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_7[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_7[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_7[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_7[9])
      treetips[j] = 9
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(1.2, 1, 1.5, 1))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf7 <- dist.topo(unroot(trees), method = "PH85")

l <- list(as.matrix(rf1), as.matrix(rf2), as.matrix(rf3), as.matrix(rf4), as.matrix(rf5), as.matrix(rf6), as.matrix(rf7))
rf <- do.call(cbind, l)
rf <- array(rf, dim=c(dim(l[[1]]), length(l)))

rownames(l[[1]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
colnames(l[[1]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
rownames(l[[2]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
colnames(l[[2]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
rownames(l[[3]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
colnames(l[[3]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
rownames(l[[4]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
colnames(l[[4]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
rownames(l[[5]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
colnames(l[[5]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
rownames(l[[6]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
colnames(l[[6]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
rownames(l[[7]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))
colnames(l[[7]]) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "SIX", "M", "NS"))

write.table(l[[1]], file = "./H3N2 1995-2004 d values_replicate 1.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[2]], file = "./H3N2 1995-2004 d values_replicate 2.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[3]], file = "./H3N2 1995-2004 d values_replicate 3.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[4]], file = "./H3N2 1995-2004 d values_replicate 4.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[5]], file = "./H3N2 1995-2004 d values_replicate 5.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[6]], file = "./H3N2 1995-2004 d values_replicate 6.csv", sep=",", quote = FALSE, row.names = T)
write.table(l[[7]], file = "./H3N2 1995-2004 d values_replicate 7.csv", sep=",", quote = FALSE, row.names = T)

rfmean <- apply(rf, c(1, 2), mean) # Compute the mean Robinson-Foulds distances across replicates.
sem <- apply(rf, c(1, 2), std.error) # Compute the standard error.
write.table(rfmean, file = "./H3N2 1995-2004 d values_mean.csv", sep=",", quote = FALSE, row.names = T)
write.table(rfsd, file = "./H3N2 1995-2004 d values_sd.csv", sep=",", quote = FALSE, row.names = T)

# Generating dendrograms from H3N2 viruses from 1995-2004.

dend <- matrix(rfmean, nrow = 8)
row.names(dend) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
colnames(dend) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
tree <- IdClusters(dend, method="UPGMA", type="dendrogram")

par(mai = c(1, 0.75, 0.75, 0.25))  # bottom, left, top, right
plot(tree, xlab = "vRNA Segment", type = "triangle", main = "1995-2004", ylim= c(0, 5))
WriteDendrogram(tree, file = "./H3N2 1995-2004 d value dendrogram")


rownames(rfmean) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(rfmean) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))


# Generating heatmaps from H3N2 trees from 1995-2004.

rfmean <- rfmean[,-8] # Remove the NS column.
rfmean <- rfmean[-1,] # Remove the PB2 row.
rfmean[1,][2:7] <- NA
rfmean[2,][3:7] <- NA
rfmean[3,][4:7] <- NA
rfmean[4,][5:7] <- NA
rfmean[5,][6:7] <- NA
rfmean[6,][7] <- NA

sem <- sem[,-8]
sem <- sem[-1,]
sem[1,][2:7] <- NA
sem[2,][3:7] <- NA
sem[3,][4:7] <- NA
sem[4,][5:7] <- NA
sem[5,][6:7] <- NA
sem[6,][7] <- NA

rownames(sem) <- as.character(c("PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(sem) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M"))


palette <- colorRampPalette(brewer.pal(7, "YlGnBu"))(20)
cols <- palette[20:1]

levelplot(rfmean, 
          at=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
          xlab = "",
          ylab = "",
          col.regions = cols)
axis(4) #  Draw the y-axis to the right of the plot area

levelplot(sem, 
          at=c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
          xlab = "",
          ylab = "",
          col.regions = cols)
axis(4)

# COMPARISON OF D VALUES 

# Null distribution of 1000 randomly generated unrooted trees:

n <- 12 # Designate the number of tips in each tree. For H1N1, the number of tips = 9.
trees <- rmtree(n, rooted = F, tip.label = 1:n, N = 1000)
rf.rm <- dist.topo(unroot(trees), method = "PH85")
hist(rf.rm, main = "Distribution of d values", xlab = "d value")

# Linear regression analysis for determination of the 95% confidence cutoff.

a <- which(rf.rm == 0)
b <- which(rf.rm == 2)
c <- which(rf.rm == 4)
d <- which(rf.rm == 6)
e <- which(rf.rm == 8)
f <- which(rf.rm == 10)
g <- which(rf.rm == 12)
h <- which(rf.rm == 14)
i <- which(rf.rm == 16)
j <- which(rf.rm == 18)

# The best method for transforming the data was Yeo-Johnson:

y <- c(length(a), length(b), length(c), length(d), length(e), 
       length(f), length(g), length(h), length(i), length(j)) # check for zeros
y2 <- yeo.johnson(y, lambda = 0) # 

x <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18)

# Linear regression indicates that y2 is the best fit for the data.
lm(y2 ~ x) # adj R-squared = 0.94 with a p-value of 2.865e-06. 

# The equation of the line is y = 0.80934x - 2.22734
((yeo.johnson((0.05*length(rf.rm)), lambda = 0) + 2.11)/0.8019) # the d value at which the 95% CI is exceeded.

plot(x, y2, pch = 16, cex = 1.5, col = "black", xlab = "d value", 
     ylab = "Adjusted Frequency", xlim = c(0, 20), ylim = c(0, 14),
     main = "leaf number = 12 (H3N2)", cex.lab = 1.2, cex.main = 1.5)
abline(lm(y2 ~ x))

# Next, calculate the Robinson-Foulds distance between the two dendrograms.

tree1 <- read.tree("./H3N2 1995-2004 d value dendrogram")
tree2 <- read.tree("./H3N2 2005-2014 d value dendrogram")
trees <- c(tree1, tree2)
dist.topo(unroot(trees), method = "PH85")

# ANALYSIS OF PROTEIN-CODING SEQUENCES

# CDS fasta files were downloaded from FluDB that correspond to all H3N2 2005-2014 strains analyzed:

files <- c("./huH3N2 2005-2014 PB1 CDS OTU 1.fasta",
           "./huH3N2 2005-2014 PB1 CDS OTU 2.fasta",
           "./huH3N2 2005-2014 PB1 CDS OTU 3.fasta",
           "./huH3N2 2005-2014 PB1 CDS OTU 4.fasta",
           "./huH3N2 2005-2014 PB1 CDS OTU 5.fasta",
           "./huH3N2 2005-2014 PB1 CDS OTU 6.fasta",
           "./huH3N2 2005-2014 PB1 CDS OTU 7.fasta",
           "./huH3N2 2005-2014 PA CDS OTU 1.fasta",
           "./huH3N2 2005-2014 PA CDS OTU 2.fasta",
           "./huH3N2 2005-2014 PA CDS OTU 3.fasta",
           "./huH3N2 2005-2014 PA CDS OTU 4.fasta",
           "./huH3N2 2005-2014 PA CDS OTU 5.fasta",
           "./huH3N2 2005-2014 PA CDS OTU 6.fasta",
           "./huH3N2 2005-2014 PA CDS OTU 7.fasta",
           "./huH3N2 2005-2014 NP CDS OTU 1.fasta",
           "./huH3N2 2005-2014 NP CDS OTU 2.fasta",
           "./huH3N2 2005-2014 NP CDS OTU 3.fasta",
           "./huH3N2 2005-2014 NP CDS OTU 4.fasta",
           "./huH3N2 2005-2014 NP CDS OTU 5.fasta",
           "./huH3N2 2005-2014 NP CDS OTU 6.fasta",
           "./huH3N2 2005-2014 NP CDS OTU 7.fasta",
           "./huH3N2 2005-2014 NA CDS OTU 1.fasta",
           "./huH3N2 2005-2014 NA CDS OTU 2.fasta",
           "./huH3N2 2005-2014 NA CDS OTU 3.fasta",
           "./huH3N2 2005-2014 NA CDS OTU 4.fasta",
           "./huH3N2 2005-2014 NA CDS OTU 5.fasta",
           "./huH3N2 2005-2014 NA CDS OTU 6.fasta",
           "./huH3N2 2005-2014 NA CDS OTU 7.fasta",
           "./huH3N2 2005-2014 PB2 CDS OTU 1.fasta",
           "./huH3N2 2005-2014 PB2 CDS OTU 2.fasta",
           "./huH3N2 2005-2014 PB2 CDS OTU 3.fasta",
           "./huH3N2 2005-2014 PB2 CDS OTU 4.fasta",
           "./huH3N2 2005-2014 PB2 CDS OTU 5.fasta",
           "./huH3N2 2005-2014 PB2 CDS OTU 6.fasta",
           "./huH3N2 2005-2014 PB2 CDS OTU 7.fasta",
           "./huH3N2 2005-2014 HA CDS OTU 1.fasta",
           "./huH3N2 2005-2014 HA CDS OTU 2.fasta",
           "./huH3N2 2005-2014 HA CDS OTU 3.fasta",
           "./huH3N2 2005-2014 HA CDS OTU 4.fasta",
           "./huH3N2 2005-2014 HA CDS OTU 5.fasta",
           "./huH3N2 2005-2014 HA CDS OTU 6.fasta",
           "./huH3N2 2005-2014 HA CDS OTU 7.fasta",
           "./huH3N2 2005-2014 M1 CDS OTU 1.fasta",
           "./huH3N2 2005-2014 M1 CDS OTU 2.fasta",
           "./huH3N2 2005-2014 M1 CDS OTU 3.fasta",
           "./huH3N2 2005-2014 M1 CDS OTU 4.fasta",
           "./huH3N2 2005-2014 M1 CDS OTU 5.fasta",
           "./huH3N2 2005-2014 M1 CDS OTU 6.fasta",
           "./huH3N2 2005-2014 M1 CDS OTU 7.fasta",
           "./huH3N2 2005-2014 M2 CDS OTU 1.fasta",
           "./huH3N2 2005-2014 M2 CDS OTU 2.fasta",
           "./huH3N2 2005-2014 M2 CDS OTU 3.fasta",
           "./huH3N2 2005-2014 M2 CDS OTU 4.fasta",
           "./huH3N2 2005-2014 M2 CDS OTU 5.fasta",
           "./huH3N2 2005-2014 M2 CDS OTU 6.fasta",
           "./huH3N2 2005-2014 M2 CDS OTU 7.fasta",
           "./huH3N2 2005-2014 NS1 CDS OTU 1.fasta",
           "./huH3N2 2005-2014 NS1 CDS OTU 2.fasta",
           "./huH3N2 2005-2014 NS1 CDS OTU 3.fasta",
           "./huH3N2 2005-2014 NS1 CDS OTU 4.fasta",
           "./huH3N2 2005-2014 NS1 CDS OTU 5.fasta",
           "./huH3N2 2005-2014 NS1 CDS OTU 6.fasta",
           "./huH3N2 2005-2014 NS1 CDS OTU 7.fasta",
           "./huH3N2 2005-2014 NS2 CDS OTU 1.fasta",
           "./huH3N2 2005-2014 NS2 CDS OTU 2.fasta",
           "./huH3N2 2005-2014 NS2 CDS OTU 3.fasta",
           "./huH3N2 2005-2014 NS2 CDS OTU 4.fasta",
           "./huH3N2 2005-2014 NS2 CDS OTU 5.fasta",
           "./huH3N2 2005-2014 NS2 CDS OTU 6.fasta",
           "./huH3N2 2005-2014 NS2 CDS OTU 7.fasta")

for (file in seq_along(files)) {
  cds <- readDNAStringSet(files[file])
  names(cds) <- gsub(".+\\|Strain Name:(.+?)\\|Protein.+", "\\1", names(cds))
  msa <- AlignTranslation(cds, type = "AAStringSet")
  writeXStringSet(msa, file =  paste(substring(files[file], 1, nchar(files[file]) - nchar("CDS sequences.fasta")), "AA MSA.fasta", sep=""))
}

size <- 1000 # number of bootstrap replicates

for (file in seq_along(files)) {
  dna <- readDNAStringSet(files[file])
  
  # remove gaps from alignment:
  
  dna2 <- DNAStringSet(character(length(dna)))
  for (i in seq_along(dna)) {
    w <- which(strsplit(as.character(dna[[i]]), "", fixed=T)[[1]] != "-")
    dna2[[i]] <- dna[[i]][w]
    names(dna2) <- names(dna)
  }
  
  dna <- AlignTranslation(dna2, type = "AAStringSet")
  d <- DistanceMatrix(dna, type="dist", correction="JC")
  tree <- IdClusters(d, method="NJ", type="dendrogram", myXStringSet=dna)   # the type of the object 'tree' is a list
  
  f <- function(x) {
    if (is.null(attributes(x)$leaf)) {                 # If the contents of x$leaf are NULL, do this:
      x0 <- paste(sort(unlist(x)), collapse=" ")       # Create a new list containing all members of x sorted and pasted together.
      x1 <- f(x[[1]])                                 
      x2 <- f(x[[2]])
      return(list(x0, x1, x2))
    } else {
      return(NULL)
    }
  }
  
  pBar <- txtProgressBar(style=3)             # A progress bar for the later functions.
  bootstraps <- list()
  l <- unique(width(dna))
  for (i in seq_len(size)) { # size iterations
    r <- sample(l, replace=TRUE)              # Sample the length of the alignment with replacement (meaning you can resample sample value multiple times).
    at <- IRanges(r, width=1)                 # Create an IRanges object with r ranges.
    dna2 <- extractAt(dna, at)                # Extract subsequence from alignment at the ranges of positions specified in at (essentially randomize the genome).
    dna2 <- lapply(dna2, unlist)
    dna2 <- AAStringSet(dna2)
    
    d <- DistanceMatrix(dna2, type="dist", correction="JC", verbose=FALSE)
    temp <- IdClusters(d, method="NJ", type="dendrogram", myXStringSet=dna2, verbose=FALSE)   # Create OTUs from each iteration of the sampling/extraction. Now have 100 resampled dendrograms.
    bootstraps[[i]] <- unlist(f(temp))        # Create a new list of OTUs for each iteration that is sorted and pasted together.
    setTxtProgressBar(pBar, i/size)
  }
  
  bootstraps <- table(unlist(bootstraps))     
  original <- unlist(f(tree))                 # Use original tree to set indices of each member of the tree.
  hits <- bootstraps[original]                # Determine how many of the bootstrap replicates match the original tree.
  names(hits) <- original
  w <- which(is.na(hits))
  if (length(w) > 0)                          # Set NAs to zero, then convert to percentages.
    hits[w] <- 0
  hits <- round(hits/size*100)
  
  f <- function(x) {
    if (is.null(attributes(x)$leaf)) {
      attr(x, "edgetext") <- as.character(hits[paste(sort(unlist(x)), collapse=" ")])
    }
    return(x)
  }
  d <- dendrapply(tree, f)                   # Dendrapply applies a function to all nodes of a dendrogram.
  attr(d, "edgetext") <- NULL                # Now have d, an object of type list and class dendrogram. 
  # Should be able to pass to tanglegram function. 
  # May need to use as.phylo to pass to write.tree, 
  # but could instead use WriteDendrogram function to bypass this issue.
  
  WriteDendrogram(d, file = paste(substring(files[file], 
                                            1, 
                                            nchar(files[file]) - nchar(".fasta.gz")), 
                                  "_AA NJtree", 
                                  sep=""))
  par(mai= c(1, 0.5, 0.25, 0.1)) # bottom, left, top, right
  plot(d, edgePar=list(t.cex=0.5), nodePar=list(lab.cex=0.7, pch=NA), edge.root=FALSE)
  
}

# The Robinson-Foulds distance between all trees can now be calculated using the script beginning on line 509.









