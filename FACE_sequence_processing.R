### Microbiomes associated with wheat grown at the FACE site
# 2017 samples used v5_6 and ITS1b primers, on a combined sequencing run (2x250 kit; MiSeq Run130)
# 2018 samples used v4 and the 2x250 sequencing kit (bacteria; MiSeq Run159) / ITS2 primers and the 2x300 kit (fungi; MiSeq Run160)
# Note that the SRA has merged what were originally separate file sets (2018 samples) for bacteria and for fungi
# Since there is going to be an overwhelming effect of target amplicon confounded with year, split and process separately by year

rm(list = ls())

# Load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

.cran_packages <- c("Biostrings", "car", "ggplot2", "gridExtra", "hexbin", "knitr", "lme4", 
                    "lmerTest", "pheatmap", "picante", "plyr", "RColorBrewer", 
                    "reshape2", "ShortRead", "tidyr", "vegan")
.bioc_packages <- c("apeglm", "BiocStyle", "dada2", "DECIPHER", "decontam", "DESeq2", 
                    "phangorn", "phyloseq","vsn")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# Package versions:
sessionInfo()
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)

# Matrix products: default

# locale:
# [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252    LC_MONETARY=English_Canada.1252
# [4] LC_NUMERIC=C                    LC_TIME=English_Canada.1252    

# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] vsn_3.58.0                  phangorn_2.6.3              DESeq2_1.30.1               decontam_1.10.0            
# [5] DECIPHER_2.18.1             RSQLite_2.2.7               dada2_1.18.0                Rcpp_1.0.6                 
# [9] BiocStyle_2.18.1            apeglm_1.12.0               tidyr_1.1.3                 ShortRead_1.48.0           
# [13] GenomicAlignments_1.26.0    SummarizedExperiment_1.20.0 Biobase_2.50.0              MatrixGenerics_1.2.1       
# [17] matrixStats_0.58.0          Rsamtools_2.6.0             GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
# [21] BiocParallel_1.24.1         reshape2_1.4.4              RColorBrewer_1.1-2          plyr_1.8.6                 
# [25] picante_1.8.2               nlme_3.1-152                vegan_2.5-7                 lattice_0.20-41            
# [29] permute_0.9-5               ape_5.5                     pheatmap_1.0.12             lmerTest_3.1-3             
# [33] lme4_1.1-26                 Matrix_1.3-2                knitr_1.33                  hexbin_1.28.2              
# [37] gridExtra_2.3               ggplot2_3.3.3               Biostrings_2.58.0           XVector_0.30.0             
# [41] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1         phyloseq_1.34.0   

path.0.project <- file.path("C:/Seq_Data/FACE")
setwd(path.0.project)
set.seed(1024)
alpha = 0.05

# download accession list from the SRA Run Selector web interface ("SRR_Acc_List.txt")
# in terminal:
# ~/sratoolkit.2.9.0-ubuntu64/bin/prefetch --option-file SRR_Acc_List.txt
# sra_list=({9124593..9124870})
# for sra_id in ${sra_list[@]}
# do
# ~/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump --split-files ~/ncbi/public/sra/SRR${sra_id}.sra
# done
# gzip /data/Research/FACE/1.data/*.fastq

path.1.data <- file.path(path.0.project, "1.data")
list.files(path.1.data)

# Sort ensures forward/reverse reads are in same order
names.1.data.R1 <- sort(list.files(path.1.data, pattern = "_1.fastq.gz"))
names.1.data.R2 <- sort(list.files(path.1.data, pattern = "_2.fastq.gz"))

get.sample.name <- function(names) strsplit(basename(names), "_")[[1]][1]
sample.names <- unname(sapply(names.1.data.R1, get.sample.name))

metadata <- read.csv("Metadata.csv", header = TRUE)
metadata <- merge(metadata, data.frame(sample.names), by.x = "Label.SRA", by.y = "sample.names")
rownames(metadata) <- metadata$Label.SRA
# 278 observations X 10 variables

accessions <- data.frame(names.1.data.R1)
names(accessions) <- "Label.SRA"
accessions$Label.SRA <- substr(accessions$Label.SRA, 1, 10)
metadata.Y17 <- merge(metadata[metadata$year == "Y17", ], accessions)
rownames(metadata.Y17) <- metadata.Y17$Label.SRA
names.1.data.Y17.R1 <- paste0(metadata.Y17$Label.SRA, "_1.fastq.gz")
names.1.data.Y17.R2 <- paste0(metadata.Y17$Label.SRA, "_2.fastq.gz")

metadata.Y18 <- merge(metadata[metadata$year == "Y18", ], accessions)
rownames(metadata.Y18) <- metadata.Y18$Label.SRA
names.1.data.Y18.R1 <- paste0(metadata.Y18$Label.SRA, "_1.fastq.gz")
names.1.data.Y18.R2 <- paste0(metadata.Y18$Label.SRA, "_2.fastq.gz")

# expand names to include full path
names.1.data.Y17.R1 <- file.path(path.1.data, names.1.data.Y17.R1)
names.1.data.Y17.R2 <- file.path(path.1.data, names.1.data.Y17.R2)
names.1.data.Y18.R1 <- file.path(path.1.data, names.1.data.Y18.R1)
names.1.data.Y18.R2 <- file.path(path.1.data, names.1.data.Y18.R2)

# Identify primer sequences
v4.f <- "GTGCCAGCMGCCGCGGTAA"
v4.f.rc <- dada2:::rc(v4.f)
v4.r <- "GGACTACHVGGGTWTCTAAT"
v4.r.rc <- dada2:::rc(v4.r)

v5_6.f <- "AACMGGATTAGATACCCKG"
v5_6.f.rc <- dada2:::rc(v5_6.f)
v5_6.r <- "GGGTTGCGCTCGTTGC"
v5_6.r.rc <- dada2:::rc(v5_6.r)

# R1 seems to lead off with my reverse primers
ITS1.f <- "GCTGCGTTCTTCATCGATGC"
ITS1.f.rc <- dada2:::rc(ITS1.f)
ITS1.r <- "CTTGGTCATTTAGAGGAAGTAA"
ITS1.r.rc <- dada2:::rc(ITS1.r)

ITS2.f <- "GATGAAGAACGYAGYRAA"
ITS2.f.rc <- dada2:::rc(ITS2.f)
ITS2.r <- "CTBTTVCCKCTTCACTCG"
ITS2.r.rc <- dada2:::rc(ITS2.r)

path.2.cut.16S <- file.path(path.0.project, "2.cut.16S")
names.2.cut.16S.Y17.R1 <- file.path(path.2.cut.16S, basename(names.1.data.Y17.R1))
names.2.cut.16S.Y17.R2 <- file.path(path.2.cut.16S, basename(names.1.data.Y17.R2))
names.2.cut.16S.Y18.R1 <- file.path(path.2.cut.16S, basename(names.1.data.Y18.R1))
names.2.cut.16S.Y18.R2 <- file.path(path.2.cut.16S, basename(names.1.data.Y18.R2))

path.2.cut.ITS <- file.path(path.0.project, "2.cut.ITS")
if(!dir.exists(path.2.cut.ITS)) dir.create(path.2.cut.ITS)
names.2.cut.ITS.Y17.R1 <- file.path(path.2.cut.ITS, basename(names.1.data.Y17.R1))
names.2.cut.ITS.Y17.R2 <- file.path(path.2.cut.ITS, basename(names.1.data.Y17.R2))
names.2.cut.ITS.Y18.R1 <- file.path(path.2.cut.ITS, basename(names.1.data.Y18.R1))
names.2.cut.ITS.Y18.R2 <- file.path(path.2.cut.ITS, basename(names.1.data.Y18.R2))

# Use cutadapt to keep only reads with primer match, and to perform other filtering (N's, length criteria)
# Cutadapt version 3.4

# for 2017 bacteria, use primer matching to discard fungal reads
 for(i in seq_along(names.1.data.Y17.R1)) {
   system2(cutadapt, args = c(flags.v5_6.R1, flags.v5_6.R2, 
                             "-n", 2,
                             "--discard-untrimmed", 
                             "--pair-filter=any",
                             "--minimum-length", 50,
                             "--max-n", 0,
                             "-o", shQuote(names.2.cut.16S.Y17.R1[i]), 
                             "-p", shQuote(names.2.cut.16S.Y17.R2[i]), # output file names
                             shQuote(names.1.data.Y17.R1[i]), 
                             shQuote(names.1.data.Y17.R2[i]))) # input file names
}

# for 2018 bacteria, use primer matching + read length to discard fungal reads
for(i in seq_along(names.1.data.Y18.R1)) {
  system2(cutadapt, args = c(flags.v4.R1, flags.v4.R2, 
                             "-n", 2,
                             "--discard-untrimmed", 
                             "--pair-filter=any",
                             "--minimum-length", 50,
                             "--maximum-length", 255, # cull reads from the 2x300 library
                             "--max-n", 0,
                             "-o", shQuote(names.2.cut.16S.Y18.R1[i]), 
                             "-p", shQuote(names.2.cut.16S.Y18.R2[i]), # output file names
                             shQuote(names.1.data.Y18.R1[i]), 
                             shQuote(names.1.data.Y18.R2[i]))) # input file names
}

path.2.cut.ITS <- file.path(path.0.project, "2.cut.ITS")
if(!dir.exists(path.2.cut.ITS)) dir.create(path.2.cut.ITS)
names.2.cut.ITS.Y17.R1 <- file.path(path.2.cut.ITS, basename(names.1.data.Y17.R1))
names.2.cut.ITS.Y17.R2 <- file.path(path.2.cut.ITS, basename(names.1.data.Y17.R2))
names.2.cut.ITS.Y18.R1 <- file.path(path.2.cut.ITS, basename(names.1.data.Y18.R1))
names.2.cut.ITS.Y18.R2 <- file.path(path.2.cut.ITS, basename(names.1.data.Y18.R2))

names.2.cut.ITS.Y17.R1[1]
names.2.cut.ITS.Y17.R2[1]
names.1.data.Y17.R1[1]
names.1.data.Y17.R2[1]

# for 2017 fungi, use primer matching to discard bacterial reads
for(i in seq_along(names.1.data.Y17.R1)) {
  system2(cutadapt, args = c(flags.ITS1.R1, flags.ITS1.R2, 
                             "-n", 2,
                             "--discard-untrimmed", 
                             "--pair-filter=any",
                             "--minimum-length", 50,
                             "--max-n", 0,
                             "-o", shQuote(names.2.cut.ITS.Y17.R1[i]), 
                             "-p", shQuote(names.2.cut.ITS.Y17.R2[i]), # output file names
                             shQuote(names.1.data.Y17.R1[i]), 
                             shQuote(names.1.data.Y17.R2[i]))) # input file names
}

# for 2018 fungi, use primer matching + read length to discard bacterial reads
for(i in seq_along(names.1.data.Y18.R1)) {
  system2(cutadapt, args = c(flags.ITS2.R1, flags.ITS2.R2, 
                             "-n", 2,
                             "--discard-untrimmed", 
                             "--pair-filter=any",
                             "--minimum-length", 255, # cull reads from the 2x250 library
                             "--max-n", 0,
                             "-o", shQuote(names.2.cut.ITS.Y18.R1[i]), 
                             "-p", shQuote(names.2.cut.ITS.Y18.R2[i]), # output file names
                             shQuote(names.1.data.Y18.R1[i]), 
                             shQuote(names.1.data.Y18.R2[i]))) # input file names
}


# Count number of reads in which the primer is found
primerHits <- function(primer, names) {
  nhits <- vcountPattern(primer, sread(readFastq(names)), fixed = FALSE, max.mismatch = 1)
  return(sum(nhits > 0))
}

Y17.16S.R1s <- sapply(c(v4.f, v4.r.rc, v5_6.f, v5_6.r.rc), primerHits, names = names.2.cut.16S.Y17.R1)
Y17.16S.R2s <- sapply(c(v4.f.rc, v4.r, v5_6.f.rc, v5_6.r), primerHits, names = names.2.cut.16S.Y17.R2)
Y18.16S.R1s <- sapply(c(v4.f, v4.r.rc, v5_6.f, v5_6.r.rc), primerHits, names = names.2.cut.16S.Y18.R1)
Y18.16S.R2s <- sapply(c(v4.f.rc, v4.r, v5_6.f.rc, v5_6.r), primerHits, names = names.2.cut.16S.Y18.R2)
Y17.ITS.R1s <- sapply(c(ITS1.f, ITS1.r.rc, ITS2.f, ITS2.r.rc), primerHits, names = names.2.cut.ITS.Y17.R1)
Y17.ITS.R2s <- sapply(c(ITS1.f.rc, ITS1.r, ITS2.f.rc, ITS2.r), primerHits, names = names.2.cut.ITS.Y17.R2)
Y18.ITS.R1s <- sapply(c(ITS1.f, ITS1.r.rc, ITS2.f, ITS2.r.rc), primerHits, names = names.2.cut.ITS.Y18.R1)
Y18.ITS.R2s <- sapply(c(ITS1.f.rc, ITS1.r, ITS2.f.rc, ITS2.r), primerHits, names = names.2.cut.ITS.Y18.R2)

#            GTGCCAGCMGCCGCGGTAA ATTAGAWACCCBDGTAGTCC  AACMGGATTAGATACCCKG     GCAACGAGCGCAACCC 
# Y17.16S.R1s                  0                    0                    0                    0
# Y17.16S.R2s                  0                    0                    0                    0
# Y18.16S.R1s                  0                    0                   39                    0 
# Y18.16S.R2s                  0                    0                   25                    0
#             GCTGCGTTCTTCATCGATGC TTACTTCCTCTAAATGACCAAG     GATGAAGAACGYAGYRAA     CGAGTGAAGMGGBAAVAG 
# Y17.ITS.R1s                    0                      0                      0                      0
# Y17.ITS.R2s                    0                      0                      0                      0
# Y18.ITS.R1s                    0                      0                      0                      0
# Y18.ITS.R2s                    0                      0                      0                      0

### 16S; Y17 samples
## Filter & trim on quality
path.3.filt.16S <- file.path(path.0.project, "3.filt.16S")
if(!dir.exists(path.3.filt.16S)) dir.create(path.3.filt.16S)

names.3.filt.16S.Y17.R1 <- file.path(path.3.filt.16S, basename(names.2.cut.16S.Y17.R1))
names.3.filt.16S.Y17.R2 <- file.path(path.3.filt.16S, basename(names.2.cut.16S.Y17.R2))
sample.no.16S.Y17 <- length(names.3.filt.16S.Y17.R1)

filter.16S.Y17 <- filterAndTrim(names.2.cut.16S.Y17.R1, names.3.filt.16S.Y17.R1, 
                                names.2.cut.16S.Y17.R2, names.3.filt.16S.Y17.R2, 
                                maxN = 0, maxEE = c(2, 2), truncLen = c(225, 225),
                                compress = TRUE)

## Infer sequence variants and merge R1 & R2
err.16S.Y17.R1 <- learnErrors(names.3.filt.16S.Y17.R1, multithread = TRUE, randomize = TRUE)
plotErrors(err.16S.Y17.R1)
err.16S.Y17.R2 <- learnErrors(names.3.filt.16S.Y17.R2, multithread = TRUE, randomize = TRUE)
plotErrors(err.16S.Y17.R2)

derep.16S.Y17.R1 <- derepFastq(names.3.filt.16S.Y17.R1, verbose = TRUE)
names(derep.16S.Y17.R1) <- substr(names(derep.16S.Y17.R1), 1, 10)
derep.16S.Y17.R2 <- derepFastq(names.3.filt.16S.Y17.R2, verbose = TRUE)
names(derep.16S.Y17.R2) <- substr(names(derep.16S.Y17.R2), 1, 10)

dada.16S.Y17.R1 <- dada(derep.16S.Y17.R1, err = err.16S.Y17.R1, multithread = TRUE)
dada.16S.Y17.R1[[4]]
dada.16S.Y17.R2 <- dada(derep.16S.Y17.R2, err = err.16S.Y17.R2, multithread = TRUE)
dada.16S.Y17.R2[[4]]

merge.16S.Y17 <- mergePairs(dada.16S.Y17.R1, derep.16S.Y17.R1, 
                            dada.16S.Y17.R2, derep.16S.Y17.R2, 
                            maxMismatch = 2, verbose = TRUE)
head(merge.16S.Y17[[1]], 15)

# remove large objects that are not needed anymore
rm(derep.16S.Y17.R1)
rm(derep.16S.Y17.R2)

## Make sequence tables
seqtab.16S.Y17 <- makeSequenceTable(merge.16S.Y17)
# 143 samples X 2,259 sequence variants

table(nchar(getSequences(seqtab.16S.Y17)))
# cull unexpectedly short or long reads
seqtab.16S.Y17.length <- seqtab.16S.Y17[, nchar(colnames(seqtab.16S.Y17)) %in% seq(293,303)]
dim(seqtab.16S.Y17.length)
# 143 samples X 2,245 sequence variants

seqtab.16S.Y17.IgnoreEnds <- collapseNoMismatch(seqtab.16S.Y17.length, minOverlap = 200)
dim(seqtab.16S.Y17.IgnoreEnds)
# no change

## Remove chimeras
seqtab.16S.Y17.chimera <- removeBimeraDenovo(seqtab.16S.Y17.length, 
                                             multithread = TRUE, verbose = TRUE, method = "consensus")
dim(seqtab.16S.Y17.chimera)
# 143 samples X 643 sequence variants

1 - (sum(seqtab.16S.Y17.chimera) / sum(seqtab.16S.Y17))
# 0.122 proportion of chimeric reads

## Assign taxonomy
classification.ref.16S <- "C:/Seq_Data/Reference/silva_nr99_v138.1_train_set.fa.gz"
species.ref.16S <- "C:/Seq_Data/Reference/silva_species_assignment_v138.1.fa.gz"
set.seed(1024)
taxtab.16S.Y17.1 <- assignTaxonomy(seqtab.16S.Y17.chimera, refFasta = classification.ref.16S, multithread = TRUE, minBoot = 80)
set.seed(1024)
taxtab.16S.Y17 <- addSpecies(taxtab.16S.Y17.1, refFasta = species.ref.16S, verbose=TRUE)
taxtab.16S.Y17 <- data.frame(taxtab.16S.Y17)
rownames(taxtab.16S.Y17) <- colnames(seqtab.16S.Y17.chimera)

# Make phyloseq objects
bacteria.Y17.ps <- phyloseq(otu_table(seqtab.16S.Y17.chimera, taxa_are_rows = FALSE),
                        tax_table(as.matrix(taxtab.16S.Y17)), sample_data(metadata.Y17))
bacteria.Y17.ps
# 143 samples X 643 taxa

# Taxonomic filtering
table(tax_table(bacteria.Y17.ps)[, "Kingdom"], exclude = NULL)
# all are bacteria

table(tax_table(bacteria.Y17.ps)[, "Phylum"], exclude = NULL)
# 17 are NA

table(tax_table(bacteria.Y17.ps)[, "Order"], exclude = NULL)
# no Chloroplast

table(tax_table(bacteria.Y17.ps)[, "Family"], exclude = NULL)
# 8 are Mitochondria
bacteria.Y17.ps <- subset_taxa(bacteria.Y17.ps, Family != "Mitochondria" | is.na(Family))
# 143 samples X 635 taxa

sort(rowSums(otu_table(bacteria.Y17.ps)))
# range 768 to 33,421 observations per sample

# 'decontam' package for identifying and removing contaminants
sample_data(bacteria.Y17.ps)$Sample_or_Control <- sample_data(bacteria.Y17.ps)$biol
sample_data(bacteria.Y17.ps)$Sample_or_Control <- gsub("FACE", "Sample", sample_data(bacteria.Y17.ps)$Sample_or_Control)
sample_data(bacteria.Y17.ps)$Sample_or_Control <- gsub("Mock", NA, sample_data(bacteria.Y17.ps)$Sample_or_Control)
sample_data(bacteria.Y17.ps)$Sample_or_Control <- gsub("nPCR", "Control", sample_data(bacteria.Y17.ps)$Sample_or_Control)
sample_data(bacteria.Y17.ps)$Sample_or_Control <- gsub("nDNA", "Control", sample_data(bacteria.Y17.ps)$Sample_or_Control)

df.16S.Y17 <- as.data.frame(sample_data(bacteria.Y17.ps))
df.16S.Y17$LibrarySize <- sample_sums(bacteria.Y17.ps)
df.16S.Y17 <- df.16S.Y17[order(df.16S.Y17$LibrarySize),]
df.16S.Y17$Index <- seq(nrow(df.16S.Y17))
ggplot(data=df.16S.Y17, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(bacteria.Y17.ps)$is.neg <- sample_data(bacteria.Y17.ps)$Sample_or_Control == "Control"
bacteria.Y17.NoMock.ps <- prune_samples(!is.na(sample_data(bacteria.Y17.ps)$Sample_or_Control), bacteria.Y17.ps)
#bacteria.Y17.NoMock.ps <- prune_taxa(taxa_sums(bacteria.Y17.NoMock.ps) > 0, bacteria.Y17.NoMock.ps) 
contam.16S.Y17.df <- isContaminant(bacteria.Y17.NoMock.ps, method="prevalence", neg="is.neg")
table(contam.16S.Y17.df$contaminant)
# FALSE  TRUE 
# 609    26 
contam.16S.Y17.check <- otu_table(bacteria.Y17.NoMock.ps)[, contam.16S.Y17.df$contaminant]
# look at the data:
unname(contam.16S.Y17.check)
# remove putative contaminant ASVs
bacteria.Y17.decontam.ps <- prune_taxa(!contam.16S.Y17.df$contaminant, bacteria.Y17.ps)
# 143 samples X 609 taxa
# remove negative control samples
bacteria.Y17.decontam.ps <- prune_samples(sample_data(bacteria.Y17.decontam.ps)$Sample_or_Control != "Control" 
                                          | is.na(sample_data(bacteria.Y17.decontam.ps)$Sample_or_Control), 
                                          bacteria.Y17.decontam.ps)
bacteria.Y17.decontam.ps <- prune_taxa(taxa_sums(bacteria.Y17.decontam.ps) > 0, bacteria.Y17.decontam.ps)
# 139 samples X 588 taxa 

# Pull out mock community controls
bacteria.Y17.mock.ps <- subset_samples(bacteria.Y17.decontam.ps, collect == "Y17_BMock")
bacteria.Y17.mock.ps <- prune_taxa(taxa_sums(bacteria.Y17.mock.ps) > 0, bacteria.Y17.mock.ps) 
write.table(otu_table(bacteria.Y17.mock.ps), file = "bacteria.Y17.mock.OTU_table.txt")
write.table(tax_table(bacteria.Y17.mock.ps), file = "bacteria.Y17.mock.tax_table.txt")

# reduce to only experimental samples
bacteria.Y17.samples.ps <- subset_samples(bacteria.Y17.decontam.ps, Sample_or_Control == "Sample")
bacteria.Y17.samples.ps <- prune_taxa(taxa_sums(bacteria.Y17.samples.ps) > 0, bacteria.Y17.samples.ps) 
# 133 samples X 552 taxa
sort(rowSums(otu_table(bacteria.Y17.samples.ps)))
# minimum size = 1,358 reads
# retain all samples

# construct phylogenetic tree
seqs.16S.Y17 <- colnames(otu_table(bacteria.Y17.samples.ps))
names(seqs.16S.Y17) <- seqs.16S.Y17
alignment.16S.Y17 <- AlignSeqs(DNAStringSet(seqs.16S.Y17), anchor = NA)
alignment.16S.Y17.phang <- phyDat(as(alignment.16S.Y17, "matrix"), type = "DNA")
dm.16S.Y17 <- dist.ml(alignment.16S.Y17.phang)
treeNJ.16S.Y17 <- NJ(dm.16S.Y17)
fit.16S.Y17 <- pml(treeNJ.16S.Y17, data = alignment.16S.Y17.phang)
fit.16S.Y17.GTR <- update(fit.16S.Y17, k = 4, inv = 0.2)
fit.16S.Y17.GTR <- optim.pml(fit.16S.Y17.GTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

# merge tree into phyloseq object
bacteria.Y17.samples.ps = merge_phyloseq(bacteria.Y17.samples.ps, phy_tree(fit.16S.Y17.GTR$tree))
# 133 samples X 552 taxa


### 16S; Y18 samples
## Filter & trim on quality
names.3.filt.16S.Y18.R1 <- file.path(path.3.filt.16S, basename(names.2.cut.16S.Y18.R1))
names.3.filt.16S.Y18.R2 <- file.path(path.3.filt.16S, basename(names.2.cut.16S.Y18.R2))
sample.no.16S.Y18 <- length(names.3.filt.16S.Y18.R1)
sample.no.16S.Y18
# number of samples = 135

filter.16S.Y18 <- filterAndTrim(names.2.cut.16S.Y18.R1, names.3.filt.16S.Y18.R1, 
                                names.2.cut.16S.Y18.R2, names.3.filt.16S.Y18.R2, 
                                maxN = 0, maxEE = c(2, 2), truncLen = c(225, 225),
                                compress = TRUE)

# some files had no reads pass the filter; identify and exclude these.
# SRR9124623, SRR9124844
filter.16S.Y18
grep("SRR9124623", names.3.filt.16S.Y18.R1) #12
grep("SRR9124844", names.3.filt.16S.Y18.R1) #120
names.3.filt.16S.Y18.R1 <- names.3.filt.16S.Y18.R1[-c(12, 120)]
names.3.filt.16S.Y18.R2 <- names.3.filt.16S.Y18.R2[-c(12, 120)]

## Infer sequence variants and merge R1 & R2
err.16S.Y18.R1 <- learnErrors(names.3.filt.16S.Y18.R1, multithread = TRUE, randomize = TRUE)
plotErrors(err.16S.Y18.R1)
err.16S.Y18.R2 <- learnErrors(names.3.filt.16S.Y18.R2, multithread = TRUE, randomize = TRUE)
plotErrors(err.16S.Y18.R2)

derep.16S.Y18.R1 <- derepFastq(names.3.filt.16S.Y18.R1, verbose = TRUE)
names(derep.16S.Y18.R1) <- substr(names(derep.16S.Y18.R1), 1, 10)
derep.16S.Y18.R2 <- derepFastq(names.3.filt.16S.Y18.R2, verbose = TRUE)
names(derep.16S.Y18.R2) <- substr(names(derep.16S.Y18.R2), 1, 10)

dada.16S.Y18.R1 <- dada(derep.16S.Y18.R1, err = err.16S.Y18.R1, multithread = TRUE)
dada.16S.Y18.R1[[4]]
dada.16S.Y18.R2 <- dada(derep.16S.Y18.R2, err = err.16S.Y18.R2, multithread = TRUE)
dada.16S.Y18.R2[[4]]

merge.16S.Y18 <- mergePairs(dada.16S.Y18.R1, derep.16S.Y18.R1, 
                            dada.16S.Y18.R2, derep.16S.Y18.R2, 
                            maxMismatch = 2, verbose = TRUE)
head(merge.16S.Y18[[1]])

# remove large files that are not needed anymore
rm(derep.16S.Y18.R1)
rm(derep.16S.Y18.R2)

## Make sequence tables
seqtab.16S.Y18 <- makeSequenceTable(merge.16S.Y18)
# 133 samples X 3,413 sequence variants

table(nchar(getSequences(seqtab.16S.Y18)))
# cull unexpectedly short or long reads
seqtab.16S.Y18.length <- seqtab.16S.Y18[, nchar(colnames(seqtab.16S.Y18)) %in% seq(247,259)]
# 133 samples X 3,093 sequence variants

seqtab.16S.Y18.IgnoreEnds <- collapseNoMismatch(seqtab.16S.Y18.length, minOverlap = 200)
dim(seqtab.16S.Y18.IgnoreEnds)
# no change

## Remove chimeras
seqtab.16S.Y18.chimera <- removeBimeraDenovo(seqtab.16S.Y18.length, 
                                             multithread = TRUE, verbose = TRUE, method = "consensus")
dim(seqtab.16S.Y18.chimera)
# 133 samples X 944 sequence variants
1 - (sum(seqtab.16S.Y18.chimera) / sum(seqtab.16S.Y18))
# 0.096 proportion of reads that look chimeric

## Assign taxonomy
set.seed(1024)
taxtab.16S.Y18.1 <- assignTaxonomy(seqtab.16S.Y18.chimera, refFasta = classification.ref.16S, multithread = TRUE, minBoot = 80)
set.seed(1024)
taxtab.16S.Y18 <- addSpecies(taxtab.16S.Y18.1, refFasta = species.ref.16S, verbose=TRUE)
dim(taxtab.16S.Y18)
# 944 sequence variants
taxtab.16S.Y18 <- data.frame(taxtab.16S.Y18)
rownames(taxtab.16S.Y18) <- colnames(seqtab.16S.Y18.chimera)

# Make phyloseq objects
bacteria.Y18.ps <- phyloseq(otu_table(seqtab.16S.Y18.chimera, taxa_are_rows = FALSE),
                            tax_table(as.matrix(taxtab.16S.Y18)), sample_data(metadata.Y18))
bacteria.Y18.ps
# 133 samples X 944 taxa

# Taxonomic filtering
table(tax_table(bacteria.Y18.ps)[, "Kingdom"], exclude = NULL)
# all are bacteria

table(tax_table(bacteria.Y18.ps)[, "Phylum"], exclude = NULL)
# 58 are NA

table(tax_table(bacteria.Y18.ps)[, "Order"], exclude = NULL)
# 40 are Chloroplast
bacteria.Y18.ps <- subset_taxa(bacteria.Y18.ps, Order != "Chloroplast" | is.na(Order))

table(tax_table(bacteria.Y18.ps)[, "Family"], exclude = NULL)
# 153 are Mitochondria
bacteria.Y18.ps <- subset_taxa(bacteria.Y18.ps, Family != "Mitochondria" | is.na(Family))
# 133 samples X 751 taxa

sort(rowSums(otu_table(bacteria.Y18.ps)))
# range 0 to 31,928 observations per sample

# 'decontam' package for identifying and removing contaminants
sample_data(bacteria.Y18.ps)$Sample_or_Control <- sample_data(bacteria.Y18.ps)$biol
sample_data(bacteria.Y18.ps)$Sample_or_Control <- gsub("FACE", "Sample", sample_data(bacteria.Y18.ps)$Sample_or_Control)
sample_data(bacteria.Y18.ps)$Sample_or_Control <- gsub("Mock", NA, sample_data(bacteria.Y18.ps)$Sample_or_Control)
sample_data(bacteria.Y18.ps)$Sample_or_Control <- gsub("nPCR", "Control", sample_data(bacteria.Y18.ps)$Sample_or_Control)
sample_data(bacteria.Y18.ps)$Sample_or_Control <- gsub("nDNA", "Control", sample_data(bacteria.Y18.ps)$Sample_or_Control)

df.16S.Y18 <- as.data.frame(sample_data(bacteria.Y18.ps))
df.16S.Y18$LibrarySize <- sample_sums(bacteria.Y18.ps)
df.16S.Y18 <- df.16S.Y18[order(df.16S.Y18$LibrarySize),]
df.16S.Y18$Index <- seq(nrow(df.16S.Y18))
ggplot(data=df.16S.Y18, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(bacteria.Y18.ps)$is.neg <- sample_data(bacteria.Y18.ps)$Sample_or_Control == "Control"
bacteria.Y18.NoMock.ps <- prune_samples(!is.na(sample_data(bacteria.Y18.ps)$Sample_or_Control), bacteria.Y18.ps)
contam.16S.Y18.df <- isContaminant(bacteria.Y18.NoMock.ps, method="prevalence", neg="is.neg")
table(contam.16S.Y18.df$contaminant)
# FALSE  TRUE 
# 736    15 
contam.16S.Y18.check <- otu_table(bacteria.Y18.NoMock.ps)[, contam.16S.Y18.df$contaminant]
# remove putative contaminant ASVs
bacteria.Y18.decontam.ps <- prune_taxa(!contam.16S.Y18.df$contaminant, bacteria.Y18.ps)
# 133 samples X 736 taxa
# remove negative control samples
bacteria.Y18.decontam.ps <- prune_samples(sample_data(bacteria.Y18.decontam.ps)$Sample_or_Control != "Control" 
                                          | is.na(sample_data(bacteria.Y18.decontam.ps)$Sample_or_Control), 
                                          bacteria.Y18.decontam.ps)
bacteria.Y18.decontam.ps <- prune_taxa(taxa_sums(bacteria.Y18.decontam.ps) > 0, bacteria.Y18.decontam.ps)
# 130 samples X 692 taxa 


# Pull out mock community controls
bacteria.Y18.mock.ps <- subset_samples(bacteria.Y18.decontam.ps, collect == "Y18_BMock")
bacteria.Y18.mock.ps <- prune_taxa(taxa_sums(bacteria.Y18.mock.ps) > 0, bacteria.Y18.mock.ps) 
sample_data(bacteria.Y18.mock.ps)
unname(otu_table(bacteria.Y18.mock.ps))
unname(tax_table(bacteria.Y18.mock.ps))
write.table(otu_table(bacteria.Y18.mock.ps), file = "bacteria.Y18.mock.OTU_table.txt")
write.table(tax_table(bacteria.Y18.mock.ps), file = "bacteria.Y18.mock.tax_table.txt")


# reduce to only experimental samples
bacteria.Y18.samples.ps <- subset_samples(bacteria.Y18.decontam.ps, Sample_or_Control == "Sample")
bacteria.Y18.samples.ps <- prune_taxa(taxa_sums(bacteria.Y18.samples.ps) > 0, bacteria.Y18.samples.ps) 
# 128 samples X 654 taxa
sort(rowSums(otu_table(bacteria.Y18.samples.ps)))
# minimum size = 0 reads
# cull 5 samples with < 1200 reads
bacteria.Y18.samples.ps <- prune_samples(sample_sums(bacteria.Y18.samples.ps) >= 1200, bacteria.Y18.samples.ps)
bacteria.Y18.samples.ps <- prune_taxa(taxa_sums(bacteria.Y18.samples.ps) > 0, bacteria.Y18.samples.ps) 
# 123 samples X 654 taxa

# construct phylogenetic tree
seqs.16S.Y18 <- colnames(otu_table(bacteria.Y18.samples.ps))
names(seqs.16S.Y18) <- seqs.16S.Y18
alignment.16S.Y18 <- AlignSeqs(DNAStringSet(seqs.16S.Y18), anchor = NA)
alignment.16S.Y18.phang <- phyDat(as(alignment.16S.Y18, "matrix"), type = "DNA")
dm.16S.Y18 <- dist.ml(alignment.16S.Y18.phang)
treeNJ.16S.Y18 <- NJ(dm.16S.Y18)
fit.16S.Y18 <- pml(treeNJ.16S.Y18, data = alignment.16S.Y18.phang)
fit.16S.Y18.GTR <- update(fit.16S.Y18, k = 4, inv = 0.2)
fit.16S.Y18.GTR <- optim.pml(fit.16S.Y18.GTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                             rearrangement = "stochastic", control = pml.control(trace = 0))

# merge tree into phyloseq object
bacteria.Y18.samples.ps = merge_phyloseq(bacteria.Y18.samples.ps, phy_tree(fit.16S.Y18.GTR$tree))
# 123 samples X 654 taxa


### ITS; Y17 samples
## Filter & trim on quality
path.3.filt.ITS <- file.path(path.0.project, "3.filt.ITS")
if(!dir.exists(path.3.filt.ITS)) dir.create(path.3.filt.ITS)

names.3.filt.ITS.Y17.R1 <- file.path(path.3.filt.ITS, basename(names.2.cut.ITS.Y17.R1))
names.3.filt.ITS.Y17.R2 <- file.path(path.3.filt.ITS, basename(names.2.cut.ITS.Y17.R2))
sample.no.ITS.Y17 <- length(names.3.filt.ITS.Y17.R1)
# number of samples = 143

filter.ITS.Y17 <- filterAndTrim(names.2.cut.ITS.Y17.R1, names.3.filt.ITS.Y17.R1, 
                                names.2.cut.ITS.Y17.R2, names.3.filt.ITS.Y17.R2, 
                                maxN = 0, maxEE = c(2, 2), truncLen = c(225, 225), # Note, this set used a v2 sequencing kit
                                compress = TRUE)

## Infer sequence variants and merge R1 & R2
err.ITS.Y17.R1 <- learnErrors(names.3.filt.ITS.Y17.R1, multithread = TRUE, randomize = TRUE)
plotErrors(err.ITS.Y17.R1)
err.ITS.Y17.R2 <- learnErrors(names.3.filt.ITS.Y17.R2, multithread = TRUE, randomize = TRUE)
plotErrors(err.ITS.Y17.R2)

derep.ITS.Y17.R1 <- derepFastq(names.3.filt.ITS.Y17.R1, verbose = TRUE)
names(derep.ITS.Y17.R1) <- substr(names(derep.ITS.Y17.R1), 1, 10)
derep.ITS.Y17.R2 <- derepFastq(names.3.filt.ITS.Y17.R2, verbose = TRUE)
names(derep.ITS.Y17.R2) <- substr(names(derep.ITS.Y17.R2), 1, 10)

dada.ITS.Y17.R1 <- dada(derep.ITS.Y17.R1, err = err.ITS.Y17.R1, multithread = TRUE)
dada.ITS.Y17.R1[[4]]
dada.ITS.Y17.R2 <- dada(derep.ITS.Y17.R2, err = err.ITS.Y17.R2, multithread = TRUE)
dada.ITS.Y17.R2[[4]]

merge.ITS.Y17 <- mergePairs(dada.ITS.Y17.R1, derep.ITS.Y17.R1, 
                            dada.ITS.Y17.R2, derep.ITS.Y17.R2, 
                            maxMismatch = 2, verbose = TRUE)
head(merge.ITS.Y17[[1]])

# remove large files that are not needed anymore
rm(derep.ITS.Y17.R1)
rm(derep.ITS.Y17.R2)

## Make sequence tables
seqtab.ITS.Y17 <- makeSequenceTable(merge.ITS.Y17)
dim(seqtab.ITS.Y17)
# 143 samples X 483 sequence variants
unname(seqtab.ITS.Y17[, 1:10])

table(nchar(getSequences(seqtab.ITS.Y17)))
# sequence variants range in length from 225 to 438 bases

seqtab.ITS.Y17.IgnoreEnds <- collapseNoMismatch(seqtab.ITS.Y17, minOverlap = 175)
dim(seqtab.ITS.Y17.IgnoreEnds)
# merges 3 ASVs... not a substantial effect

## Remove chimeras
seqtab.ITS.Y17.chimera <- removeBimeraDenovo(seqtab.ITS.Y17, 
                                             multithread = TRUE, verbose = TRUE, method = "consensus")
dim(seqtab.ITS.Y17.chimera)
# 143 samples X 393 sequence variants

1 - (sum(seqtab.ITS.Y17.chimera) / sum(seqtab.ITS.Y17))
# 0.012 proportion of reads look chimeric

## Trim off conserved flanking regions
path.4.ITSx <- file.path(path.0.project, "4.ITSx")
if(!dir.exists(path.4.ITSx)) dir.create(path.4.ITSx)
Y17.ITSx.input <- DNAStringSet(colnames(seqtab.ITS.Y17.chimera))
Y17.ASV.names <- seq(dim(seqtab.ITS.Y17.chimera)[2])
names(Y17.ITSx.input) <- Y17.ASV.names
writeXStringSet(Y17.ITSx.input, file.path(path.4.ITSx, "Y17.ITSx.input.fasta"), format = "fasta")

# in Ubuntu:
# ITSx --license
# v. 1.1.3

# permit ITSx to flag reads as belonging to other than fungi:
# ITSx -i Y17.ITSx.input.fasta -o Y17.ITSx_out --allow_reorder T --allow_single_domain 1e-3,0 --partial 50
# output suggests many non-fungal taxa present

# pull out just the fungal reads:
# grep -A 1 "|F" Y17.ITSx_out.ITS1.full_and_partial.fasta > Y17.ITSx.fungi.fasta
# sed -i 's/|/ /' Y17.ITSx.fungi.fasta
# sed -i 's/\s.*$//' Y17.ITSx.fungi.fasta
# sed -i 's/--//' Y17.ITSx.fungi.fasta

fungi.Y17 <- readDNAStringSet(file.path(path.4.ITSx, "Y17.ITSx.fungi.fasta"))
width(fungi.Y17)
length(Y17.ITSx.input) - length(fungi.Y17)
# 198 of the original sequence variants were not recognizable as fungal ITS reads


## Assign taxonomy
classification.ref.ITS1 <- "C:/Seq_Data/Reference/UNITE_8.2_RefSingletons_ITSx.ITS1.full_and_partial.fasta"
set.seed(1024)
taxtab.ITS.Y17 <- assignTaxonomy(fungi.Y17, refFasta = classification.ref.ITS1, multithread = TRUE, minBoot = 80)
dim(taxtab.ITS.Y17)
# 195 sequence variants
taxtab.ITS.Y17 <- data.frame(taxtab.ITS.Y17)

fungi.Y17.df <- data.frame(fungi.Y17)
# row names tell which columns of seqtab.chimera should be retained, and which should be dropped as non-fungal / non-ITS2
dim(seqtab.ITS.Y17.chimera)
# 143 X 393
seqtab.ITS.Y17.fungi <- seqtab.ITS.Y17.chimera[, as.numeric(rownames(fungi.Y17.df))]
dim(seqtab.ITS.Y17.fungi)
# 143 X 195

# problem using actual sequences as ASV identifiers, because post-ITSx there are duplicate sequences...
sum(duplicated(fungi.Y17.df$fungi.Y17))
# 2
duplicated(fungi.Y17.df$fungi.Y17)
# 60 & 107 have duplicates
fungi.Y17.df$fungi.Y17[60]
fungi.Y17.df$fungi.Y17[107]
# Keep using the pre-ITSx sequences instead
rownames(taxtab.ITS.Y17) <- colnames(seqtab.ITS.Y17.fungi)

# Make phyloseq objects
fungi.Y17.ps <- phyloseq(otu_table(seqtab.ITS.Y17.fungi, taxa_are_rows = FALSE),
                         tax_table(as.matrix(taxtab.ITS.Y17)), sample_data(metadata.Y17))
fungi.Y17.ps
# 143 samples X 195 taxa


write.table(fungi.Y17, file = "ASVs_Y17_postITSx.txt")
write.table(rownames(tax_table(fungi.Y17.ps)), file = "ASVs_fungi.ps.txt")

# Taxonomic filtering
table(tax_table(fungi.Y17.ps)[, "Kingdom"], exclude = NULL)
# all are fungi

table(tax_table(fungi.Y17.ps)[, "Phylum"], exclude = NULL)
# 3 are NA

# 'decontam' package for identifying and removing contaminants
sample_data(fungi.Y17.ps)$Sample_or_Control <- sample_data(fungi.Y17.ps)$biol
sample_data(fungi.Y17.ps)$Sample_or_Control <- gsub("FACE", "Sample", sample_data(fungi.Y17.ps)$Sample_or_Control)
sample_data(fungi.Y17.ps)$Sample_or_Control <- gsub("Mock", NA, sample_data(fungi.Y17.ps)$Sample_or_Control)
sample_data(fungi.Y17.ps)$Sample_or_Control <- gsub("nPCR", "Control", sample_data(fungi.Y17.ps)$Sample_or_Control)
sample_data(fungi.Y17.ps)$Sample_or_Control <- gsub("nDNA", "Control", sample_data(fungi.Y17.ps)$Sample_or_Control)

df.ITS.Y17 <- as.data.frame(sample_data(fungi.Y17.ps))
df.ITS.Y17$LibrarySize <- sample_sums(fungi.Y17.ps)
df.ITS.Y17 <- df.ITS.Y17[order(df.ITS.Y17$LibrarySize),]
df.ITS.Y17$Index <- seq(nrow(df.ITS.Y17))
ggplot(data=df.ITS.Y17, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(fungi.Y17.ps)$is.neg <- sample_data(fungi.Y17.ps)$Sample_or_Control == "Control"
fungi.Y17.NoMock.ps <- prune_samples(!is.na(sample_data(fungi.Y17.ps)$Sample_or_Control), fungi.Y17.ps)
#fungi.Y17.NoMock.ps <- prune_taxa(taxa_sums(fungi.Y17.NoMock.ps) > 0, fungi.Y17.NoMock.ps) 
contam.ITS.Y17.df <- isContaminant(fungi.Y17.NoMock.ps, method="prevalence", neg="is.neg")
table(contam.ITS.Y17.df$contaminant)
# FALSE  TRUE 
# 194    1 
contam.ITS.Y17.check <- otu_table(fungi.Y17.NoMock.ps)[, contam.ITS.Y17.df$contaminant]
# remove putative contaminant ASVs
fungi.Y17.decontam.ps <- prune_taxa(!contam.ITS.Y17.df$contaminant, fungi.Y17.ps)
# 143 samples X 194 taxa
# remove negative control samples
fungi.Y17.decontam.ps <- prune_samples(sample_data(fungi.Y17.decontam.ps)$Sample_or_Control != "Control" 
                                          | is.na(sample_data(fungi.Y17.decontam.ps)$Sample_or_Control), 
                                          fungi.Y17.decontam.ps)
fungi.Y17.decontam.ps <- prune_taxa(taxa_sums(fungi.Y17.decontam.ps) > 0, fungi.Y17.decontam.ps)
# 139 samples X 193 taxa 


# Pull out mock community controls
fungi.Y17.mock.ps <- subset_samples(fungi.Y17.decontam.ps, collect == "Y17_FMock")
fungi.Y17.mock.ps <- prune_taxa(taxa_sums(fungi.Y17.mock.ps) > 0, fungi.Y17.mock.ps) 
sample_data(fungi.Y17.mock.ps)
unname(otu_table(fungi.Y17.mock.ps))
unname(tax_table(fungi.Y17.mock.ps))
write.table(otu_table(fungi.Y17.mock.ps), file = "fungi.Y17.mock.OTU_table.txt")
write.table(tax_table(fungi.Y17.mock.ps), file = "fungi.Y17.mock.tax_table.txt")


# reduce to only experimental samples
fungi.Y17.samples.ps <- subset_samples(fungi.Y17.decontam.ps, Sample_or_Control == "Sample")
fungi.Y17.samples.ps <- prune_taxa(taxa_sums(fungi.Y17.samples.ps) > 0, fungi.Y17.samples.ps) 
# 133 samples X 185 taxa
sort(rowSums(otu_table(fungi.Y17.samples.ps)))
# minimum size = 6 reads
# cull 33 samples with < 495 reads
fungi.Y17.samples.ps <- prune_samples(sample_sums(fungi.Y17.samples.ps) >= 495, fungi.Y17.samples.ps)
fungi.Y17.samples.ps <- prune_taxa(taxa_sums(fungi.Y17.samples.ps) > 0, fungi.Y17.samples.ps) 
# 100 samples X 165 taxa



### ITS; Y18 samples
## Filter & trim on quality
names.3.filt.ITS.Y18.R1 <- file.path(path.3.filt.ITS, basename(names.2.cut.ITS.Y18.R1))
names.3.filt.ITS.Y18.R2 <- file.path(path.3.filt.ITS, basename(names.2.cut.ITS.Y18.R2))
sample.no.ITS.Y18 <- length(names.3.filt.ITS.Y18.R1)
# number of samples = 135

filter.ITS.Y18 <- filterAndTrim(names.2.cut.ITS.Y18.R1, names.3.filt.ITS.Y18.R1, 
                                names.2.cut.ITS.Y18.R2, names.3.filt.ITS.Y18.R2, 
                                maxN = 0, maxEE = c(2, 2), truncLen = c(280, 250),
                                compress = TRUE)
# some files had no reads pass the filter; identify and exclude these.
# SRR9124860, SRR9124861
grep("SRR9124860", names.3.filt.ITS.Y18.R1) #130
grep("SRR9124861", names.3.filt.ITS.Y18.R1) #131
names.3.filt.ITS.Y18.R1 <- names.3.filt.ITS.Y18.R1[-c(130, 131)]
names.3.filt.ITS.Y18.R2 <- names.3.filt.ITS.Y18.R2[-c(130, 131)]

## Infer sequence variants and merge R1 & R2
err.ITS.Y18.R1 <- learnErrors(names.3.filt.ITS.Y18.R1, multithread = TRUE, randomize = TRUE)
plotErrors(err.ITS.Y18.R1)
err.ITS.Y18.R2 <- learnErrors(names.3.filt.ITS.Y18.R2, multithread = TRUE, randomize = TRUE)
plotErrors(err.ITS.Y18.R2)

derep.ITS.Y18.R1 <- derepFastq(names.3.filt.ITS.Y18.R1, verbose = TRUE)
names(derep.ITS.Y18.R1) <- substr(names(derep.ITS.Y18.R1), 1, 10)
derep.ITS.Y18.R2 <- derepFastq(names.3.filt.ITS.Y18.R2, verbose = TRUE)
names(derep.ITS.Y18.R2) <- substr(names(derep.ITS.Y18.R2), 1, 10)

dada.ITS.Y18.R1 <- dada(derep.ITS.Y18.R1, err = err.ITS.Y18.R1, multithread = TRUE)
dada.ITS.Y18.R1[[4]]
dada.ITS.Y18.R2 <- dada(derep.ITS.Y18.R2, err = err.ITS.Y18.R2, multithread = TRUE)
dada.ITS.Y18.R2[[4]]

merge.ITS.Y18 <- mergePairs(dada.ITS.Y18.R1, derep.ITS.Y18.R1, 
                            dada.ITS.Y18.R2, derep.ITS.Y18.R2, 
                            maxMismatch = 2, verbose = TRUE)
head(merge.ITS.Y18[[1]])

# remove large files that are not needed anymore
rm(derep.ITS.Y18.R1)
rm(derep.ITS.Y18.R2)

## Make sequence tables
seqtab.ITS.Y18 <- makeSequenceTable(merge.ITS.Y18)
dim(seqtab.ITS.Y18)
# 133 samples X 3,954 sequence variants
unname(seqtab.ITS.Y18[, 1:10])

table(nchar(getSequences(seqtab.ITS.Y18)))
# sequence variants range in length from 285 to 518 bases

seqtab.ITS.Y18.IgnoreEnds <- collapseNoMismatch(seqtab.ITS.Y18, minOverlap = 175)
dim(seqtab.ITS.Y18.IgnoreEnds)
# merges 1 ASV... not a substantial effect

## Remove chimeras
seqtab.ITS.Y18.chimera <- removeBimeraDenovo(seqtab.ITS.Y18, 
                                             multithread = TRUE, verbose = TRUE, method = "consensus")
dim(seqtab.ITS.Y18.chimera)
# 133 samples X 1,689 sequence variants

1 - (sum(seqtab.ITS.Y18.chimera) / sum(seqtab.ITS.Y18))
# 0.094 proportion of ASVs that appear to be chimeric

## Trim off conserved flanking regions
Y18.ITSx.input <- DNAStringSet(colnames(seqtab.ITS.Y18.chimera))
Y18.ASV.names <- seq(dim(seqtab.ITS.Y18.chimera)[2])
names(Y18.ITSx.input) <- Y18.ASV.names
writeXStringSet(Y18.ITSx.input, file.path(path.4.ITSx, "Y18.ITSx.input.fasta"), format = "fasta")

# in Ubuntu:
# permit ITSx to flag reads as belonging to other than fungi:
# ITSx -i Y18.ITSx.input.fasta -o Y18.ITSx_out --allow_reorder T --allow_single_domain 1e-3,0 --partial 50
# output suggests many non-fungal taxa present

# pull out just the fungal reads:
# grep -A 1 "|F" Y18.ITSx_out.ITS2.full_and_partial.fasta > Y18.ITSx.fungi.fasta
# sed -i 's/|/ /' Y18.ITSx.fungi.fasta
# sed -i 's/\s.*$//' Y18.ITSx.fungi.fasta
# sed -i 's/--//' Y18.ITSx.fungi.fasta

fungi.Y18 <- readDNAStringSet(file.path(path.4.ITSx, "Y18.ITSx.fungi.fasta"))
width(fungi.Y18)
length(Y18.ITSx.input) - length(fungi.Y18)
# 299 of the original sequence variants were not recognizable as fungal ITS reads

## Assign taxonomy
classification.ref.ITS2 <- "C:/Seq_Data/Reference/UNITE_8.2_RefSingletons_ITSx.ITS2.full_and_partial.fasta"
set.seed(1024)
taxtab.ITS.Y18 <- assignTaxonomy(fungi.Y18, refFasta = classification.ref.ITS2, multithread = TRUE, minBoot = 80)
dim(taxtab.ITS.Y18)
# 1,390 sequence variants
taxtab.ITS.Y18 <- data.frame(taxtab.ITS.Y18)

fungi.Y18.df <- data.frame(fungi.Y18)
# row names tell which columns of seqtab.chimera should be retained, and which should be dropped as non-fungal / non-ITS2
dim(seqtab.ITS.Y18.chimera)
# 133 X 1,689
seqtab.ITS.Y18.fungi <- seqtab.ITS.Y18.chimera[, as.numeric(rownames(fungi.Y18.df))]
dim(seqtab.ITS.Y18.fungi)
# 133 X 1,390

# problem using actual sequences as ASV identifiers, because post-ITSx there are duplicate sequences...
sum(duplicated(fungi.Y18.df$fungi.Y18))
# 162
# Keep using the pre-ITSx sequences instead
rownames(taxtab.ITS.Y18) <- colnames(seqtab.ITS.Y18.fungi)

# Make phyloseq objects
fungi.Y18.ps <- phyloseq(otu_table(seqtab.ITS.Y18.fungi, taxa_are_rows = FALSE),
                         tax_table(as.matrix(taxtab.ITS.Y18)), sample_data(metadata.Y18))
fungi.Y18.ps
# 133 samples X 1,390 taxa

# Taxonomic filtering
table(tax_table(fungi.Y18.ps)[, "Kingdom"], exclude = NULL)
# all are fungi

table(tax_table(fungi.Y18.ps)[, "Phylum"], exclude = NULL)
# 170 are NA

# 'decontam' package for identifying and removing contaminants
sample_data(fungi.Y18.ps)$Sample_or_Control <- sample_data(fungi.Y18.ps)$biol
sample_data(fungi.Y18.ps)$Sample_or_Control <- gsub("FACE", "Sample", sample_data(fungi.Y18.ps)$Sample_or_Control)
sample_data(fungi.Y18.ps)$Sample_or_Control <- gsub("Mock", NA, sample_data(fungi.Y18.ps)$Sample_or_Control)
sample_data(fungi.Y18.ps)$Sample_or_Control <- gsub("nPCR", "Control", sample_data(fungi.Y18.ps)$Sample_or_Control)
sample_data(fungi.Y18.ps)$Sample_or_Control <- gsub("nDNA", "Control", sample_data(fungi.Y18.ps)$Sample_or_Control)

df.ITS.Y18 <- as.data.frame(sample_data(fungi.Y18.ps))
df.ITS.Y18$LibrarySize <- sample_sums(fungi.Y18.ps)
df.ITS.Y18 <- df.ITS.Y18[order(df.ITS.Y18$LibrarySize),]
df.ITS.Y18$Index <- seq(nrow(df.ITS.Y18))
ggplot(data = df.ITS.Y18, aes(x = Index, y = LibrarySize, color = Sample_or_Control)) + geom_point()

sample_data(fungi.Y18.ps)$is.neg <- sample_data(fungi.Y18.ps)$Sample_or_Control == "Control"
fungi.Y18.NoMock.ps <- prune_samples(!is.na(sample_data(fungi.Y18.ps)$Sample_or_Control), fungi.Y18.ps)
#fungi.Y18.NoMock.ps <- prune_taxa(taxa_sums(fungi.Y18.NoMock.ps) > 0, fungi.Y18.NoMock.ps) 
contam.ITS.Y18.df <- isContaminant(fungi.Y18.NoMock.ps, method="prevalence", neg="is.neg")
table(contam.ITS.Y18.df$contaminant)
# FALSE  TRUE 
# 1390   0 
# remove negative control samples
fungi.Y18.decontam.ps <- prune_samples(sample_data(fungi.Y18.ps)$Sample_or_Control != "Control" 
                                       | is.na(sample_data(fungi.Y18.ps)$Sample_or_Control), 
                                       fungi.Y18.ps)
fungi.Y18.decontam.ps <- prune_taxa(taxa_sums(fungi.Y18.decontam.ps) > 0, fungi.Y18.decontam.ps)
# 130 samples X 1,387 taxa 
# lose 3 ASVs present only in the negative controls


# Pull out mock community controls
fungi.Y18.mock.ps <- subset_samples(fungi.Y18.decontam.ps, collect == "Y18_FMock")
fungi.Y18.mock.ps <- prune_taxa(taxa_sums(fungi.Y18.mock.ps) > 0, fungi.Y18.mock.ps) 
sample_data(fungi.Y18.mock.ps)
unname(otu_table(fungi.Y18.mock.ps))
unname(tax_table(fungi.Y18.mock.ps))
write.table(otu_table(fungi.Y18.mock.ps), file = "fungi.Y18.mock.OTU_table.txt")
write.table(tax_table(fungi.Y18.mock.ps), file = "fungi.Y18.mock.tax_table.txt")


# reduce to only experimental samples
fungi.Y18.samples.ps <- subset_samples(fungi.Y18.decontam.ps, Sample_or_Control == "Sample")
fungi.Y18.samples.ps <- prune_taxa(taxa_sums(fungi.Y18.samples.ps) > 0, fungi.Y18.samples.ps) 
# 128 samples X 1,375 taxa
sort(rowSums(otu_table(fungi.Y18.samples.ps)))
# minimum size = 12,132 reads
# retain all samples


### The above has produced processed datasets for each year & locus
# bacteria.Y17.samples.ps
# bacteria.Y18.samples.ps
# fungi.Y17.samples.ps
# fungi.Y18.samples.ps

save(bacteria.Y17.samples.ps, file = "FACE_unfiltered_b17.ps.RData")
write.csv(otu_table(bacteria.Y17.samples.ps), file = "bacteria.Y17.ASVs.csv")
write.csv(tax_table(bacteria.Y17.samples.ps), file = "bacteria.Y17.names.csv")
save(bacteria.Y18.samples.ps, file = "FACE_unfiltered_b18.ps.RData")
write.csv(otu_table(bacteria.Y18.samples.ps), file = "bacteria.Y18.ASVs.csv")
write.csv(tax_table(bacteria.Y18.samples.ps), file = "bacteria.Y18.names.csv")
save(fungi.Y17.samples.ps, file = "FACE_unfiltered_f17.ps.RData")
write.csv(otu_table(fungi.Y17.samples.ps), file = "fungi.Y17.ASVs.csv")
write.csv(tax_table(fungi.Y17.samples.ps), file = "fungi.Y17.names.csv")
save(fungi.Y18.samples.ps, file = "FACE_unfiltered_f18.ps.RData")
write.csv(otu_table(fungi.Y18.samples.ps), file = "fungi.Y18.ASVs.csv")
write.csv(tax_table(fungi.Y18.samples.ps), file = "fungi.Y18.names.csv")



