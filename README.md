# Greer_2021
Transcriptome and ribosome profiling analysis in C. elegans  

**Prerequisites:**  
[cutadapt 1.18](https://cutadapt.readthedocs.io/en/stable/index.html)  
[STAR-2.7.2b](https://github.com/alexdobin/STAR)  
[gffread utility](http://ccb.jhu.edu/software/stringtie/gff.shtml)  
Transcriptome samples were sequenced in paired-end 150 nt mode on Illumina sequencer.
Ribosome profiling samples were prepared with Illumina Small RNA TruSeq kit and sequenced in single-end 50 nt mode on Illumina sequencer.
Raw sequencing files are available from [GEO]().


### Preparing genome annotation and index files
C. elegans genomic sequences and annotation files (WS268) were downloaded from the [Wormbase](https://wormbase.org/).

| files                                       | MD5 check sum (unzipped)         | Description              |
| ------------------------------------------- |:--------------------------------:| -------------------------|
| c_elegans.PRJNA13758.WS268.annotations.gff3 | 2b353175bf6e8410815aede3a77a8a62 | annotation               |
| c_elegans.PRJNA13758.WS268.genomic.fa       | d570defcdc006a7c2859fc92dbb21bc4 | Genome sequence          |

<details><summary><b>Prepare custom genomic annotation</b></summary>
Keep only 'Wormbase' feature types for C. elegans (manually curated). Discard other types (usually predicted or related to other nematode species). Drop annotation of non-coding features such as miRNA and pseudogenes.  
     
```R
library(data.table)
library(magrittr)
library(rstudioapi)
library(stringr)
setwd(dirname(getActiveDocumentContext()$path))

#------------------------------------------ Define some useful functions -------------------------------------------------------------
# creates a 2-column table with children->parent linkages. Takes original gff annotation as its argument.
linkage <- function(gff) { 
  output <- apply(gff, 1, function(x) {
     type   <- x[3]
     id     <- str_match(x[9], 'ID=([^;]+)')[[2]]
     parent <- str_match(x[9], 'Parent=([^;,]+)')[[2]]
     return(c(type, id, parent))
  })
  
  output <- t(output)
  colnames(output) <- c('Type', 'ID', 'Parent')
  output[is.na(output[, 'Parent']), 'Parent'] <- 'Primary'
  
  # sometimes feature have no ID. In that case generate a unique ID as 'generatedID' + [line number]
  extendedID <- sapply(1:nrow(output), function(x) {
            if(is.na(output[x, 'ID'])) { new_id <- paste0('generatedID_',x); return(new_id)}
            else {return(output[x, 'ID'])}
  })
  
  output <- data.table('type' = output[, 'Type'], 'ID' = extendedID, 'Parent' = output[, 'Parent'], 'status' = rep('keep', nrow(output)), stringsAsFactors = FALSE)
  return(output)
}

removeFeatures <- function(gff, parentsTable, featureType){
  #find a top parent (usually a gene name) of every feature you want to remove
  topParents <- c() 
  features <- parentsTable[featureType, on = 'type', ID]
  while(TRUE) {
    candidates <- parentsTable[features, on = 'ID', Parent]
    topParents <- union(topParents, features[candidates == 'Primary'])
    candidates <- candidates[candidates != 'Primary' & !is.na(candidates)]
    if(length(candidates) == 0) { break }
    features   <- parentsTable[candidates, on = 'ID', ID]
  }
  
  # remove all children of the corresponding top parents
  parentsTable[topParents, on = 'ID', status := 'remove']
  children  <- parentsTable[topParents, on = 'Parent', ID]
  while(length(children) > 0) {
    parentsTable[children, on = 'ID', status := 'remove']
    children <- parentsTable[parentsTable$Parent %in% children, ID]
  }
  
  annotation <- gff[parentsTable[, status] %in% 'keep', ]
  return(annotation)
}
#-------------------------------------------------------------------------------------------------------------------------------------
# Load GFF annotation file
gff  <- fread(file="../Original/c_elegans.PRJNA13758.WS268.annotations.gff3", skip = 8, stringsAsFactors = F, header = F, fill = T, na.strings = c("", "NA"), sep="\t") %>% na.omit() #deals with unwanted #comment lines in the gff
gff  <- gff[grepl('ID=|Parent=', gff$V9), ]  # discard feature with no ID in the attributes field. These are typically things like SNPs, TF-binding sites and other genomic features.
gff  <- gff[gff$V2 == 'WormBase', ]          # discard predictions and non-curated junk

# con <- file("../Original/c_elegans.PRJNA13758.WS268.annotations.gff3", "r")
# header <- readLines(con, n = 8)
# write.table(header, file = "Wormbase.gff", col.names = F, row.names = F, quote = F)
# write.table(gff, file = "Wormbase.gff", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
# close(con); rm(con)

parentsTable <- linkage(gff)
setindex(parentsTable, 'ID')
setindex(parentsTable, 'Parent')
setindex(parentsTable, 'type')

# Remove non-coding features
gff2 <- removeFeatures(gff, featureType = c( 'antisense_RNA','nc_primary_transcript','snRNA','lincRNA','ncRNA','tRNA','pre_miRNA','miRNA','scRNA','snoRNA', 'pseudogenic_tRNA', 'piRNA', 'pseudogenic_transcript', 'pseudogenic_rRNA','rRNA'), parentsTable = parentsTable)

# con <- file("../Original/c_elegans.PRJNA13758.WS268.annotations.gff3", "r")
# header <- readLines(con, n = 8)
# write.table(header, file = "WS268_Wormbase_coding.gff3", col.names = F, row.names = F, quote = F)
# write.table(gff2, file = "WS268_Wormbase_coding.gff3", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
# close(con); rm(con)
```
</details>

<details><summary><b>Convert annotation from GFF3 to GTF format</b></summary>  
     
```bash
gffread WS279_Wormbase_coding.gff3 -T -o WS279_Wormbase_coding.gtf
# -T          - convert gff/gtf
```
</details>

```bash  
# rDNA indexing for Bowtie
bowtie-build Elegans_rRNA.fa ./Elegans_indices/Elegans_rRNA  
# Genome indexing for STAR
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir ./Elegans_index/ --genomeFastaFiles ./c_elegans.PRJNA13758.WS268.genomic.fa --sjdbGTFfile ./WS268_Wormbase_coding.gtf
```
</details>



### mRNA-seq sequencing reads filtering and mapping   
<details><summary><b>Illumina adapters trimming</b></summary>

```bash
# trim illumina adapters
cutadapt -j 20 -m 75 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -o trimmed_1.fq.gz -p trimmed_2.fq.gz read.1.fq.gz read.2.fq.gz
# -j      - number of threads
# -m      - discard read pair if any of the mates if shorter than 50 nucleotides after adapter trimming

#automate with bash for loop if needed. In the folder containing sample subfolders (R1 and R2 files in each subfolder) run this:
for dir in */; do
file1=$(find "$dir" -name "*_1.fq.gz");
file2=$(find "$dir" -name "*_2.fq.gz");
echo "$(cutadapt -j 20 -m 75 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -o "$dir""trimmed_1.fq.gz" -p "$dir""trimmed_2.fq.gz" "$file1" "$file2")";
done
```
</details>

