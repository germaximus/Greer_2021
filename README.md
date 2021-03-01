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

<details><summary><b>Read mapping and counting with STAR</b></summary>
     
```bash
STAR --genomeLoad LoadAndExit --genomeDir ../STAR-2.7.2b/Mouse_male_index/ 	# load genome once in the shared memory
STAR --runThreadN 40 --outSAMtype BAM Unsorted --outSAMmultNmax 1 --quantMode GeneCounts --genomeLoad LoadAndKeep --genomeDir ../STAR-2.7.2b/Mouse_male_index/ --readFilesCommand gunzip -c --readFilesIn trimmed_1.fastq.gz trimmed_2.fastq.gz --outFileNamePrefix ./OUT_folder/ 
STAR --genomeLoad Remove 	# remove loaded genome from shared memory
# ipcs - check shared memory consumption
# ipcrm - remove object from shared memory

#automate with bash for loop if needed. 
#females
for dir in *[1-4]/; do
echo "$(STAR --runThreadN 40 --outSAMtype BAM Unsorted --outSAMmultNmax 1 --quantMode GeneCounts --genomeLoad LoadAndKeep --genomeDir ../../STAR-2.7.2b/Mouse_female_index/ --readFilesCommand gunzip -c --readFilesIn "$dir""trimmed_1.fastq.gz" "$dir""trimmed_2.fastq.gz" --outFileNamePrefix ./"$dir")";
done
#males
for dir in *[6-9]/; do
echo "$(STAR --runThreadN 40 --outSAMtype BAM Unsorted --outSAMmultNmax 1 --quantMode GeneCounts --genomeLoad LoadAndKeep --genomeDir ../../STAR-2.7.2b/Mouse_male_index/ --readFilesCommand gunzip -c --readFilesIn "$dir""trimmed_1.fastq.gz" "$dir""trimmed_2.fastq.gz" --outFileNamePrefix ./"$dir")";
done
#rename generic file names
for dir in */; do
basename="${dir%/}";
mv ./"$dir"ReadsPerGene.out.tab ./"$dir"/"$basename".counts
done
#copy files in the same directory for easier transfer
for dir in */; do
basename="${dir%/}";
cp ./"$dir"/"$basename".counts ./;
done
#transfer all files to a local windows machine (geneCounts folder)
pscp -pw password germax@aging:/PATH/diets/mrna/*.counts .
```
</details>
