# Greer_2021
Transcriptome and ribosome profiling analysis in C. elegans  




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
