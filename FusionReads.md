# These codes were used to identify vector-genome fusion reads, through three rounds of mapping
There are 7 types of reads: 
|Yes|Yes|Yes|Yes|Yes|Yes|Yes|
|---|---|---|---|---|---|---|
|Junk|Vector|Genome|Vector-junk|Genome-junk|Genome-vector|Vector-vector|

*Assuming that junk may not be able to map to either genome and vector*

## Round 1: Get reads that contain vector sequences
### Extract the 20bp at 5' end
*Base pair we extract can be changed
```
#updated trimfastq.py for python3
python trimfastq.py alltrimmedfastq 20 > 20fiveprime.fastq
```
### Map to vector and extract mapped reads
```
bowtie2 Original_vector -p 8 -p 8 --very-sensitive -k 1 --time -q --al alltrimmedVectorfastq -U 20fiveprime.fastq -S 20fiveprime.sam
```
### Retrieve the original full length of the mapped reads
```
grep -f <(cat alltrimmedVectorfastq | paste - - - - | cut -f 1) \
<(cat alltrimmedfastq | sed 's/ /_/g' | paste - - - - ) | tr "\t" "\n" \
> alltrimmedVectorFullLengthfastq
```
Now we have:
|No|Yes|No|Yes|No|Yes|Yes|
|---|---|---|---|---|---|---|
|~Junk~|Vector|~Genome~|Vector-junk|~Genome-junk~|Genome-vector|Vector-vector|

## Round 2: Get reads that aren't just about vector sequences and vector-vector sequences
### Extract the 20bp at 3' end
```
paste <(cat alltrimmedVectorFullLengthfastq | paste - - - - | cut -f1 | tr "\t" "\n") \
    <(cat alltrimmedVectorFullLengthfastq | paste - - - - | cut -f2 | tr "\t" "\n" | grep -o '.\{20\}$' )\
    <(cat alltrimmedVectorFullLengthfastq | paste - - - - | cut -f3 | tr "\t" "\n") \ 
    <(cat alltrimmedVectorFullLengthfastq | paste - - - - | cut -f4 | tr "\t" "\n" | grep -o '.\{20\}$' ) | tr "\t" "\n" > 20threeprime.fastq
```
### Map to vector index and get unmapped reads
```
bowtie2 Original_vector -p 8 --very-sensitive -k 1 --time -q --un 3primeVectorUnmapfastq -U 20threeprime.fastq -S /dev/null
```
We take the unmapped reads.
Now we have:
|No|No|No|Yes|No|Yes|
|---|---|---|---|---|---|
|~~Junk~~|~~Vector~~|~~Genome~~|Vector-junk|~~Genome-junk~~|Genome-vector|

## Round 3: Get reads that contain genome sequences
### Extract the 20bp at 3' end
```
paste <(cat alltrimmedVectorUnmapfastq | paste - - | cut -f1 | tr "\t" "\n") \
<(cat alltrimmedVectorUnmapfastq | paste - - | cut -f2 | tr "\t" "\n" | grep -o '.\{'20'\}$' ) | tr "\t" "\n" > 20threeprime.fastq
 ```
### Map to dm3 genome
```
bowtie genome/bowtie-indexes/dm3 -p 8 -v 0 -k 1 -m 1 -t --sam-nh --best --strata -y -q --sam 20threeprime.fastq | \
samtools-0.1.8/samtools view -bT genome/bowtie-indexes/dm3.fa - | \
samtools-0.1.16/bin/samtools sort - 20threeprime.dm3.bam
samtools-0.1.16/bin/samtools index 20threeprime.dm3.bam.bam
```
Now we have:
|No|No|No|No|No|Yes|
|---|---|---|---|---|---|
|~~Junk~~|~~Vector~~|~~Genome~~|~~Vector-junk~~|~~Genome-junk~~|Genome-vector|

Then the bam files can be loaded onto genome browser to see whether the fusion reads land 
