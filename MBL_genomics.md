### Unix Tips

[00b_Unix_Intro.pdf](https://github.com/LeslieDay/UsefulScripts/files/9095439/00b_Unix_Intro.pdf)

### Genome Assembly 
[01_Genome_Assembly.pdf](https://github.com/LeslieDay/UsefulScripts/files/9095551/01_Genome_Assembly.pdf)

### Login 

```bash
ssh -Y -i ~/.ssh/MyPrivateKey.pem leslie@18.219.176.206

#enter directory for MBL things 
cd MBL/ExampleAssembly

#copy R1 and R2 reads in for assembly 
cp /datahaus/genomes_for_assembly/ACB* .

#count number of sequences in file 
zcat ACB001_GAGAATGGTT-TCGGCAGCAA_L002_R1_001.fastq.gz | wc -l
>5480456

zcat ACB001_GAGAATGGTT-TCGGCAGCAA_L002_R2_001.fastq.gz | wc -l
>5480456

#run fastqc on both files 
mkdir fastqc_out

fastqc ACB001_GAGAATGGTT-TCGGCAGCAA_L002_R1_001.fastq.gz ACB001_GAGAATGGTT-TCGGCAGCAA_L002_R2_001.fastq.gz -o fastqc_out

#use multiqc to aggregate fastqc outputs 
multiqc fastqc_out

#trim reads
trimmomatic PE -baseout YOURFILE.trimmed.fastq.gz \
ACB001_GAGAATGGTT-TCGGCAGCAA_L002_R1_001.fastq.gz ACB001_GAGAATGGTT-TCGGCAGCAA_L002_R2_001.fastq.gz \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/combined.fa:2:30:10 \
SLIDINGWINDOW:4:15 LEADING:2 TRAILING:2 MINLEN:100

#assemble using spades
spades.py -o ACB_spades_assembly -1 YOURFILE.trimmed_1P.fastq.gz \
-2 YOURFILE.trimmed_2P.fastq.gz -t 2 -m 20

quast.py -o quast_output -t 2 /ACB_spades_assembly/scaffolds.fasta
```
