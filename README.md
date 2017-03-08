# Functional_Genomics_GP1
#Code used for Auburn University Functional Genomics Spring 2017 class 
#Contributors 

Data:Dr.Tonia Schwatz, Auburn University

Code: Chase Rushton , Auburn University , Stephen Gardner , Auburn University 

Resources: Alabama Super Computer, Huntsville Alabama 

Project:The epigenetic signature of developmentally induced desiccation resistance in a lizard.

Tonia Schwartz, Corey Cates, Nicole Riddle, David Allison, Dan Warner



###Goals
Analyse Methyl DIPSEQ data from Anolis measung thermal resistance differences. 

####Samples:   
Female cuban brown anoles, Anolis sagrei, were collected from their field sites in Florida (islands in the Tomoka River in Florida, an intertidal estuary with high salinity fluctuations) and brought them back to the lab, along with soil samples from each of the islands. The “arid” islands have more shell-based soil compare to the soil on the “wet” islands. The lizards were randomly assigned to cages and mated. Throughout the spring/summer these lizards produce one egg every few days and the eggs were continuously rolled into the experiment. Eggs from a particular female were distributed across four treatments: Dry soil (-600 kpa), Dry shell (-600kpa), Wet soil (-30 kpa), Wet shell (-30 kpa).

We used ten male hatchlings from the wet and dry soil treatments (5 from each treatment: 12-19 days old) that were frozen after their test for desiccation resistance.  The skin on ventral side of torso was removed and minced.  The ventral skin was chosen as this is the location where (these) lizards lose water through evaporation (ref).  DNA was isolated using DNeasy (Qiagen) columns including a DTT treatment during digestion to improve the digestion of the skin tissue, and a RNAse treatment prior to loading on the spin column. 
DNA from 10 individuals was sent for MeDIP-Seq and Illumina library preparation through EpiGenTek.  MeDIP-Seq uses an antibody to capture DNA fragments containing 5-mC. This technique uses chromatin immunoprecipitation to enrich for regions of the genome that are highly methylated. DNA from two individuals were also used as input controls (Individual 1 and 10) where no antibody was added. The libraries (12 total) were barcoded and sequenced on the Illumina Hi-Seq 2500, for 100 bp Single end sequencing.

#Genomic References
The Anolis carolinensis genome and annotation
Assembly
ftp://ftp.ensembl.org/pub/release-87/fasta/anolis_carolinensis/dna/

Annotation
ftp://ftp.ensembl.org/pub/release-87/gff3/anolis_carolinensis/
ftp://ftp.ensembl.org/pub/release-87/gtf/anolis_carolinensis/


###Programs used

Picard 1.79

GATK 3.4-46

gnu_parallel/201612222

fastqc/0.10.1

trimmomatic/0.35

bwa/0.7.12

samtools/1.2

###Step1- Run FastQC to check data quality

```Shell
#Get list of files to run
ls | grep ".fastq$" > list.txt
#Define loop
for x in {1..36}  #36 comes form wc -l of list.txt

do

echo "#! /bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh">$x.fastqc.sh

    echo "module load fastqc/0.10.1">>$x.fastqc.sh
    input=$(awk "NR==$x {print \$0}" ./list.txt)
    echo "cd /scratch/aubcls07/trim1" >>$x.fastqc.sh
    echo "fastqc ${input}">>$x.fastqc.sh
    chmod a+x $x.fastqc.sh
done
```

>To run
>IF you don't want to manually input in the queue parameters each time edit your .asc_queue fil

```Shell
emacs .asc_queue
```

>Run this file on the command line to generate each input.fastqc.sh file then  run the following on the command line to submit all 36 at once

```Shell
for x in {1..36}
do
run_script $x.fastqc.sh
done
```
###Step2 - Trimmomatic 
```Shell
for i in files 
do
############ Trimmomatic #############

java -Xms2g -Xmx16g -jar /opt/asn/apps/trimmomatic_0.35/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 6 -phred33 "$i"_All_R1.fastq.gz "$i"_All_R1_trimmed.fastq.gz ILLUMINACLIP:AdaptersToTrim.fa:2:30:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36 
java -Xms2g -Xmx16g -jar /opt/asn/apps/trimmomatic_0.35/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 6 -phred33 "$i"_CTTGTA_All_R1.fastq.gz "$i"_CTTGTA_All_R1_trimmed.fastq.gz ILLUMINACLIP:AdaptersToTrim.fa:2:30:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36

done

```

Step 3 - Index with with BWA

```
bwa index -p  Anolis_carolinensis_genome -a bwtsw  Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa
```
Step 4- Mapping with BWA
```Shell
or x in {1..12} #12 comes from wc -l of list.txt
do
        echo "module load bwa/0.7.12">$x.alignment.sh
        echo "module load samtools/1.2">>$x.alignment.sh
    input=$(awk "NR==$x {print \$0}" ./list.txt)
    echo "cd /scratch/aubcls07/group/trim/" >>$x.alignment.sh
        echo "bwa mem -M -v 2 -t 12 Anolis_carolinensis_genome ${input} > ${input}.bam" >>$x.alignment.sh
        echo "samtools view -Sbu ${input}.bam | samtools sort - -@ 12 ${input}.sort" >>$x.alignment.sh
    chmod a+x $x.alignment.sh
done
```
Step 5 - Build Bam index with picard 

```Shell
java -Xms2g -Xmx14g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/BuildBamIndex.jar I= Anolis_carolinensis_genome.fa O= Anolis_carolinensis_genome.fa.fai
```

Step 6 - Build dict with picard 

```Shell
java -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/CreateSequenceDictionary.jar  R= Anolis_carolinensis_genome.fa  O= Anolis_carolinensis_genome.dict
```

Step 7 - Marking Duplicates with Picard tools 
```Shell
or x in {1..12} #12 comes from wc -l of list.txt
do
        echo "module load picard">$x.dup.sh

    input=$(awk "NR==$x {print \$0}" ./list.txt)
    echo "cd /scratch/aubcls07/group/trim/sortedbam" >>$x.dup.sh
    echo "java -Xms2g -Xmx14g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/MarkDuplicates.jar INPUT=${input} OUTPUT= ${input}.mkdup.bam METRICS_FILE=Anol.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500" >>$x.dup.sh

    chmod a+x $x.dup.sh
done
```


Step 8 - Run stats on the output with Samtools

```Shell
ls | grep ".bam$" > list.txt
files=`cat list.txt`
ref=Anolis_carolinensis_genome.fa

for i in ${files}
do
samtools index ${i}
samtools flagstat ${wd}${i} > ${wd}${i}.flagstat.txt
samtools stats ${wd}${i} > ${wd}${i}.stats.txt
samtools depth ${wd}${i} > ${wd}${i}.depth.txt
samtools idxstats ${wd}${i}  > ${input}_Counts.txt
done
```

Step 9 - Average depth 
```Shell
ls | grep "depth.txt$" > list.txt
files=`cat list.txt`
for i in ${files}
do
cat ${i} awk | '{sum+=$3} END {print sum/NR}' > ${i}.avg.depth.txt
done
```
Step 10 - Graphics in R on statistics 

```R
```

