# Functional_Genomics_GP1
#Code used for Auburn University Functional Genomics Spring 2017 class 
#Contributors 
Data:Dr.Tonia Schwatz, Auburn University
Code: Chase Rushton , Auburn University 
Resources: Alabama Super Computer, Huntsville Alabama 

 

###Goals
Analyse Methyl DIPSEQ data from Anolis measung thermal resistance differences. 

####Samples:   
10 Anolis RNA seq data sets provided by Dr.Schwartz
reference= Anolis_carolinensis_genome.fa

###Programs used
Picard
GATK
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
or x in {1..12} #125 comes from wc -l of Table.txt
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

Step 6 - Marking Duplicates with Picard tools 
```Shell
or x in {1..12} #125 comes from wc -l of Table.txt
do
        echo "module load picard">$x.dup.sh

    input=$(awk "NR==$x {print \$0}" ./list.txt)
    echo "cd /scratch/aubcls07/group/trim/sortedbam" >>$x.dup.sh
    echo "java -Xms2g -Xmx14g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/MarkDuplicates.jar INPUT=${input} OUTPUT= ${input}.mkdup.bam METRICS_FILE=Anol.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500" >>$x.dup.sh

    chmod a+x $x.dup.sh
done
```

Step 7 - Run stats on the output with Samtools

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


