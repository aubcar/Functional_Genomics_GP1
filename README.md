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


###Programs used
FastQC
Samtools
Trimmomatic
BWA
Picard
GATK

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
