###assessing raw reads quality
module add UHTS/Quality_control/fastqc/0.11.2  ##load the module 

fastqc  read_R1.fastq.gz read_R2.fastq.gz


#load module from University cluster.
module add UHTS/Analysis/trimmomatic/0.33 

#provide sequences of adapters, indexes, flow cell attached DNA fragments, and each of their reverse complementary sequences.
ADAPTERS="~/software/Trimmomatic-0.33/adapters/TruSeq2-PE.fa" 

#for loop to trim raw reads for each RNA sample using Trimmomatic v0.33 with default parameters.

for f in *_R1.fastq.gz

do
     trimmomatic PE -phred33 \
     ${f%%_R1.fastq.gz}_R1.fastq.gz \
     ${f%%_R1.fastq.gz}_R2.fastq.gz \
     ${f%%_R1.fastq.gz}_pairedR1.fastq.gz \
     ${f%%_R1.fastq.gz}_unpairedR1.fastq.gz \
     ${f%%_R1.fastq.gz}_pairedR2.fastq.gz \
     ${f%%_R1.fastq.gz}_unpairedR2.fastq.gz \
     ILLUMINACLIP:$ADAPTERS:2:30:10 \
     LEADING:3 \
     TRAILING:3 \
     SLIDINGWINDOW:4:15 \
     MINLEN:36 &> ${f%%R1combined.fastq.gz}.trim.log

done

###assessing trimmed reads quality
module add UHTS/Quality_control/fastqc/0.11.2

fastqc  read_trim_pairedR1.fastq.gz read_trim_pairedR2.fastq.gz