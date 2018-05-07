###quantify abundance of transcripts with Kallisto

module add UHTS/Analysis/kallisto/0.43.0

for f in *_pairedR1.fastq
do
kallisto quant -i tv_final_transcripts.idx -o ./output/${f%%_pairedR1.fastq} -b 1000 ${f%%_pairedR1.fastq}_pairedR1.fastq ${f%%_pairedR1.fastq}_pairedR2.fastq
done
