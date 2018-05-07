#loading module from university cluster, Trinity version: v2.4.0
module add UHTS/Assembler/trinityrnaseq/2.4.0

#this transcriptome assembly is based on five development stages of samples from clutch1, G43 samples of clutch6, and G43 sample of clutch3.
#this includes one female, one male from each stage, plus 2 individuals per type (XX females, XY0 females, XX males) each for metamorphosis stage, and one XX female and one XY0 male for froglet stage.
#this is to include samples from all stages, meanwhile we can keep the haplotypes as low number as possible. 

Trinity --seqType fq --max_memory 300G \
--left Tv1_231_L7_pairedR1.fastq,Tv1_232_L1_pairedR1.fastq,Tv1_271_L1_pairedR1.fastq,Tv1_272_L2_pairedR1.fastq,Tv1_314_L6_pairedR1.fastq,Tv1_317_L3_pairedR1.fastq,Tv1_431_L3_pairedR1.fastq,Tv1_432_L2_pairedR1.fastq,Tv1_462_L7_pairedR1.fastq,Tv1_463_L6_pairedR1.fastq,Tv3_431_L7_pairedR1.fastq,Tv6_431_L8_pairedR1.fastq,Tv6_432_L3_pairedR1.fastq,Tv6_433_L2_pairedR1.fastq \
--right Tv1_231_L7_pairedR2.fastq,Tv1_232_L1_pairedR2.fastq,Tv1_271_L1_pairedR2.fastq,Tv1_272_L2_pairedR2.fastq,Tv1_314_L6_pairedR2.fastq,Tv1_317_L3_pairedR2.fastq,Tv1_431_L3_pairedR2.fastq,Tv1_432_L2_pairedR2.fastq,Tv1_462_L7_pairedR2.fastq,Tv1_463_L6_pairedR2.fastq,Tv3_431_L7_pairedR2.fastq,Tv6_431_L8_pairedR2.fastq,Tv6_432_L3_pairedR2.fastq,Tv6_433_L2_pairedR2.fastq \
--SS_lib_type RF --normalize_reads --CPU 30 --min_kmer_cov 2 --output Trinity_out_Rttv 2>&1 | tee run.log