################################################
################################################
########Step1 Transcriptome assembly############

## clutch1 (Tv1) all 5 developmental stages, G43 samples of clutch6, and G43 sample of clutch3 ###
#(this includes one female, one male from each stage, plus 2 individuals per type (XX females, XY0 females, XX males) each for metamorphosis stage, and one XX female and one XY0 male for froglet stage.)#
Trinity --seqType fq --max_memory 300G \
--left Tv1_231_L7_pairedR1.fastq,Tv1_232_L1_pairedR1.fastq,Tv1_271_L1_pairedR1.fastq,Tv1_272_L2_pairedR1.fastq,Tv1_314_L6_pairedR1.fastq,Tv1_317_L3_pairedR1.fastq,Tv1_431_L3_pairedR1.fastq,Tv1_432_L2_pairedR1.fastq,Tv1_462_L7_pairedR1.fastq,Tv1_463_L6_pairedR1.fastq,Tv3_431_L7_pairedR1.fastq,Tv6_431_L8_pairedR1.fastq,Tv6_432_L3_pairedR1.fastq,Tv6_433_L2_pairedR1.fastq \
--right Tv1_231_L7_pairedR2.fastq,Tv1_232_L1_pairedR2.fastq,Tv1_271_L1_pairedR2.fastq,Tv1_272_L2_pairedR2.fastq,Tv1_314_L6_pairedR2.fastq,Tv1_317_L3_pairedR2.fastq,Tv1_431_L3_pairedR2.fastq,Tv1_432_L2_pairedR2.fastq,Tv1_462_L7_pairedR2.fastq,Tv1_463_L6_pairedR2.fastq,Tv3_431_L7_pairedR2.fastq,Tv6_431_L8_pairedR2.fastq,Tv6_432_L3_pairedR2.fastq,Tv6_433_L2_pairedR2.fastq \
--SS_lib_type RF --normalize_reads --CPU 30 --min_kmer_cov 2 --output Trinity_out_Rttv 2>&1 | tee run.log




################################################
################################################
########Step2 Evaluating transcriptome############

######## 2.1) global access quality###############
grep '>' Tv_transc.fasta | wc -l
******558745******

########2.2) examine transcriptome stats#########
/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/TrinityStats.pl Tv_transc.fasta 
****************************************
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  272330
Total trinity transcripts:      558745




################################################
################################################
########Step3 Filtering transcriptome############

############3.1) remove transcripts shorter than 300bp####
perl remove_small.pl 300 Tv_transc.fasta > Tv_transc_300bp.fasta #use custom perl scripts for remove smaller transcripts.
*******349137*********

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  160642
Total trinity transcripts:      349137

################3.2) select highest expressed transcript per gene, and TPM>3 ##############
###quantify abundance of transcripts with Kallisto
module add UHTS/Analysis/kallisto/0.43.0

for f in Am2_23*_pairedR1.fastq
do
kallisto quant -i tv_final_transcripts.idx -o ./output/${f%%_pairedR1.fastq} -b 1000 ${f%%_pairedR1.fastq}_pairedR1.fastq ${f%%_pairedR1.fastq}_pairedR2.fastq
done

### construct an expression matrix 
module add UHTS/Analysis/kallisto/0.43.0
% $TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
    --est_method kallisto  --out_prefix Trinity_tv \
    --name_sample_by_basedir \
    Tv1_231/abundance.tsv \
    Tv1_232/abundance.tsv \
    ......

###based on expression level, filtering transcriptome by selecting the highest expressed transcript per gene.
module add UHTS/Analysis/kallisto/0.43.0
%  $TRINITY_HOME/util/filter_low_expr_transcripts.pl --matrix Trinity_tmm_tv_norm.txt --transcripts Tv_transc_300bp.fasta \
--min_expr_any 1 \
--highest_iso_only 

### the output file
tv_tmm_mostexp.fasta 
Total transcripts
********159109******

################3.3) collapsed haplotypes with 90% identity ##############
module add UHTS/Analysis/cd-hit/4.6.1
cd-hit-est -i tv_tpm3_mostexp.fasta -out tv_tpm3_mostexp_tmm_90iden.fasta -c 0.9 -n 9

**************************************
***tv_tpm3_mostexp_tmm_90iden_cdhit.fasta***67331
**************************************

################3.4) remove ERCC RNA spike-in control transcripts and ribosome transcripts with custom bash scripts##############
### final transcriptome file
tv_mostexp_tmm_90iden.pep
transcript number 67288




################################################
################################################
########Step4 Evaluation the quality of de novo assembled transcriptome###################################

#first need to install BUSCO, then use the following commands.
python3 BUSCO.py -o $RUNNAME.tetrapoda -i $RUNNAME.fasta -l tetrapoda_odb9 -m transcriptome -c 4 -f

#######4.1) checking complete single-copy gene number using BUSCO.

Summarized benchmarks in BUSCO notation:
        C:90%[D:23%],F:1.8%,M:7.9%,n:429
Representing:
        285     Complete Single-copy BUSCOs
        102     Complete Duplicated BUSCOs
        8       Fragmented BUSCOs
        34      Missing BUSCOs
        429     Total BUSCO groups searched

#######4.2) Using bowtie2 to mapping reads of a subset of samples in this population to the Tvedora transcriptome
### Tv reads of 1 clutch with 5 stages)###
351049098 reads; of these:
  351049098 (100.00%) were paired; of these:
    86357889 (24.60%) aligned concordantly 0 times
    260185617 (74.12%) aligned concordantly exactly 1 time
    4505592 (1.28%) aligned concordantly >1 times
    ----
    86357889 pairs aligned concordantly 0 times; of these:
      19008865 (22.01%) aligned discordantly 1 time
    ----
    67349024 pairs aligned 0 times concordantly or discordantly; of these:
      134698048 mates make up the pairs; of these:
        99770902 (74.07%) aligned 0 times
        26082375 (19.36%) aligned exactly 1 time
        8844771 (6.57%) aligned >1 times
85.79% overall alignment rate
