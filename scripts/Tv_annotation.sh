#####################################################
#####################################################
####Functional annotation using Trinotate###
### Trinotate V2.0 ###
**note that Trinotate V3 has a automatic script to blast and laoding the files to databses***
##for toubleshooting, visit this website:http://informatics.fas.harvard.edu/trinotate-workflow-example-on-odyssey.html
###1) Software requirements:
Trinotate
NCBI-Blast
Hmmer
hmmscan
tmhmm
rnammer

###all softwares are installed at Vital-it, hence I use it to run the jobs, otherwise download them from websites, see instructions below.
###2) databases requirements:
#Uniprot
wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/uniprot_sprot.trinotate_v2.0.pep.gz
#Pfam
wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/Pfam-A.hmm.gz
#Trinotate SQLite database
wget "https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/Trinotate.sprot_uniref90.20150131.boilerplate.sqlite.gz" -O Trinotate.sqlite.gz

###3) de novo assembled transcriptome: 
tv_tpm3_mostexp_tmm_90iden.pep

###4) Run the blast, Pfam, rnammer, tmhmm

module add Blast/ncbi-blast/2.6.0+
blastp -query tv_tpm3_mostexp_tmm_90iden.pep -db uniprot_sprot.pep -num_threads 10 -max_target_seqs 1 -outfmt 6 > blastp_tv_tpm3_mostexp_tmm_90iden_pep.txt

module add Blast/ncbi-blast/2.6.0+
blastx -query tv_tpm3_mostexp_tmm_90iden.fasta -db uniprot_sprot.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 > blastx_tv_transcriptome.txt

module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2
hmmscan --cpu 10 --domtblout tv_tpm3_mostexp_tmm_90iden_pep.out Pfam-A.hmm tv_tpm3_mostexp_tmm_90iden.pep > pfam_tv_tpm3_mostexp_tmm_90iden_pep.log

module add SequenceAnalysis/StructurePrediction/signalp/4.1
signalp -f short -n tv_tpm3_mostexp_tmm_90iden_pep.out tv_tpm3_mostexp_tmm_90iden.pep > tv_tpm3_mostexp_tmm_90iden_signalp.out

module add UHTS/Analysis/trinotate/2.0.1
module add UHTS/Analysis/rnammer/1.2
perl /software/UHTS/Analysis/trinotate/2.0.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome tv_tpm3_mostexp_tmm_90iden.fasta --path_to_rnammer /software/UHTS/Analysis/rnammer/1.2/rnammer

module add SequenceAnalysis/StructurePrediction/tmhmm/2.0
tmhmm --short < tv_tpm3_mostexp_tmm_90iden.pep > tv_tpm3_mostexp_tmm_90iden_tmhmm.out

module add UHTS/Analysis/trinotate/2.0.1
perl /software/UHTS/Assembler/trinityrnaseq/2.1.1/util/support_scripts/get_Trinity_gene_to_trans_map.pl tv_tpm3_mostexp_tmm_90iden.fasta > tv_tpm3_mostexp_tmm_90iden_gene_trans_map.txt

###5) Combine search results and create annotation report
module add UHTS/Analysis/trinotate/2.0.1
Trinotate Trinotate.sqlite init --gene_trans_map tv_tpm3_mostexp_tmm_90iden_gene_trans_map.fa --transcript_fasta tv_tpm3_mostexp_tmm_90iden.fasta --transdecoder_pep tv_tpm3_mostexp_tmm_90iden.pep

Load Blast results:
TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_swissprot_blastp  blastp_tv_tpm3_mostexp_tmm_90iden_pep.txt
TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_swissprot_blastx  blastx_tv_transcriptome.txt

Load Pfam results:
RINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_pfam  pfam_tv_tpm3_mostexp_tmm_90iden_pep.log

Load results of additional processing (tmHMM, SignalP, RNAMMER):
TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_tmhmm  tv_tpm3_mostexp_tmm_90iden_tmhmm.out
TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_signalp  tv_tpm3_mostexp_tmm_90iden_signalp.out
TRINOTATE_HOME/Trinotate Trinotate.sqlite LOAD_rnammer  tv_tpm3_mostexp_tmm_90iden.rnammer.gff

###Each of these database loading commands takes about a minute (or less) to run.
###Create annotation report table; a parameter that can be specified is "-E" (maximum E-value for reporting best blast hit and associated annotations); the default value (0.00001) is stringent enough, so that is what is used here:
module add UHTS/Analysis/trinotate/2.0.1
TRINOTATE_HOME/Trinotate Trinotate.sqlite report -E 0.00001 > tv_tpm3_mostexp_tmm_90iden_annotation_report.xls


###6)For annotation of the genome, we use custom python scripts (see blast_RBH.py) to perform reciprocal best blast hits to search for one-to-one orthologs between Rana temporaria and Xenopus tropicalis.
module add Blast/ncbi-blast/2.6.0+
python blast_RBH.py -a prot -t blastp -i 30 -c 50 -o tv_mostexp_tpm3_tmm_90ident.tsv tv_tpm3_mostexp_tmm_90iden.pep Xtropicalisv9.0.Named.primaryTrs.pep.fa

### Finally, combining the Trinotate annotation and Xenopus orthologs, using custom bash scripts, I make a clean annotation file for downstream analysis.
tv_annotation.txt

###7)Using custom bash and python scripts to generate gene ontology annotation of the transcriptome file
GO_tv.txt