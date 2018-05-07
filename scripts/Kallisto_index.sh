#loading necessary modules and databases
module add UHTS/Analysis/kallisto/0.43.0

kallisto index -i tv_transcripts.idx /scratch/beegfs/monthly/wjma/tv/kallisto/tv_tpm3_mostexp_tmm_90iden.fasta
