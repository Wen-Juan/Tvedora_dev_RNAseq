#before running Prank, need to prepare the input file.
#for codon alignment, it requires the sequences to be length which can multiply by 3.

module add UHTS/Assembler/TransDecoder/2.0.1
TransDecoder.LongOrfs -t file.fasta > file_output.cds
TransDecoder.Predict -t target_transcripts.fasta #the final cds or .pep file are not uniqe, need to remove duplicates with shorter sequence lengths with custom scripts

# note TransDecoder normally generates more than one open reading frames per transcript.
#select the longest ORF per transcript
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  input.fasta  # unwrap fasta sequences
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' #calculate size
sort -t ' ' -k2,2 -k1,1nr  #sort on name, inverse length
sort -k1,1 -u -s) #sort on name, unique, stable sort
sed 's/    /./'  #restore name
cut -f 1,2 #cut name, sequence
tr "\t" "\n" < linearized.fasta #go back to fasta
tr "\t" "\n" < linearized.fasta | fold -w 60 #pretty fasta

# run prank with codon function, with script "for_loop_codeml.sh".

# run PAML with codeml function
