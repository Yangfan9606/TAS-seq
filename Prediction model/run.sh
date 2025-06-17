#!/bin/bash

#--- Random select X number for training (100000 POS + 100000 NEG required ~3h)
shuf -n100000 Positive_A_sites_from_HEK293T_NES|awk -F "\t" '{i=$3-10;j=$3+10;print "Chr"$2"\t"i"\t"j"\t"$1"\t1\t"$4}' > POS.bed
shuf -n100000 Negative_A_sites_from_HEK293T_NES|awk -F "\t" '{i=$3-10;j=$3+10;print "Chr"$2"\t"i"\t"j"\t"$1"\t0\t"$4}' > NEG.bed
#--- Downlad the reference fasta and rename Chr
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
perl fna.NCBI.Chr_Rename.pl GCF_000001405.40_GRCh38.p14_genomic.fna.gz -o GRCh38.p14.fa
#--- Find the 21bp sequence
perl GetSeqFromFa.pl GRCh38.p14.fa -f POS.bed
perl GetSeqFromFa.pl GRCh38.p14.fa -f NEG.bed
cat POS.bed.fa NEG.bed.fa|sed ":a;N;s/\n/\t/;ta"|sed "s/>/\n/g"|sed "s/:/\t/g"|awk -F "\t" '$1!=""{print $8" "$5}' > Train_input.txt
#--- Train model: Output is 'rna_classifier_model.pth'
python Sequence_classifier.py Train_input.txt
#--- Predict
python Sequence_predictor.py -m rna_classifier_model.pth -s ATCAATATCTATAACTGCATT
python Sequence_predictor.py -m rna_classifier_model.pth -i input.sequence.txt -o output.txt
