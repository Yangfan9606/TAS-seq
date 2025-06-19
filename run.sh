#!/bin/bash
#------------------------------
Reference_fasta="GRCh38.p14.fa"
Perl_bin="./Perl_scripts"
# SNP: 
weget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz
zcat common_all_20180418.vcf.gz | awk -F "\t" '$1!~/^#/ {print $1"\t"$2"\t"$2"\t"$3"\t"$4"-"$5"\t+"}'|gzip > human.SNP.bed.gz
SNP_file="human.SNP.bed.gz"
#--------------------------

# Preprocessing:
sample="Ctl.Rep1 Ctl.Rep2 Treat.Rep1 Treat.Rep2"
for s in $sample
do
samtools mpileup -q 25 -Q 20 -f $Reference_fasta $s.uniq.bam | gzip > $s.uniq.bam.Q25bq20.depth.gz
perl $Perl_bin/Bam.RNA_Stranded.AT_sites.depth.pl $s.uniq.bam.Q25bq20.depth.gz -o $s.uniq.Q25bq20d1.AT
perl $Perl_bin/Base.site_type.depth.pl $s.uniq.Q25bq20d1.AT.gz -o $s.uniq.Q25bq20d1.AT.sites -chr 6
perl $Perl_bin/Base.mutation.d1.pl $s.uniq.Q25bq20d1.AT.sites -o $s.uniq.Q25bq20d1.AT.sites.mutation
read_num=$(samtools flagstat $s.uniq.bam|grep "properly paired"|awk 'print $1')
perl $Perl_bin/Mutation_score_cal.pl $s.uniq.Q25bq20d1.AT.sites.mutation $read_num -o $s.uniq.Q25bq20d1.AT.sites.mutation.score
perl $Perl_bin/Snp.filter.pl $SNP_file $s.uniq.Q25bq20d1.AT.sites.mutation.score $s.uniq.Q25bq20d1.AT.sites.mutation.score.noSNP
mu1=uniq.Q25bq20d1.AT.sites.mutation.score.noSNP.A2G_mu1
cut -f1-5 $s.uniq.Q25bq20d1.AT.sites.mutation.score.noSNP|awk 'NR!=1' > $s.$mu1.sites
perl $Perl_bin/A2G.RNA_straned.mutation_summary.pl $s.$mu1.sites -o $s.Filter.A2G.sites.bed
done
# Comparison:
F=Filter.A2G.sites.bed
Ctl_R1_depth=(samtools flagstat Ctl.Rep1.uniq.bam|grep "properly paired"|awk 'print $1')
Ctl_R2_depth=(samtools flagstat Ctl.Rep2.uniq.bam|grep "properly paired"|awk 'print $1')
Treat_R1_depth=(samtools flagstat Treat.Rep1.uniq.bam|grep "properly paired"|awk 'print $1')
Treat_R2_depth=(samtools flagstat Treat.Rep2.uniq.bam|grep "properly paired"|awk 'print $1')
perl $Perl_bin/Delta_rate.cal.pl Treat.Rep1.$F Ctl.Rep1.$F -treat $Treat_R1_depth -control $Ctl_R1_depth -o Compare.Rep1.delta_rate
perl $Perl_bin/Delta_rate.cal.pl Treat.Rep1.$F Ctl.Rep1.$F -treat $Treat_R1_depth -control $Ctl_R1_depth -o Compare.Rep1.delta_rate
# Statistic test:
perl $Perl_bin/Fishers.exact.test.pl Compare.Rep1.delta_rate -h 1 -s1 7 -t1 8 -s2 10 -t2 11 -name Treat.R1.A2G_count,Treat.R1.total_count,Ctl.R1.A2G_count,Ctl.R1.total_count -o Compare.Rep1.delta_rate.p
perl $Perl_bin/Fishers.exact.test.pl Compare.Rep2.delta_rate -h 1 -s1 7 -t1 8 -s2 10 -t2 11 -name Treat.R2.A2G_count,Treat.R2.total_count,Ctl.R2.A2G_count,Ctl.R2.total_count -o Compare.Rep2.delta_rate.p
