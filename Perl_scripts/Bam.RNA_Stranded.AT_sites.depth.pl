#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-18
# version:   3.0
# license:   MIT
# brief:     统计bam文件每个A/T的depth。
#------------------------------------------#
use strict;
use Getopt::Long;
use lib './';
use yangfan;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : 统计 每个[A/T]位点的reads (输出RNA strand信息)

Usage: .pl [.bam.Q25bq20.depth.gz] [.uniq.bam] -o [out.name]
        (gz) header =  1
        depth >= 1
        only A/T
        MAPQ 25
        base quility 20
        -o                      out name
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "h:s", "c:s", "f:s", "nf:s", "qf:s", "help!");
##########
die $usage if ( @ARGV!=2 || defined($opts{"help"}));

###################################################
#                  START                          #
###################################################
my $optko;
foreach my $opt(keys %opts){
        $optko .= " -$opt $opts{$opt}";
}
print "##########Start############ perl $0 @ARGV ($optko)\n";
Ptime("Start");
my $infile=shift;
my $bam_file=shift;
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";
open(OUT, "|gzip > $outname.gz") or die $!;
my @low20 = ("b_!","b_\"","b_#","b_\$","b_%","b_&","b_'","b_(","b_)","b_*","b_+","b_,","b_-","b_.","b_/","b_0","b_1","b_2","b_3","b_4","b_5");
my %low;
foreach my $l2(@low20){
        $low{$l2} = 1;
}
#########
my @OB=qw(A C G T);
my $cn = 100;
my @order=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M 45s_NR_046235.3 U1_NR_004430 U2_NR_002716.3 U4_NR_003925.1 U4atac_NR_023343.1 U5_NR_002756.2 U6_NR_004394.1 U6atac_NR_023344.1 U11_NR_004407.1 U12_NR_029422.2);
my (%all,%all1);
for (my $cr=0;$cr<@order;$cr+=$cn){
        my %Flag;
        for (my $i=0;$i<$cn&&($cr+$i)<@order;$i++){
                $Flag{$order[$cr+$i]} = 1;
        }
#############
open(IN,"zcat $infile|") or die $!;
open(BAM,"samtools view $bam_file|sort -k1|")or die $!;
my (%Lake,%R,%ST,%reads);
while(<IN>){
        chomp;
        my @c = split/\t/;
        my $chr = $c[0];
        next unless $Flag{$chr} == 1;
        my $site = $c[1];
        my $base = uc $c[2];
        next unless ($base eq "A" || $base eq "T");
        my $depth = $c[3];
        next unless $depth >= 1;
        $Lake{$chr}{$site} = "$base\t$depth";
}
close IN;
print "read done...\n\tprcessing bam...\n";
while(<BAM>){
        my $L1 = $_;
        my $L2 = <BAM>;
        chomp $L1;
        chomp $L2;
        my @reads = ($L1,$L2);
        my $mF = 0;
        my $chr;
        my ($MO1,$MO2);
        for (my $n=0;$n<2;$n++){
                my $R = $reads[$n];
                my @c = split/\t/,$R;
                my $flag = $c[1];
                my $MAPQ = $c[4];
#################################   MAPQ filter
                next if $MAPQ < 25;
                my $strand = "-";
                my $read = "1";
                if (($flag&32)==32){
                        $strand = "+";
                }
                if (($flag&128)==128){
                        $read = "2";
                }
                $chr = $c[2];
                next unless $Flag{$chr} == 1;
                my $s = $c[3];
                my @seq = split//,$c[9];
                my @qual = split//,$c[10];
                my $CIGAR = $c[5];
                my $len=0;
                if ($CIGAR =~ /^(\d+)([A-Z])/){
                        if ($2 eq "M"){
                                for (my $i=0;$i<$1;$i++){
                                        my $site = $s+$i;
                                        if ($Lake{$chr}{$site} ne ""){
                                                $mF = 1;
                                                my $b = $seq[$i];
                                                my $q = $qual[$i];
                                                if ($low{"b_$q"} == 1){
                                                }else{
                                                        my $ori = "$read,$b,$strand";
                                                        my $R_strand = $strand;
                                                        if ($read eq "1"){
                                                                if ($strand eq "+"){
                                                                        $b =~ tr/ATGC/TACG/;
                                                                        $R_strand = "-";
                                                                }else{
                                                                        $R_strand = "+";
                                                                }
                                                        }else{
                                                                if ($strand eq "-"){
                                                                        $b =~ tr/ATGC/TACG/;
                                                                }
                                                        }
                                                        print OUT "$chr\t$site\t$Lake{$chr}{$site}\t$b\t$R_strand\t$ori\t$c[0]\n";
                                                }
                                        }
                                }
                                $s+=$1;
                                $len+=$1;
                        }
                        if ($2 eq "D"){
                                $s+=$1;
                        }
                        if ($2 eq "I"){
                                $len+=$1;
                        }
                        if ($2 eq "N"){
                                $s+=$1;
                        }
                }
                while($CIGAR =~ /(?=[A-Z](\d+)([A-Z]))/g){
                        if ($2 eq "M"){
                                my $j=0;
                                for (my $i=$len;$i<$len+$1;$i++){
                                        my $site = $s+$j;
                                        $j++;
                                        if ($Lake{$chr}{$site} ne ""){
                                                $mF = 1;
                                                my $b = $seq[$i];
                                                my $q = $qual[$i];
                                                if ($low{"b_$q"} == 1){
                                                }else{
                                                        my $ori = "$read,$b,$strand";
                                                        my $R_strand = $strand;
                                                        if ($read eq "1"){
                                                                if ($strand eq "+"){
                                                                        $b =~ tr/ATGC/TACG/;
                                                                        $R_strand = "-";
                                                                }else{
                                                                        $R_strand = "+";
                                                                }
                                                        }else{
                                                                if ($strand eq "-"){
                                                                        $b =~ tr/ATGC/TACG/;
                                                                }
                                                        }
                                                        print OUT "$chr\t$site\t$Lake{$chr}{$site}\t$b\t$R_strand\t$ori\t$c[0]\n";
                                                }
                                        }
                                }
                                $s+=$1;
                                $len+=$1;
                        }
                        if ($2 eq "D"){
                                $s+=$1;
                        }
                        if ($2 eq "I"){
                                $len+=$1;
                        }
                        if ($2 eq "N"){
                                $s+=$1;
                        }
                }
        }
}
close BAM;
print "bam read done...\n";
}
#############

close OUT;

Ptime("End");
print "##########End############\n";
