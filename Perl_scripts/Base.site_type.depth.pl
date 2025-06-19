#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-18
# version:   1.0
# license:   MIT
# brief:     将A/T的depth整合。
#------------------------------------------#
use strict;
use Getopt::Long;
use lib './';
use yangfan;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : Processing the file

Usage: .pl [K562.8eRNANLS.uniq.Q25bq20d10.AT.gz] -o [out.name]
        -chr                    chr number [100]
        -o                      out name
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "h:s", "c:s", "f:s", "chr:s", "help!");
##########
die $usage if ( @ARGV!=1 || defined($opts{"help"}));

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
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";

open(OUT, ">$outname") or die $!;

my $cn = $opts{chr} eq ""?100:$opts{chr};
#my @order=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M 45s_NR_046235.3 U1_NR_004430 U2_NR_002716.3 U4_NR_003925.1 U4atac_NR_023343.1 U5_NR_002756.2 U6_NR_004394.1 U6atac_NR_023344.1 U11_NR_004407.1 U12_NR_029422.2);
my @order=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M);
###45s_NR_046235.3 U1_NR_004430 U2_NR_002716.3 U4_NR_003925.1 U4atac_NR_023343.1 U5_NR_002756.2 U6_NR_004394.1 U6atac_NR_023344.1 U11_NR_004407.1 U12_NR_029422.2);
#############
for (my $cr=0;$cr<@order;$cr+=$cn){
        my %Flag;
        for (my $i=0;$i<$cn&&($cr+$i)<@order;$i++){
                $Flag{$order[$cr+$i]} = 1;
        }
my (%Lake,%reads);
open(IN,"zcat $infile|") or die $!;
while(<IN>){
        chomp;
        my @c = split/\t/;
        my $chr = $c[0];
        next unless $Flag{$chr} == 1;
        my $site = "$c[0]\t$c[1]\t$c[2]\t$c[3]";
        my $type = "$c[4],$c[5]";
        $reads{$site}{$type}{$c[-1]} = 1;
}
close IN;
foreach my $site(keys %reads){
        print OUT "$site\t";
        foreach my $type(sort {$a cmp $b} keys %{$reads{$site}}){
                my @Rs = keys %{$reads{$site}{$type}};
                my $n = @Rs;
                print OUT "$type,$n;";
        }
        print OUT "\n";
}
}
#############
close OUT;

Ptime("End");
print "##########End############ perl $0 @ARGV ($optko)\n";
