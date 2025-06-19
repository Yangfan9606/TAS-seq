#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-18
# version:   1.0
# license:   MIT
# brief:     计算A/T的mutation rate。
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

Usage: .pl [K562.8eRNANES.uniq.Q25bq20d10.AT.sites.mutation] [depth] -o [out.name]
        -o                      out name
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "h:s", "c:s", "f:s", "help!");
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
open(IN, $infile) or die $!;
my $total_reads = shift;
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";

open(OUT, ">$outname") or die $!;

#############
my %Lake;

while(<IN>){
        chomp;
        if ($. == 1){
                print OUT "$_\tA2G_rate\tA2G_cpm\tA2G_Mscore\n";
                next;
        }
        my @c = split/\t/;
        my $dp;
        if ($c[4] =~ /A2G/){
                my @D = split/;/,$c[3];
                foreach my $d(@D){
                        my @I = split/,/,$d;
                        $dp += $I[2];
                }
                my @T = split/;/,$c[4];
                foreach my $t(@T){
                        my @I = split/,/,$t;
                        next unless $I[0] eq "A2G";
                        my $m = $I[2];
                        my $rate = $m/$dp;
                        my $cpm = ($m*1000000)/$total_reads;
                        my $R = ($rate)*($cpm);
                        print OUT "$_\t$rate\t$cpm\t$R\n";
                        last;
                }
        }else{
                print OUT "$_\t0\t0\t0\n";
        }
}
close IN;

#############
close OUT;

Ptime("End");
print "##########End############ perl $0 @ARGV ($optko)\n";
