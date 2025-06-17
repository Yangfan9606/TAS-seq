#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-17
# version:   2.0
# license:   MIT
# brief:     修改Refseq Fasta染色体名。
#------------------------------------------#

use strict;
use Getopt::Long;
use lib './';
use yangfan;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : Output chr name as Chr1 Chr2 ... / Only output chr with

Usage: .pl [GCF_000001405.40_GRCh38.p14_genomic.fna.gz] -o [out.name]
        -o                      out name
        -help                   output help information

USAGE
GetOptions(\%opts,"a:s","b:s","c:s","d:s","e:s","f:s","g:s","h:s","i:s","j:s","k:s","l:s","m:s","n:s","o:s","p:s","q:s","r:s","s:s","t:s","u:s","v:s","w:s","x:s","y:s","z:s","head:s", "help!");
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
open(IN, "zcat $infile|") or die $!;
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";

open(OUT, ">$outname") or die $!;

#############
my %Lake;
my $flag = 0;
while(<IN>){
        chomp;
    if ($_ =~ /^>/){
        $flag = $_ =~ /^>NC_/?1:0;
        my @c = split/chromosome /,$_;
        my @b = split/,/,$c[1];
        if ($flag == 1){
            if ($_ =~ /mitochondrion/){
                print OUT ">ChrM\n";
            }else{
               print OUT ">Chr$b[0]\n";;
           }
       }
    }else{
        print OUT "$_\n" if $flag == 1;
    }
}
close IN;
close OUT;
#############
Ptime("End");
print "##########End############ perl $0 @ARGV ($optko)\n";
