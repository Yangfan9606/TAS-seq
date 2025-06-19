#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-18
# version:   3.0
# license:   MIT
# brief:     Filter SNP
#------------------------------------------#
use strict;
use Getopt::Long;
use lib './';
use yangfan;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : Stranded    A2G

Usage: .pl [human.SNP_2018.bed.gz] [all.sample.site_rate.filter_Q25.out] -o [out.name]
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
open(IN,"zcat $infile|") or die $!;
my $infile=shift;
open(IN1, $infile) or die $!;
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";

open(OUT, ">$outname") or die $!;
open(OUT1, ">$outname.snp") or die $!;

#############
my %Lake;

while(<IN>){
        chomp;
        my @c = split/\t/;
        my @S = split/-/,$c[4];
        if ($S[0] =~ /A/){
                if ($S[1] =~ /G/){
                        $Lake{"$c[0],$c[1],+"} = $_;
                }
        }
        if ($S[0] =~ /T/){
                if ($S[1] =~ /C/){
                        $Lake{"$c[0],$c[1],-"} = $_;
                }
        }
}
close IN;
while(<IN1>){
        chomp;
        if ($.==1){
                print OUT "$_\n";
                next;
        }
        my @c = split/\t/;
        my $flag = 0;
        if ($c[4] =~ /A2G/){
                my @T = split/;/,$c[4];
                foreach my $t(@T){
                        my @I = split/,/,$t;
                        if ($I[0] eq "A2G"){
                                if ($Lake{"$c[0],$c[1],$I[1]"} ne ""){
                                        $flag = 1;
                                        print OUT1 "$_\t".$Lake{"$c[0],$c[1],$I[1]"}."\n";
                                }
                        }
                }
        }
        print OUT "$_\n" if $flag == 0;
}
close IN1;

#############
close OUT;

Ptime("End");
print "##########End############ perl $0 @ARGV ($optko)\n";
