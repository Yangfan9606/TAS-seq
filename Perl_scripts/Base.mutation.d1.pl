#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-18
# version:   1.0
# license:   MIT
# brief:     统计A/T突变类型。
#------------------------------------------#
use strict;
use Getopt::Long;
use lib './';
use yangfan;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : 统计 每个 点的 detph 及存在的(RNA)突变种类及个数

Usage: .pl [Q25bq20d10.AT.sites] -o [out.name]
        -o                      out name
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "h:s", "c:s", "f:s", "help!");
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
open(IN, $infile) or die $!;
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";
open(OUT,">$outname") or die $!;
print OUT "chr\tsite\tref\treads\tmutation\n";
close OUT;
open(OUT, "|sort -k1,1 -k2,2n >>$outname") or die $!;
open(OUT1, ">$outname.summary") or die $!;
#############
my %Lake;
while(<IN>){
        chomp;
        my @c = split/\t/;
        my $site = "$c[0]\t$c[1]\t$c[2]";
        my @M = split/;/,$c[4];
        my $depth;
        my $o;
        foreach my $m(@M){
                my $ref = $c[2];
                my @T = split/,/,$m;
                $depth += $T[2];
                if ($T[1] eq "-"){
                        $ref =~ tr/AT/TA/;
                }
                $o .= $ref."2".$T[0].",".$T[1].",".$T[2].";" if ($ref ne $T[0]);
        }
        next if $depth < 1;
        $o = "NA" if $o eq "";
        if ($o eq "NA"){
                my @nm = split/;/,$c[4];
                foreach my $NM(@nm){
                        my @on = split/,/,$NM;
                        $Lake{"$on[0],$on[1]"}{num} += $on[2];
                        $Lake{"$on[0],$on[1]"}{sites} += 1;
                }
        }else{
                my @mu = split/;/,$o;
                foreach my $MU(@mu){
                        my @om = split/,/,$MU;
                        $Lake{"$om[0],$om[1]"}{num} += $om[2];
                        $Lake{"$om[0],$om[1]"}{sites} += 1;
                }
        }
        print OUT "$site\t$c[4]\t$o\n";
}
close IN;
foreach my $k(sort {$a cmp $b} keys %Lake){
        print OUT1 "$k\t$Lake{$k}{sites}\t$Lake{$k}{num}\n";
}
#############
close OUT;

Ptime("End");
print "##########End############ perl $0 @ARGV ($optko)\n";
