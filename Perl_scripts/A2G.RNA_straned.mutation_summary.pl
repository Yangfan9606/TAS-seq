#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-18
# version:   1.0
# license:   MIT
# brief:     A2G mutation summary
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

Usage: .pl [K562.CtrlRNA.uniq.Q25bq20d1.AT.sites.mutation.score.noSNP.A2G_mu1.sites] -o [out.name]
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

open(OUT, ">$outname") or die $!;

#############
my %Lake;
# 10      1000822 A       A,+,271;G,+,1;  A2G,+,1;        0.00367647058823529     0.0515375134860788      0.000189476152522349    +
while(<IN>){
        chomp;
        my @c = split/\t/;
        my $site = "$c[0]\t$c[1]";
        my @M = split/;/,$c[3];
        my ($p,$n,$pa,$na)=(0,0,0,0);
        foreach my $m(@M){
                my @N = split/,/,$m;
                $p += $N[2] if $N[1] eq "+";
                $n += $N[2] if $N[1] eq "-";
        }
        my @S = split/;/,$c[4];
        foreach my $s(@S){
                my @N = split/,/,$s;
                if ($N[0] eq "A2G"){
                        $pa += $N[2] if $N[1] eq "+";
                        $na += $N[2] if $N[1] eq "-";
                }
        }
        if ($p != 0){
        my $rp = $pa/$p;
        print OUT "$site\t$c[1]\t$pa\t$rp\t+\t$p\t$c[3]\t$c[4]\n";
        }
        if ($n!=0){
        my $rn = $na/$n;
        print OUT "$site\t$c[1]\t$na\t$rn\t-\t$n\t$c[3]\t$c[4]\n";
        }
}
close IN;

#############
close OUT;

Ptime("End");
print "##########End############ perl $0 @ARGV ($optko)\n";
