#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-18
# version:   1.0
# license:   MIT
# brief:     Calculate delta mutation rate
#------------------------------------------#
use strict;
use Getopt::Long;
use lib './';
use yangfan;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : 根据mutation site，统计点上CTL的reads数, mutation rate 为 raw rate

Usage: .pl [K562.8eRNANES.Filter.A2G.sites.bed] [K562.CtrlRNA.Filter.A2G.sites.bed] -o [out.name]
        -o                      out name
        -control                [] (CTL uniq reads number)
        -treat                  [] (NES/NLS uniq reads number)
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "h:s", "c:s", "f:s","control:s", "treat:s", "help!");
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
if ($infile =~ /\.gz$/){open(IN, "zcat $infile|") or die $!;
}else{open(IN, $infile) or die $!;}
my $infile=shift;
if ($infile =~ /\.gz$/){open(IN1, "zcat $infile|") or die $!;
}else{open(IN1, $infile) or die $!;}
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";

open(OUT, ">$outname") or die $!;

#############
my %Lake;
#Usage: .pl [K562.8eRNANES.Filter.A2G.sites.bed] [K562.CtrlRNA.Filter.A2G.sites.bed] -o [out.name]
while(<IN>){
        chomp;
        my @c = split/\t/;
        my $site = "$c[0]\t$c[1]\t$c[5]";
        $Lake{$site}{rate_T} = $c[4];
        $Lake{$site}{mu_T} = $c[3];
        $Lake{$site}{total_T} = $c[6];
        $Lake{$site}{flag} = 1 if $c[8] =~ /A2G/;
}
close IN;
while(<IN1>){
        chomp;
        my @c = split/\t/;
        my $site = "$c[0]\t$c[1]\t$c[5]";
        $Lake{$site}{rate_C} = $c[4];
        $Lake{$site}{mu_C} = $c[3];
        $Lake{$site}{total_C} = $c[6];
        $Lake{$site}{flag} = 1 if $c[8] =~ /A2G/;
}
close IN1;
my $psudo=1e-8;
my $K562_CTL=19403342;
my $K562_NES=12109824;
my $K562_NLS=10409173;
my $T293_CTL=28663205;
my $T293_NLS=13740344;
my $treat = $opts{treat};
my $control = $opts{control};
print "normalize use treat = $treat\n";
print OUT "chr\tsite\tstrand\tT_rate\tT_m\tT_total\tNor_m\tNor_total\tC_rate\tC_m\tC_total\tdelta_rate\tlog2fc_rate\n";
foreach my $site(keys %Lake){
        next if $Lake{$site}{flag} eq "";
        my $T_r = $Lake{$site}{rate_T}+0;
        my $C_r = $Lake{$site}{rate_C}+0;
        my $T_m = $Lake{$site}{mu_T}+0;
        my $C_m = $Lake{$site}{mu_C}+0;
        my $T_t = $Lake{$site}{total_T}+0;
        my $C_t = $Lake{$site}{total_C}+0;
        next if ($T_m == 0 && $C_m == 0);
        my $Nor_m = ($T_m*$control)/$treat;
        my $Nor_t = ($T_t*$control)/$treat;
        my $delta = $T_r - $C_r;
        my $log2fc = log(($T_r+$psudo)/($C_r+$psudo))/log(2);
        print OUT "$site\t$T_r\t$T_m\t$T_t\t$Nor_m\t$Nor_t\t$C_r\t$C_m\t$C_t\t$delta\t$log2fc\n";
}
############
close OUT;

Ptime("End");
print "##########End############ perl $0 @ARGV ($optko)\n";
