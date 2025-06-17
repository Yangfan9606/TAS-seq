#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-17
# version:   2.0
# license:   MIT
# brief:     使用bed文件，提取Fasta中的序列。
#------------------------------------------#

use strict;
use Getopt::Long;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********# 

Program : Get Sequence(site/region/.bed file) from file.FA file

Usage: .pl [Genome.fa/RNA.fa] [-t chr/-t gene)] [chr/gene_id] [Pos(100)/Pos1-Pos2(100-200)]
           
		  /Data/Database/hg38/GRCh38.fa
	-s		(0/1) default [1] 使用bed文件时，-链是否自动转成正链 [即是否考虑strand信息]
	-f		file.bed (use -f will output a new FA file)
	-t		chr/gene (chr:use Genome.fa as FA_file; gene: use .RNA.fa as FA_file[default:chr])
	-help		output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "t:s", "f:s", "s:s", "help!");
##########
die $usage if ( @ARGV==0 || defined($opts{"help"}));

###################################################
#                  START                          #
###################################################
my $infile=shift;
open(INfile, $infile) or die $!;
my $Qname=shift;
print "Qname $Qname\n";
my $Qsite=shift;
print "Qsite $Qsite\n";
my @Qtest=split/-/,$Qsite;
my $Qflag=@Qtest;
my $S1;
my $S2;
if ($Qflag==2){
	$S1=$Qtest[0];
	$S2=$Qtest[1];
	print "Qflag $Qflag\tS1 $S1\tS2 $S2\n";
	die "Wrong Start-End \n" if $S1 > $S2;
}elsif($Qflag != 1 && $opts{f} eq ""){
	die "Site was wrong! Use onePos(100) or region(100-200)\n";
}
if ($opts{f} ne ""){
open(Qfile, $opts{f}) or die $!;
open(OUT, ">$opts{f}.fa") or die $!;
}
my $Tflag;
if ($opts{t} eq "gene"){
	$Tflag=0;
}elsif($opts{t} eq "chr" || $opts{t} eq ""){
	$Tflag=1;
}else{
	die "Wrong Type! use -t chr or -t gene \n";
}
my $Sflag = $opts{s} eq ""?1:$opts{s};

my $optko;
foreach my $opt(keys %opts){
	$optko .= " -$opt $opts{$opt}";
}
print "##########Start############ perl $0 @ARGV ($optko)\n";
Ptime("Start");
#############
my %Lake;

if ($Tflag == 1){
	print "find chr pos~~~:\n";
	my $ChrSame=0;
		my $count=0;
	if ($opts{f} ne ""){
		my %BedChr;
		while(<Qfile>){
			chomp;
			my @temp=split/\t/;
			$BedChr{$temp[0]}{'flag'}=1;
		}
		my $CHR;
		my $dflag;
		while(<INfile>){
			chomp;
			if (/^>/){
			my @temp=split/>/;
			$dflag=0;
			next unless $BedChr{$temp[1]}{'flag'}==1;
			my @temp = split/>/;
			$CHR=$temp[1];
			$dflag=1;
			next;
			}
			next unless $dflag==1;
			$Lake{$CHR} .= $_;
		}
		close Qfile;
		open (Qfile, $opts{f}) or die $!;
		while(<Qfile>){
			chomp;
			my @temp=split/\t/;
			my $st = $temp[1]-1;
			my $len = $temp[2]-$st;
			my $strand = $temp[5];
			my $OutHead=">";
			foreach my $c(@temp){
				$OutHead .= "$c:";
			}
			my $OutSeq= substr($Lake{$temp[0]},$st,$len);
			if ($strand eq "-" && $Sflag == 1){
				$OutSeq = reverse $OutSeq;
				$OutSeq =~ tr/ATGCatgc/TACGtacg/;
			}
			print OUT "$OutHead\n$OutSeq\n";
		}
		close Qfile;
		close OUT;
	}else{
while(<INfile>){
	last if $ChrSame==2;
	chomp;
	my @temp = split/>/;
	my $tempnum = @temp;
	if ($tempnum == 2){
		$ChrSame=2 if $ChrSame==1;
		next unless $temp[1] eq $Qname;
		$ChrSame=1;
		print "Chr\tSite\tSeq\n";
		next;
	}
	next unless $ChrSame==1;
	$Lake{$Qsite} .= $_;
}
if ($Qflag==1){
	my $st = $Qsite-1;
	my $OutSeq= substr($Lake{$Qsite},$st,1);
	print "$Qname\t$Qsite\t$OutSeq\n";}
elsif($Qflag==2){
	my $st = $S1-1;
	my $len =$S2-$st;
	my $OutSeq= substr($Lake{$Qsite},$st,$len);
	print "$Qname\t$Qsite:\n\t\t$OutSeq\n";}
	}
}else{
	print "find gene pos~~~:\n";
	my $GeneSame;
	my $LakeFlag;
	if ($opts{f} ne ""){
		my %BedChr;
		while(<Qfile>){
			chomp;
			my @temp=split/\t/;
			$BedChr{$temp[0]}{'flag'}=1;
		}
		my $CHR;
		my $dflag;
		while(<INfile>){
			chomp;
			if (/^>/){
			my @temp=split/>/;
			$dflag=0;
			next unless $BedChr{$temp[1]}{'flag'}==1;
			my @temp = split/>/;
			$CHR=$temp[1];
			$dflag=1;
			next;
			}
			next unless $dflag==1;
			$Lake{$CHR} .= $_;
		}
		close Qfile;
		open (Qfile, $opts{f}) or die $!;
		while(<Qfile>){
			chomp;
			my @temp=split/\t/;
			my $st = $temp[1]-1;
			my $len = $temp[2]-$st;
			my $OutHead=">";
			foreach my $c(@temp){
				$OutHead .= "$c:";
			}
			my $OutSeq= substr($Lake{$temp[0]},$st,$len);
			print OUT "$OutHead\n$OutSeq\n";
		}
		close Qfile;
		close OUT;
	}else{
while(<INfile>){
	last if $GeneSame == 2;
	chomp;
	my @temp = split/>/;
	my $tempnum = @temp;
	if ($tempnum == 2){
		$GeneSame =2 if $GeneSame==1;
		next unless $temp[1] eq $Qname;
		$GeneSame=1;
		print "GeneID\tSite\tSeq\n";
		next;
	}
	next unless $GeneSame==1;
	$Lake{$Qsite} .= $_;
}
if ($Qflag==1){
	my $Qt=$Qsite-1;
	my $OutSeq= substr($Lake{$Qsite},$Qt,1);
	print "$Qname\t$Qsite\t$OutSeq\n";}
elsif($Qflag==2){
	my $Qt=$S1-1;
	my $len =$S2-$S1+1;
	my $OutSeq= substr($Lake{$Qsite},$Qt,$len);
	print "$Qname\t$Qsite:\n\t\t$OutSeq\n";}
}
}

#############

close INfile;
#close OUTlog;

sub Ptime {
	my $time = localtime;
	my ($msg) = @_;
	print "$msg at $time\n";
}
Ptime("End");
print "##########End############\n";
