#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
# $Bin为当前路径
use lib '/home/yangfan/Data/Bin/perl_script_my/final/';
use yangfan;
use Statistics::Descriptive;
use List::Util qw(shuffle);
use List::MoreUtils qw(uniq);
#fuzzy_pattern($x,1);
my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : seqlogo (Bits or Prob)

Usage: .pl [.fa] -o [out.name].seqlogo.pdf
        -w                      width
        -h                      height
        -type                   [RNA] (RNA/DNA)
        -m                      methods [P] (P/B) probability / Bits
        -o                      out name
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "h:s", "c:s", "f:s", "w:s", "h:s", "type:s", "m:s", "help!");
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

open(OUT, ">$outname.R") or die $!;
my $W = $opts{w} eq ""?12:$opts{w};
my $H = $opts{h} eq ""?4:$opts{h};
my $type = $opts{type} eq ""?"RNA":$opts{type};
$type = $type eq "RNA"?"rna":"dna";
my $method = $opts{m} eq ""?"P":$opts{m};
$method = $method eq "P"?"prob":"bits";
#############
my %Lake;

while(<IN>){
        chomp;
        my @c = split/\t/;
}
close IN;

print OUT "
require(ggplot2)
require(ggseqlogo)
require(Biostrings)

fasta = '$infile'
fasta_input = readBStringSet(fasta,format='fasta', nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=FALSE)
fasta_input <- as.vector(fasta_input)

pdf('$outname.seqlogo.pdf',wi=$W, he = $H)
ggseqlogo(fasta_input, seq_type='$type', method = '$method')+
#  annotate('rect', xmin = 4.5, xmax = 6.5, ymin = 0, ymax = 1, alpha = 0.5, fill = 'grey50')+
  theme_classic()+
  theme(axis.text = element_text(size=36),
        axis.title.y = element_text(size=40),
        axis.line = element_line(linewidth = 1),
        plot.margin=unit(c(1, 1, 1, 1),'cm')
        )+
  scale_x_continuous(breaks=seq(1,11,1),expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))
dev.off()
";
system("Rscript $outname.R");
#############

close OUT;

Ptime("End");
print "##########End############ perl $0 @ARGV ($optko)\n";
