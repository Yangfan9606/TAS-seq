#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-18
# version:   3.0
# license:   MIT
# brief:     Cal Fisher's exact test p-value and p-adj
#------------------------------------------#
use strict;
use Getopt::Long;
use lib './';
use yangfan;
use Statistics::Descriptive;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : Fisher's exact test , 指定4列输入, 依次计算每行的 Fisher's Pvalue 和 P ajd(BH or FDR方法)   !!!!  注意 title 不要有一样name的
        For each reigion/site/point...
         input sample1_edited_reads + sample2_edited_reads
         input sample1_total_reads + sample2_total_reads
         当真0假设(true null hypotheses)比例较高时，FDR方法更好一点(FDR容许一定的假阳性率的存在, BH 更保守)

Usage: .pl [] -o [out.name]
        -h                      header or not [0] (0/1) (head 中不可出现 ; 号，会被自动替换为 _)
        -s1                     column of sample 1 edited value [2]
        -t1                     column of sample 1 total value  [3]
        -s2                     column of sample 2 edited value [4]
        -t2                     column of sample 2 total value  [5]
        -name                   name of s1,t1,s2,t2 [T_count,T_total,C_count,C_total]
        -o                      out name
        -p                      padjust method [BH] (BH/FDR)
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "h:s", "c:s", "f:s", "s1:s", "t1:s", "s2:s", "t2:s","name:s", "help!");
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
#open(OUT, ">$outname") or die $!;
my $head = $opts{h} == 1?"T":"F";
if ($head eq "T"){
        open(HOUT,">$outname.head") or die $!;
}
#######################################
my $Tc = $opts{s1};
my $Tt = $opts{t1};
my $Cc = $opts{s2};
my $Ct = $opts{t2};
my $dfname = "T_count,T_total,C_count,C_total";
my @name = split/,/,$dfname;
if ($opts{name} ne ""){
        @name = split/,/,$opts{name};
        for(my $i=0;$i<@name;$i++){
                if ($name[$i] =~ /^[0-9]/){
                        $name[$i] = "x.$name[$i]";
                }
        }
}
my $PM = $opts{p} eq "FDR"?"FDR":"BH";
my ($col_name,$out_name);
while(<IN>){
        chomp;
        my @c = split/\t/;
        for(my $i=0;$i<@c;$i++){
                my $j = $i+1;
                if ($head eq "T"){
                        $c[$i] =~ s/;/_/g;
                        if ($c[$i] =~ /^[0-9]/){
                                print "change name of $c[$i] to x.$c[$i]\n";
                                $c[$i] = "x.$c[$i]";
                        }
                }
                my $tempN = $head eq "T"?$c[$i]:"V$j";
                $tempN = $name[0] if $j == $Tc;
                $tempN = $name[1] if $j == $Tt;
                $tempN = $name[2] if $j == $Cc;
                $tempN = $name[3] if $j == $Ct;
                if ($head eq "T"){
                        print HOUT "$tempN\t";
                }
                $col_name .= "\"$tempN\",";
                $out_name .= "$tempN=data\$$tempN,";
        }
        if ($head eq "T"){
                print HOUT "PValue\tP.adj\n";
        }
        last;
}
close IN;
print "-------$col_name\n-------$out_name\n";
chop ($col_name,$out_name);
#$out_name .= "PValue=p_value\$PValue,P_adj=adjusted_p_value";
# Write the R script to perform Fisher's exact test, Bonferroni correction, and FDR control
my $r_script = "$outname.fisher_test.R";
open(my $fh, ">", $r_script) or die "Cannot open R script file: $!";
print $fh "
data=read.table(\"$infile\",header = $head, sep = \"\\t\")
colnames(data)=c($col_name)
p_value=data.frame($out_name,PValue = NA)
for(i in 1:nrow(data)){
        contingency_table = matrix(c(data\$$name[0]\[i],data\$$name[1]\[i],data\$$name[2]\[i],data\$$name[3]\[i]),nrow = 2,byrow = T)
        fisher_result = fisher.test(contingency_table)
        p_value\$PValue[i] = fisher_result\$p.value
}
adjusted_p_value = p.adjust(p_value\$PValue,method = \"$PM\")
results=data.frame($out_name,PValue=p_value\$PValue,P_adj.$PM=adjusted_p_value)
write.table(results, \"$outname.temp\",sep=\"\\t\", row.names = F, quote = F)
";
print "Running...R\n";
# Run the R script
system("Rscript $r_script");
if ($head == 1){
        system("cat $outname.head $outname.temp > $outname");
        system("rm -rf $outname.head $outname.temp");
}else{
        system("mv $outname.temp $outname");
}
#############
close OUT;

Ptime("End");
print "##########End############ perl $0 @ARGV ($optko)\n";
