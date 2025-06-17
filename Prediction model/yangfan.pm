package yangfan;

use strict;
require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 1.00;
@ISA = qw(Exporter);
#@EXPORT = qw(Ptime Overlap Merge median max min avg sum codon2aa readFa GC tmCal interval);
#@EXPORT_OK = qw(Ptime Overlap Merge median max min avg sum codon2aa readFa GC tmCal interval);
@EXPORT = qw(Ptime mean median max min avg sum qianfen GC Tm interval find_diffs fuzzy_pattern combin getnum array_overlap read_fasta_file read_MDZ read_CIGAR);
@EXPORT_OK = qw(Ptime mean median max min avg sum qianfen GC Tm interval find_diffs fuzzy_pattern combin getnum array_overlap read_fasta_file read_MDZ read_CIGAR);

sub Ptime {
	my $time = localtime;
	my ($msg) = @_;
	print "$msg at $time\n";
}

sub mean {
	my (@array) = @_;
	my $sum=0;
	foreach my $i(@array){
		$sum+=$i/@array;
	}
	return $sum;
}

sub median {
	my (@array) = @_;
	my @infoS = sort{$a <=> $b} @array;
	if(@infoS%2 != 0){
		return $infoS[int(@infoS/2)];
	}else{
		return ($infoS[@infoS/2-1]+$infoS[@infoS/2])/2;
	}
}

sub avg{
	my (@array) = @_;
	my $sum=0;
	foreach my $i(@array){
		$sum+=$i/@array;
	}
	return $sum;
}

sub sum{
	my (@array) = @_;
	my $sum=0;
	foreach my $i(@array){
		$sum+=$i;
	}
	return $sum;
}

sub min{
	my (@array) = @_;
	my @infoS = sort {$a<=>$b} @array;
	return $infoS[0];
}

sub max{
	my (@array) = @_;
	my @infoS = sort {$a<=>$b} @array;
	return $infoS[-1];
}

sub qianfen {
	my ($num) = @_;
	my @sp = split/\./,$num;
	while($sp[0] =~ s/(\d)(\d{3})((,\d\d\d)*)$/$1,$2$3/){};
	my $R = $sp[1] eq ""?$sp[0]:"$sp[0].$sp[1]";
	return $R;
}

sub GC{
	my($seq) = @_;
	my @S = split//,$seq;
	my %count;
	$count{'A'}+=0;
	$count{'T'}+=0;
	$count{'G'}+=0;
	$count{'C'}+=0;
	$count{'N'}+=0;
	my $n = 0;
	foreach my $s(@S){
		$count{uc $s}++;
		$n++;
	}
	my $g = $count{'G'};
	my $c = $count{'C'};
	my $gc = ($g+$c)/$n;
	return $gc;
}

###http://biotools.nubic.northwestern.edu/OligoCalc.html
sub Tm{
	my($seq,$na)=@_;
	my @S = split//,$seq;
	$na = 0.05 if($na eq "");
	my $tm=0;
	my %Count;
	my ($a,$t,$g,$c)=(0,0,0,0);
	foreach my $s(@S){
		$Count{uc $s}++;
	}
	$a = $Count{'A'}+0;
	$t = $Count{'T'}+0;
	$g = $Count{'G'}+0;
	$c = $Count{'C'}+0;
	my $len = scalar @S;
#	print "$seq\t$len\tA:$a\tT:$t\tG:$g\tC:$c\n";
	if(@S<13){
		$tm = ($a+$t)*2+($g+$c)*4+16.6*(log($na/0.05)/log(10));
	}elsif(@S<50){
		$tm = 100.5+(41*($g+$c)/($a+$t+$g+$c))-(820/($a+$t+$g+$c))+16.6*(log($na)/log(10));
	}elsif(@S>=50){
		$tm = 81.5+(41*($g+$c)/($a+$t+$g+$c))-(500/($a+$t+$g+$c))+16.6*(log($na)/log(10));
	}
	return $tm;
}

################  查看点是否落在区间内, ($chr,$dot,@b); 格式"1,123,(1-100-101,2-105-107)"
sub interval{
	my ($chr,$dot,@p) = @_;
	my $R;
	foreach my $P(@p){
		my @T = split/\-/,$P;
		my $c = $T[0];
		next if $c ne $chr;
		my $s = $T[1];
		next if $dot < $s;
		my $e = $T[2];
		next if $dot > $e;
		$R .= $R eq ""?"$P":",$P";
	}
	return $R;
}

##################### 寻找2个字符串之间的差异,并输出位置及mutation(1 base) out @ ("1,$M1","7,$M2"...)
sub find_diffs{
	my ($s1,$s2) = @_;
	my @S1 = split//,$s1;
	my @S2 = split//,$s2;
	my @R;
	for (my $i=0;$i<@S1;$i++){
		my $p = $i+1;
		my $b1 = $S1[$i];
		my $b2 = $S2[$i];
		if ($b1 ne $b2){
			my $m = "$p,$b2";
			push @R,$m;
		}
	}
	return @R;
}
###### 返回 允许x个错配的匹配pattern
# my $new_TSO = "TTTCTTATATGGG";
# my $new_TSO_mis1 = fuzzy_pattern($new_TSO,1);
# if ($seq =~ /$new_TSO_mis1/)...
sub fuzzy_pattern {
	my ($original_pattern, $mismatches_allowed) = @_;
	$mismatches_allowed >= 0 or die "Number of mismatches must be greater than or equal to zero\n";
	my $new_pattern = make_approximate($original_pattern, $mismatches_allowed);
	return qr/$new_pattern/;
}
sub make_approximate {
	my ($pattern, $mismatches_allowed) = @_;
	if ($mismatches_allowed == 0) { return $pattern }
	elsif (length($pattern) <= $mismatches_allowed){ $pattern =~ tr/ACTG/./; return $pattern }
	else {
		my ($first, $rest) = $pattern =~ /^(.)(.*)/;
		my $after_match = make_approximate($rest, $mismatches_allowed);
		if ($first =~ /[ACGT]/) {
			my $after_miss = make_approximate($rest, $mismatches_allowed-1);
			return "(?:$first$after_match|.$after_miss)";
		}else { return "$first$after_match" }
	}
}
###  计算组合数
sub combin{
	my ($a,$b) = @_;
	open(R,">temp.R") or die $!;
	print R "library(gtools)\nnrow(combinations($a, $b))\n";
	my $o = `Rscript temp.R`;
	chomp $o;
	my @N = split/ /,$o;
	system("rm -rf temp.R");
	return $N[1];
}
### 获取字符串中的 number ，用于排序 sort { getnum($a) <=> getnum($b) }
sub getnum {
	my $v = shift;
	return( ($v =~ /(\d+)/)[0] || 0);
}

### 数组overlap , 输出uniq 和overlap, 返回值需要用 @{   }; 进行数组化    
#   array_overlap(\@a,\@b);
sub array_overlap{
	my ($first, $second) = @_;
	my @array1 = @{ $first };
	my @array2 = @{ $second };
	my (%count1,%count2);
	$count1{$_}++ for @array1;
	$count2{$_}++ for @array2;
	my (@unique1,@unique2,@overlap);
	foreach my $elem (@array1) {
		if (exists $count2{$elem}) {
			push @overlap, $elem unless grep { $_ eq $elem } @overlap;
		}else{
			push @unique1, $elem;
		}
	}
	foreach my $elem (@array2) {
		if (!exists $count1{$elem}) {
			push @unique2, $elem;
		}
	}
	return (\@unique1, \@unique2, \@overlap);
}

#######################
sub read_fasta_file{
	my ($file) = @_;
	my %sequences;
	my $seq_id;
	open my $fh, '<', $file or die "Cannot open fasta file: $file";
	while (<$fh>) {
		chomp;
		if (/^>(\S+)/) {
			$seq_id = $1;
		}else{
			$sequences{$seq_id} .= uc $_;
		}
	}
	close $fh;
	return \%sequences;
}
###############################
################# read MD:Z: flag output Mutation/Deletion ref bases
################# my $mdz = "10A5^AC4T"; my $start_pos = 1;  my $mdz_results = parse_mdz($mdz, $start_pos); foreach my $result (@$mdz_results) { if ($result->{type} eq 'mutation') { print "Mutation at position $result->{pos}: $result->{base}\n"; } elsif ($result->{type} eq 'deletion') { print "Deletion at position $result->{pos}: $result->{del_seq}\n"; } }
sub read_MDZ{
	my ($mdz_tag, $start_pos) = @_;
	my @results;
	my $pos = $start_pos;
	while ($mdz_tag =~ /(\d+)([A-Z]|\^[A-Z]+)/g) {
		my $match_len = $1;
		my $mutation = $2;
		$pos += $match_len;
		if ($mutation =~ /^[A-Z]$/) {
			push @results, {
				type => 'mutation',
				pos => $pos,
				base => $mutation
			};
			$pos++;
		}elsif($mutation =~ /^\^([A-Z]+)$/) {
			my $del_bases = $1;
			my $del_len = length $del_bases;
			$pos+=$del_len;
			push @results, {
				type => 'deletion',
				pos => $pos,
				del_seq => $del_bases
			};
		}
	}
	return \@results;
}
######################
################# reads CIGAR and output position -> base (no Insertion and Deletion)
################# my $cigar_results = read_CIGAR($CIGAR,$start_pos,$seq);   foreach my $pos (sort { $a <=> $b } keys %$cigar_results) ...  my $base = $cigar_results->{$pos};    
sub read_CIGAR{
	my ($cigar, $start_pos, $sequence) = @_;
	my %pos_map;
	my $ref_pos = $start_pos; # Assuming start_pos is 1-based
	my $seq_idx = 0;
####### Parse CIGAR operations
	my @cig_ops = $cigar =~ /(\d+)([MIDNSHPX=])/g;
	while (my ($len_str, $op) = splice @cig_ops, 0, 2){
		my $len = int($len_str);
		if ($op eq 'M' || $op eq '=' || $op eq 'X') {
			# Handle operations that consume both reference and query
			for (my $i = 0; $i < $len; $i++) {
				my $current_ref_pos = $ref_pos + $i;
				my $current_seq_base = substr($sequence, $seq_idx + $i, 1);
				$pos_map{$current_ref_pos} = $current_seq_base;
			}
			$ref_pos += $len;
			$seq_idx += $len;
		}elsif($op eq 'I' || $op eq 'S') {
			# Insertion or Soft clip: consume query only
			$seq_idx += $len;
		}elsif($op eq 'D' || $op eq 'N') {
			# Deletion or Reference skip: consume reference only
			$ref_pos += $len;
		}
		# Other operations (H, P, etc.) are ignored
	}
	return \%pos_map;
}
