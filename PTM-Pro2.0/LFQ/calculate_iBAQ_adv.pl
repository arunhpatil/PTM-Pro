#!/usr/bin/perl
use Math::Complex;

$input =$ARGV[0];
$output = $ARGV[1];
$tryptic_file = $ARGV[2];
open (TRY, "$tryptic_file") or die $!;
open (FILE, "$input") or die $!;
open (OUT, ">$output") or die $!;
print OUT "Gene Symbol\tNumber of trptic peptides\tNumber of PSMs\tsummed Intensities\tiBAQ\tlog(iBAQ,2)\n";
while (<FILE>)
{
	chomp ($_);
	$_ =~ s/\r//g;
	$_ =~ s/\"//g;
	@lines = split (/\t/, $_);
	#if ($lines[1] !~ /;/)
	#{
		#@dots = split (/\./,$lines[0]);
		#$pep=$dots[1];
		$pep=$lines[0];
		if ($lines[1] =~ /\#/)
		{
			@genesymbol = split (/\#/, $lines[1]);
			$gene  = "$genesymbol[1]";
		}
		else {$gene = $lines[1]}
		$scan = "$lines[2]";
		$int = "$lines[3]";
		$var = "$pep\*$gene\*$scan";
		#print "$var\t$int\n";
		$preBAQ{$var} = "$int";
	#}

}
while ($tr = <TRY>)
{
	chomp ($tr);
	$tr =~ s/\r//g;
	@hagene = split (/\#/, $tr);
	$tryptic{$hagene[0]} = "$hagene[1]";
	$glength{$hagene[0]} = "$hagene[2]";
}
while (($a, $b) = each %preBAQ)
{
	#print "$a\t$b\n";
	@under = split (/\*/, $a);
	$sum{$under[1]} .= "$b;";
	$psm{$under[1]}++;
}
#NSAF Pre-Calculation
$NSAFsumedInten=0;
while (($x, $y) = each %tryptic)
{
	if (exists $sum{$x})
	{
		@semi_nsaf = split (/\;/, $sum{$x});
		foreach $NSAFsumIn (@semi_nsaf)
		{
			$NSAFsumedInten = $NSAFsumedInten+$NSAFsumIn;
		}
		$NSAFdivision = ($NSAFsumedInten/$glength{$x});
		$preNSAF{$x} = "$NSAFdivision";
		push @sigma, $NSAFdivision;
		$NSAFsumedInten=0;
	}
}
foreach $nsaf_sumVal (@sigma)
{
	$sumedNSAFInten = $sumedNSAFInten+$nsaf_sumVal;
}
foreach $e (sort (keys (%tryptic)))
{
	if (exists $sum{$e})
	{
		@semi = split (/\;/, $sum{$e});
		foreach $sumIn (@semi)
		{
			$sumedInten = $sumedInten+$sumIn;
		}
		$iBAQ = ($sumedInten/$tryptic{$e});
		##$NSAF_final = ($preNSAF{$e}/$sumedNSAFInten);
		if ($iBAQ > 0)
		{
			$log2 = logn($iBAQ, 2);
		}
		else
		{
			$log2 = 0;
		}
		$log2round = sprintf "%.2f", $log2;
		print OUT "$e\t$tryptic{$e}\t$psm{$e}\t$sumedInten\t$iBAQ\t$log2round\n";
		#print OUT "$e\t$tryptic{$e}\t$psm{$e}\t$sumedInten\t$iBAQ\t$log2round\t$NSAF_final\n";
		#print "$e\t$sum{$e}\t-->$sumedInten-->$f\t$iBAQ\n";
		$iBAQ=0; $sumedInten=0;
	}
	else
	{
		print OUT "$e\t-\t-\t-\t0\t0\n";
	}
	
}
#"Annotated Sequence"	"Master Protein Accessions"	"Protein Descriptions"	"First Scan"	"Intensity"
