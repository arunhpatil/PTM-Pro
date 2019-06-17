#!/usr/bin/perl
$input = $ARGV[0];
open (FILE,"$input") or die $!;
#open (FILE,"HsRefSeq87.fasta") or die $!;
#open (FILE,"HsRefSeq_RS89.fasta") or die $!;
while (<FILE>)
{
	chomp($_);
	$_ =~ s/\r//g;
	if ($_ =~ /^>/)
	{
		@pipes = split (/\|/, $_);
		@hashes = split (/\#/,$pipes[-1]);
		$gene{$pipes[1]} = $hashes[1];
	}
	else
	{
		$seq{$pipes[1]} .= $_;
	}
}
while (($a, $b) = each %seq)
{
	#print "$a\t$gene{$a}\n$b\n";
	$seqL = length($b);
	#print "$a\t$seqL\n";
	$gene_length{$gene{$a}} .= "$seqL;";
	$b =~ s/([RK])/$1=/g;
	$b =~ s/=P/P/g;
	@split_seq = split (/=/, $b);
	foreach $pep (@split_seq)
	{
		$pep_len = length ($pep);
		if ($pep_len >= 7 && $pep_len <=35)
		{
			$count{$a}++;
		}
		#print "$pep\n";
	}
}
while (($c, $d) = each %count)
{
	#print "$c\t$d\t$gene{$c}\n";
	$ibaq{$gene{$c}} .= "$d;";
#	$gin{$gene{$c}} = 
}
while (($e, $f) = each %ibaq)
{
	@peptides = split (/\;/, $f);
	@sorted = sort {$a <=> $b} @peptides;
	
	@prtlen = split (/\;/, $gene_length{$e});
	@sorted_prtlen = sort {$a <=> $b} @prtlen;
	#print "$e\t@sorted\n";
	#print "$e\t@sorted\t-->$sorted[-1]\n";
	print "$e#$sorted[-1]#$sorted_prtlen[-1]\n";
	#print "Gene Symbol#number_tryptic_peptides#protein_length\n";
	#print "$e\n";
}



# >gi|53828740|ref|NP_001005484.1| NP_001005484.1#OR4F5#79501#olfactory receptor 4F5 [Homo sapiens]
