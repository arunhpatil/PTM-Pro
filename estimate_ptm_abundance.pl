#!/usr/bin/perl
open (FILE, "ptm_abundance_per_gene2.txt") or die $!;
open (GENE, "gene_list_ptm_abundance.txt") or die $!;
while (<FILE>)
{
	chomp ($_);
	$_ =~ s/\r//g;
	@tabs = split (/\t/, $_);
	$gene = $tabs[2];
	$gene =~ s/\|\|.*//;
	$md = $tabs[1];
	$md =~ s/(\d+)//g;
	$md =~ s/ //g;
	@semis = split (/\;/, $md);
	$siz = @semis -1;
	if ($tabs[0] =~ /k$/ || $tabs[0] =~ /r$/)
	{
		if ($siz != 0)
		{
			$newsiz = $siz-1;
			for ($i =0; $i<= $newsiz; $i++)
			{
				#print "$gene\t$semis[$i]\n";
				$ptm_masterlist{$semis[$i]}="1";
				$var = $gene."#".$semis[$i];
				$estimate{$var}++;
				#print "$i-->$_\t$gene\t$md\t$siz\t$semis[$i]\n";
			}
		}
	}
	else
	{
		for ($j =0; $j<= $siz; $j++)
		{
			#print "$gene\t$semis[$j]\n";
			$var = $gene."#".$semis[$j];
			$estimate{$var}++;
			#print "$j-->$_\t$gene\t$md\t$siz\t$semis[$j]\n";
		}
	}
}
@modifications = keys %ptm_masterlist;
$head = join ("\t", @modifications);
print "Gene\t$head\n";
#while (($x, $y) = each %estimate)
#{
#	print "$x\t$y\n";
#}

while ($gf=<GENE>)
{
	chomp ($gf);
	$gf =~ s/\r//g;
	for ($ptm = 0; $ptm<= $#modifications; $ptm++)
	{
		$req_key = $gf."#".$modifications[$ptm];
		if (exists $estimate{$req_key})
		{
			$present = $estimate{$req_key};
		}
		else
		{
			$present =0;
		}
		push @row_count, $present;
	}
	$join_row_count = join ("\t", @row_count);
	print "$gf\t$join_row_count\n";
	$present=0; @row_count=();
}
#R_Deamidated;
#K_Acetyl
#
