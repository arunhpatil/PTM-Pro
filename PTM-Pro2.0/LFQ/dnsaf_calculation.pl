#!/usr/bin/perl
#echo shell_exec("perl dnsaf_calculation.pl $input $workingDir/$newoutfile $organism");
open (UNG, "$ARGV[0]") or die $!;
open (OUT, ">$ARGV[1]") or die $!;
open (FASTA, "$ARGV[2]") or die $!;
print OUT "Id\tdSAF\tdNSAF\n";
while ($x = <FASTA>)
{
	chomp ($x);
	$x =~ s/\r//g;
	if ($x =~ /^>/)
	{
		@pipes = split (/\|/, $x);
		$gid =$pipes[1];
	}
	else
	{
		$seq{$gid} .= $x;
	}
}
while (($l, $m) = each %seq)
{
	$seqLen = length($m);
	#print "$l\t$m\t$seqLen\n";
	$prtLenSeq{$l} = "$seqLen";
}
while (<UNG>)
{
	chomp ($_);
	$_ =~ s/\r//g;
	$_ =~ s/\"//g;
	@lines = split (/\t/, $_);
	if ($_ =~ /Annotated [S|s]equence/)
	{
		for ($i=0; $i<= $#lines; $i++)
		{
			if ($lines[$i] =~ /PSM count/) { $psmCount = $i;} # print "$lines[$i]\t$i\n";
			if ($lines[$i] =~ /Master Protein Accessions/) { $masterPrtAcc = $i;} # print "$lines[$i]\t$i\n";
		}	
		#splice @lines, 2, 0, 'PSM count';
		$header_join = join ("\t", @lines);
	}
	else
	{
		# yelli jaarito manavu 
		if ($lines[$masterPrtAcc] =~ /\;/)
		{
			$lines[$masterPrtAcc] =~ s/ //g;
			@semi = split (/\;/, $lines[$masterPrtAcc]);
			foreach $j (@semi)
			{
				$shareSpectralCounts{$j} .= "$lines[$psmCount];";
				$proteinList{$j} = "1";
				foreach $z (@semi)
				{
					if ($z ne "$j")
					{
						$proteinM{$j} .= "$z;";
					}
				}
			}
		}
		else
		{
			$proteinList{$lines[$masterPrtAcc]} = "1";
			$uniqSpectralCounts{$lines[$masterPrtAcc]} .= "$lines[$psmCount];";
			#print "$_\n";	
		}
	}
}
$sum_dSAF =0;
while (($abc, $def) = each %proteinList)
{
	$preS = $pre2S = 0;$sum=0; $udSAF=0; $uSpCm = 0; $uSpCi =0; $sSpCi =0; 
	if (exists $shareSpectralCounts{$abc})
	{
		@shared_proteinM = split (/\;/, $proteinM{$abc});
		foreach $denom (@shared_proteinM)
		{
			push @uniqMspectralValues, $uniqSpectralCounts{$denom};
		}
		@uniqMspectralValues_cols = split (/\;/, $uniqMspectralValues);
		foreach $pre2S (@uniqMspectralValues_cols)
		{
			#print "-->@uniqMspectralValues\n";
			$uSpCm += $pre2S;
		}
		@uniqMspectralValues="";
		$lenT_Sc = $prtLenSeq{$abc};
		@uSpCi_array = split (/\;/, $uniqSpectralCounts{$abc});
		@sSpCi_array = split (/\;/, $shareSpectralCounts{$abc});
		foreach $uSi (@uSpCi_array)
		{
			$uSpCi += $uSi;
		}
		foreach $sSi (@sSpCi_array)
		{
			$sSpCi += $sSi;
		}
		if ($lenT_Sc >0 )
		{
			if ($uSpCm == 0) 
			{
				$dSAFi = ( $uSpCi + 0)/ $lenT_Sc;
			}
			else
			{
				$dSAFi = ( $uSpCi + (($uSpCi/$uSpCm)*$sSpCi) ) / $lenT_Sc;
			}
			#if ($abc eq "10800130")
			#{
			#print "$abc\t$dSAFi\t-->$uSpCi-->$uSpCm-->$sSpCi\t$lenT_Sc\n";
			#}
			#print OUT "$abc\t$dSAFi\n";
			$dSAF_gene{$abc} = $dSAFi;
			$sum_dSAF += $dSAFi;
		}
		
	}
	else
	{
		@Uniqu_counts = split (/\;/, $uniqSpectralCounts{$abc});
		$lenT = $prtLenSeq{$abc};
		foreach $preS (@Uniqu_counts) 
		{
			$sum+= $preS;
		}
		if ($lenT > 0)
		{
			$udSAF = ($sum)/$lenT;
			#print "$abc\t$uniqSpectralCounts{$abc}\t$sum\t$lenT\n";
			#print OUT "$abc\t$udSAF\n";
			$dSAF_gene{$abc} = $udSAF;
			$sum_dSAF += $udSAF;
		}
	}
}
while (($geneS, $s) = each %dSAF_gene)
{
	if ($s >0 && $sum_dSAF >0)
	{
		$var_temp = $s/$sum_dSAF;
		print OUT "$geneS\t$s\t$var_temp\n";
		$var_temp = "";
	}
}
#"Checked"	"Confidence"	PSM count	"Identifying Node"	"PSM Ambiguity"	"Annotated Sequence"	"Modifications"	"# Protein Groups"	"# Proteins"	"Master Protein Accessions"	"Protein Accessions"	"Protein Descriptions"	"# Missed Cleavages"	"Charge"	"DeltaScore"	"DeltaCn"	"Rank"	"Search Engine Rank"	"m/z [Da]"	"MH+ [Da]"	"DeltaM [ppm]"	"Deltam/z [Da]"	"Intensity"	"Activation Type"	"MS Order"	"Isolation Interference [%]"	"Ion Inject Time [ms]"	"RT [min]"	"First Scan"	"Spectrum File"	"Ions Matched"	"XCorr"	"Ions Score"	"Identity Strict"	"Identity Relaxed"	"Expectation Value"	"Percolator q-Value"	"Percolator PEP"
