#!/usr/bin/perl
$file = $ARGV[0];
$outfile = $ARGV[1];
$organism = $ARGV[2];
$probScore = $ARGV[3];
open (FILE, "$file") or die $!;
#open (FILE, "IL-33_pST_pY_sites_corrected_SILAC_at_prt_level_051914.txt") or die $!;
#open (FILE, "Phospho_sites_at_prt_level_051714.txt") or die $!;
#open (FILE, "Phospho_sites_at_prt_level_051914.txt") or die $!;
while (<FILE>)
{
	chomp ($_);
	$_ =~ s/\r//g;
	@lines = split (/\t/, $_);
	#$lines[0] =~ tr/[a-r]/[A-R]/;
	#$lines[0] =~ tr/[u-x]/[U-X]/;
	#$lines[0] =~ tr/z/Z/;	
	$pepBackup = $lines[0];
	#$pepBackup =~ tr/[a-z]/[A-Z]/;
	$uniq_phospho_peps{$lines[0]} = "1";
	#$geneSym = $lines[3];
	#$site = $lines[-3];
	$geneSym = $lines[-4];
	#print "$geneSym\n";
	$site = $lines[-7];
	$modification = $lines[-8];
	$geneSym =~ s/\s//g;
	$geneSym =~ s/\,//g;
	if ($geneSym =~ /\|\|/)
	{
		@genes = split (/\|\|/, $geneSym);
		$geneSym = $genes[0];
	}
	if ($site =~ /\|\|/)
	{
		@pipeSplit= split (/\|\|/, $site);
		#foreach $ys (@pipeSplit)
		#{
		$ys = $pipeSplit[0];	
		$ys =~ s/\s//g;
		if ($ys =~ /\;/)
		{
			@modiSplit = split (/\;/, $modification);
			@sorted_modiSplit = sort @modiSplit;
			$space_attached_modifin = join (":", @sorted_modiSplit);
			$peptide_modification_count{$space_attached_modifin} .= "$pepBackup;";
			@semi = split (/\;/, $ys);
			for ($ms=0; $ms<=$#modiSplit; $ms++)
			{
				$varMod = "$modiSplit[$ms]\#$semi[$ms]";
				$gene{$geneSym} .= "$varMod*";
			}
			#foreach $yz (@semi)
			#{	
				#$gene{$geneSym} .= "$yz*";
				#print "$geneSym\_$yz\n";	
			#}
		} 
		else {
			
			$varMod = "$modification\#$ys";
			$peptide_modification_count{$modification} .= "$pepBackup;";
			$gene{$geneSym} .= "$varMod*";
			#$gene{$geneSym} .= "$ys*";
				#print "$geneSym\_$ys\n";	
		}
		#}
	}
	elsif ($site =~ /\;/) 
	{
		$site =~ s/\s//g;
		#print "$_\t$geneSym\n";
		@modiSplit_1 = split (/\;/, $modification);
		@sorted_modiSplit = sort @modiSplit_1;
		$space_attached_modifin = join (":", @sorted_modiSplit);
		$peptide_modification_count{$space_attached_modifin} .= "$pepBackup;";
		@semi_1 = split (/\;/, $site);
		#print "$_\t$geneSym\t$site\n";
		#@semi_1 = split (/\;/, $site);
		#foreach $s (@semi_1)
		for ($ms1=0; $ms1<=$#modiSplit_1; $ms1++)
		{
			$varMod_1 = "$modiSplit_1[$ms1]\#$semi_1[$ms1]";
			#print "$varMod_1\n";
			$gene{$geneSym} .= "$varMod_1*";	
			#print "$geneSym\_$s\n";	
			#print "$geneSym\_$varMod_1\n";	
		}
	}
	else 
	{
		$varMod_2 = "$modification\#$site";
		$peptide_modification_count{$modification} .= "$pepBackup;";
		$varMod_2 =~ s/\s//g;
		$gene{$geneSym} .= "$varMod_2*";
		#print "$geneSym\_$site\n";	
	}
			#print "$geneSym\n";
}
$totalPS=0;
$size=0;
$SSite=$TSite=$YSite=0;
$s=$t=$y=0;
while (($k,$v)= each %gene)
{
	#print "$k\t$v\n";
	@psites = split (/\*/,$v);
	foreach $star (@psites)
	{
		$unique_hash{$star} = 1;
	}
	@sites_unique = keys(%unique_hash);
	foreach $sitz (@sites_unique)
	{
		#print "$k\_$sitz\n";
		($modiSum, $siteSum) = split (/\#/, $sitz);
		my $color  = substr $siteSum, 0, 1;
		#print "$color<--\n";
		$color =~ tr/[a-z]/[A-Z]/;
		$indvis = "$modiSum\#$color";
		#print "$indvis\n";
		$IndividualUniqueModifications{$indvis}++;
		$TotalUniqueModifications{$modiSum}++;
	}
	%unique_hash=();
#	print "$k\t@unique\n";
	$size = @sites_unique;
#	print "$k\t$v\t$s\t$t\t$y\t$size\n";
	$s=$t=$y=0;
	$totalPS += $size;
}
@unique_gene= keys(%gene);
@sites_unique_pep = keys(%uniq_phospho_peps);
$pcount  = @sites_unique_pep;
$gene_count  = @unique_gene;
#print OUT "1) Total number of unique phosphopeptides\*: $pcount corresponding to $gene_count proteins.\n";
open OUT, ">$outfile" or die $!;
print OUT "Post translational modifications (PTM) data analysis summary:\n\n";
print OUT "SUMMARY:\n";
print OUT "\nOrganism Selected: $organism\n";
print OUT "PTM probability cutoff: $probScore\n";
$single_sdt=$double_sdt=$multi_sdt=0;
#open (FILE, "$file") or die $!;
$serial_nu=1;
while (($pepa, $pepb) = each %peptide_modification_count)
{
	#print "$pepa\t$pepb\n";
	@semi_uniq_peps = split (/\;/, $pepb);
	foreach $findUniq (@semi_uniq_peps)
	{
		$onlyUniqPeps{$findUniq}++;
		$total_onlyUniqPeps{$findUniq}++;
	}
	@semiUniqValue = keys %onlyUniqPeps;
	$PTM_peptides_counts = @semiUniqValue;
	$new_peptide_modification_count{$pepa} = "$PTM_peptides_counts";
	#%onlyUniqPeps=@semiUniqValue=$PTM_peptides_counts="";
	%onlyUniqPeps=();
	#$peptide_modification_count_unique{$pepa} = "";
	#$PTM_peptides_counts = @semi_uniq_peps;
	#print "$pepa\t$PTM_peptides_counts\n";
}
$Total_PTM_peptides_counts = (@semi_total_uniq_peps = keys %total_onlyUniqPeps);
print OUT "\nTotal number of unique peptides with PTM modificatoins: $Total_PTM_peptides_counts corresponding to $gene_count proteins.\n\n";
while (($ttmod, $sumMods) = each %TotalUniqueModifications)
{
	print OUT "\n$serial_nu\). Total unique $ttmod site modificatioins are: $sumMods\n";
	foreach $aminoAcids (A..Z)
	{
		$temp_var = "$ttmod\#$aminoAcids";
		if (exists $IndividualUniqueModifications{$temp_var})
		{
			print OUT "\t$ttmod modificatioin at $aminoAcids is: $IndividualUniqueModifications{$temp_var}\n";
		}
	}
	$serial_nu++;
}
print OUT "\n";
while (($newpepa, $newpepb) = each %new_peptide_modification_count)
{
	#print "$newpepa\t$newpepb\n";
	$percent = ($newpepb/$Total_PTM_peptides_counts)*100;
	$apercent = sprintf("%0.2f", $percent);
	print OUT "$serial_nu\). Peptides with modification \($newpepa\) = $newpepb \($apercent\%\)\n";
	$serial_nu++;
}
 
print OUT "\n\*NOTE: Phosphopeptides containing the same phosphosite is considered as unique. For example: \"VFDKDGNGyISAAELR, VFDKDGNGyISAAELR\" are 
two phosphopeptides but will be considered as unique. On the other hand, if your peptides are phosphorylated at
different sites for example: \"AADEEPDSPsGALQTAAEEEETK, AADEEPDSPSGALQtAAEEEETK\" will be considered as two phosphopeptides instead of one.\n
If you have any furthe queries, please contact us at arun\@ibioinformatics.org.\n\nOnce again Thank you for using \"PTM-pro (Product of In-House scripts)\".\nClick here to analyze new experiment: http://ptm-pro.inhouseprotocols.com/\n"; 
