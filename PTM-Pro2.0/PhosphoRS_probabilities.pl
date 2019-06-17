#!/usr/bin/perl
#open (FILE, "Gastric_PST_R2_Sequest_Mascot_All_phosphopep.csv") or die $!;
#perl PhosphoRS_probabilities.pl $filePointer $ptmRsprob $workingDir/$ts $workingDir
$x = $ARGV[0];
$ptmProb = $ARGV[1];
$time_stamp = $ARGV[2];
$path = $ARGV[3];
$mapped = "$time_stamp"."_temp_pepSite_out.txt";
$not_mapped = "$time_stamp"."_temp_pepSite_Not_mapped.txt";
open (FILE, "$x") or die $!;
open (OUT, ">$mapped") or die $!;
open (NOT_OUT, ">$not_mapped") or die $!;
#open (FILE, "Gastric_pST_PSM_table_080414.csv") or die $!;
#open (FILE, "Gastric_PST_R1_Sequest_MAscot_All_phosphopep.csv") or die $!;
while (<FILE>)
{
	chomp ($_);
	@lines = split (/\t/,$_);
	$_ =~ s/\r//g;
	if ($lines[0] !~ /[S|s]equence/)
	{
		#$arr_sizeIn = @lines;
		#print "$arr_sizeIn\n";
		$lines[-1] =~ s/\s//g;
	#	$lines[-1] =~ s/[P|p]hospho//g; #Added this line 071115 to suite PD-2.0
	#	#The following modifications are added on 031017 to suite PD-2.1
	#	$lines[-1] =~ s/[C|c]arbamidomethyl//g; $lines[-1] =~ s/[D|d]eamidated//g; $lines[-1] =~ s/[O|o]xidation//g; 
	#	$lines[-1] =~ s/[A|a]cetyl//g;$lines[-1] =~ s/[A|a]midated//g;$lines[-1] =~ s/[B|b]iotin//g;$lines[-1] =~ s/Carbamyl//g;$lines[-1] =~ s/Carboxymethyl//g;$lines[-1] =~ s/Crotonyl//g;
	#	$lines[-1] =~ s/Dehydrated//g;$lines[-1] =~ s/Dehydro//g;$lines[-1] =~ s/Dioxidation//g;$lines[-1] =~ s/Ethanolyl//g;$lines[-1] =~ s/Formyl//g;
	#	$lines[-1] =~ s/Guanidinyl//g;$lines[-1] =~ s/Methyl//g;$lines[-1] =~ s/Methylthio//g;$lines[-1] =~ s/Myristoyl//g;$lines[-1] =~ s/Palmitoyl//g;$lines[-1] =~ s/Propionamide//g;
	#	$lines[-1] =~ s/Pyridylethyl//g;$lines[-1] =~ s/Pyro-carbamidomethyl//g;$lines[-1] =~ s/Sulfo//g;
		$lines[-1] =~ s/\(/_/g;
		$lines[-1] =~ s/\)//g;
		$lines[0] =~ tr/[a-z]/[A-Z]/;
		if ($lines[-1] ne "")
		{
			#print "$_\n";
			@sites = split (/\;/, $lines[-1]);
			foreach $a (@sites)
			{
				#print "***$a\n";
				@scolon = split (/\:/, $a);
				$modfn = $scolon[0];
				$modfn =~ s/.*\_//;
				int ($scolon[-1]);
				#if ($scolon[-1] >= 10)
				#if ($scolon[-1] >= 1)
				if ($scolon[-1] >= $ptmProb)
				#if ($scolon[-1] >= 75)
				{
					push @only_site, $scolon[0].";";
					push @reqsites, $a.";";
					$modfn_pep = "$modfn"."\#$lines[0]";
					$unambi{$modfn_pep}++;
					$goodPepSeq{$lines[0]} = "1";
				}
				else
				{
					if ($modfn !~ /isoform/ && $modfn !~ /-/ && $modfn !~ /[I|i]nconclusive/)
					{
						$modfn_pep_ambi = "$modfn"."\#$lines[0]";
						$ambi{$modfn_pep_ambi} = "1";
						#$ambi{$modfn_pep_ambi} .= "$lines[0];";
					}
				}
			}
			$rsize = @reqsites;
			#print ">>>$rsize\n";
			if ($rsize != 0)
			{
				
				print OUT "$_\t@only_site\n";
			}
			else 
			{
				#print "$_\t@sites\tFALSE\n";
				print NOT_OUT "$_\t-\n";
			}
			@reqsites=(); @only_site =();
		}
		else
		{
			print NOT_OUT "$_\t-\n";
		}
	}
	#print "@sites\n";
	#print "$lines[8]\n";
}
#isoform
open BAROUT, ">$path/bar_inc.py" or die $!;
while (($unama, $unamb) = each %unambi)
{
	@hash_una = split (/\#/, $unama);
	$unambiModPepCount{$hash_una[0]}++;
}
while (($a, $b) = each %ambi)
{
	@semiSplitamb = split (/\#/, $a);
	if (!exists $goodPepSeq{$semiSplitamb[1]})
	{
		$ambiguous{$semiSplitamb[0]}++;
	}
#unambiguous_peptides = (20, 35, 30, 35, 27) 
#ambiguous_peptides = (25, 32, 34, 20, 25)
#ptms = ('G1', 'G2', 'G3', 'G4', 'G5')
}
while (($final_uambi, $final) = each %ambiguous)
{
#	print BAROUT "#----> $final_uambi\t$final\n";
	if (!exists $unambiModPepCount{$final_uambi})
	{
		$unambiCountVal = "0";
	}
	else
	{
		$unambiCountVal = "$unambiModPepCount{$final_uambi}";
	}
	$pre_ptmList .= "\'$final_uambi\'\,";
	$pre_ambiguous_peptides .= "$final\,";
	$pre_unambiguous_peptides  .= "$unambiCountVal\,";
}
@ambi_keys = keys %unambiModPepCount;
foreach $forgotten_keys (@ambi_keys)
{
	if (!exists $ambiguous{$forgotten_keys})
	{
		$pre_ptmList .= "\'$forgotten_keys\'\,";
		$pre_ambiguous_peptides .= "0\,";
		$pre_unambiguous_peptides  .= "$unambiModPepCount{$forgotten_keys}\,";	
	}
}
if ($pre_ptmList =~ /^\,/) {$pre_ptmList =~ s/^,//;}
if ($pre_ptmList =~ /\,$/) {$pre_ptmList =~ s/,$//;}
if ($pre_ambiguous_peptides =~ /^\,/){$pre_ambiguous_peptides =~ s/^,//;}
if ($pre_ambiguous_peptides =~ /\,$/){$pre_ambiguous_peptides =~ s/,$//;}
if ($pre_unambiguous_peptides =~ /^\,/){$pre_unambiguous_peptides =~ s/^,//;}
if ($pre_unambiguous_peptides =~ /\,$/){$pre_unambiguous_peptides =~ s/,$//;}

print BAROUT "ptms = \($pre_ptmList\)\n";
print BAROUT "unambiguous_peptides = \($pre_unambiguous_peptides\)\n";
print BAROUT "ambiguous_peptides= \($pre_ambiguous_peptides\)\n";

print NOT_OUT "\n\nThose peptides did not map and you should be concerned about:\n\n";
