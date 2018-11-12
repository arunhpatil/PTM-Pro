#!/usr/bin/perl
#open (FILE, "Gastric_PST_R2_Sequest_Mascot_All_phosphopep.csv") or die $!;
$x = $ARGV[0];
$ptmProb = $ARGV[1];
$time_stamp = $ARGV[2];
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
	$_ =~ s/\r//g;
	@lines = split (/\t/,$_);
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
	if ($lines[-1] ne " ")
	{
		@sites = split (/\;/, $lines[-1]);
		foreach $a (@sites)
		{
			#print "***$a\n";
			@scolon = split (/\:/, $a);
			int ($scolon[-1]);
			#if ($scolon[-1] >= 10)
			#if ($scolon[-1] >= 1)
			if ($scolon[-1] >= $ptmProb)
			#if ($scolon[-1] >= 75)
			{
				push @only_site, $scolon[0].";";
				push @reqsites, $a.";";
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
	#print "@sites\n";
	#print "$lines[8]\n";
}
print NOT_OUT "\n\nThose peptides did not map and you should be concerned about:\n\n";
