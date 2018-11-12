#!/usr/bin/perl
$input_file = $ARGV[0];
$organism = $ARGV[1];
$output_file = $ARGV[2];
#print "Input file name: $input_file\tOrganism: $organism\tOutfile name: $output_file\n";
open (FILE, "$input_file") or die $!;
open (PSP, "PSP_datasets.txt") or die $!;
#open (PSP, "$organism") or die $!;
while (<PSP>)
{
	chomp ($_);
	@lines = split (/\t/, $_);
	if ($lines[7] eq "$organism")
	{
		$psp_gene = $lines[1];
		$psp_modif = $lines[0];
		$psp_win = $lines[10];
		$psp_win =~ tr/[a-z]/[A-Z]/;
		$pspKey = "$psp_gene\*$psp_modif\*$psp_win";
		$pspWindow{$pspKey}="1";
		if ($psp_modif eq "Acetylation")
		{
			$pspKey = "$psp_gene\*Acetyl\*$psp_win";
			$pspWindow{$pspKey}="1";
		}
		if ($psp_modif eq "Methylation")
		{
			$pspKey = "$psp_gene\*Methyl\*$psp_win";
			$pspWindow{$pspKey}="1";
		}
		if ($psp_modif eq "Phosphorylation")
		{
			$pspKey = "$psp_gene\*Phospho\*$psp_win";
			#print "####$pspKey###\n";
			$pspWindow{$pspKey}="1";
		}
		if ($psp_modif eq "Sumoylation")
		{
			$pspKey = "$psp_gene\*Sumoyl\*$psp_win";
			$pspWindow{$pspKey}="1";
		}
		if ($psp_modif eq "Ubiquitination")
		{
			$pspKey = "$psp_gene\*Ubiquitin\*$psp_win";
			$pspWindow{$pspKey}="1";
		}
		if ($psp_modif eq "O-GlcNAc")
		{
			$pspKey = "$psp_gene\*HexNac\*$psp_win";
			$pspWindow{$pspKey}="1";
		}
		if ($psp_modif eq "O-GalNAc")
		{
			$pspKey = "$psp_gene\*HexN\*$psp_win";
			$pspWindow{$pspKey}="1";
		}
	}
}
#MODIFICATION	GENE	PROTEIN	ACC_ID	HU_CHR_LOC	MOD_RSD	SITE_GRP_ID	ORGANISM	MW_kD	DOMAIN	SITE_+/-7_AA	LT_LIT	MS_LIT	MS_CST	CST_CAT#
#Acetylation	YWHAB	14-3-3 beta	P31946	20q13.12	K5-ac	33347661	human	28.08		___MtMDksELVQkA		2		
open OUT, ">$output_file" or die $!;
while ($x = <FILE>)
{
	chomp($x);
	$x =~ s/\r//g;
	$noSemi= 0;$pipeYesSemi=0;
	@il33_temp = split (/\t/, $x);
	#Initializing the array
	$il33[0] = $il33_temp[-4]; #Gene 
	$il33[1] = $il33_temp[-8]; #Modification type
	#$il33[1] = $il33_temp[-7]; #Site
	$il33[2] = $il33_temp[-1]; #Window
	$il33[2] =~ tr/[a-z]/[A-Z]/;
	$il33[0] =~ s/\s//g;
	$il33[1] =~ s/\s//g;
	$il33[2] =~ s/\s//g;	
	$il33[0] =~ s/\,//g;
	#$pspWindow{$pspKey};
	if ($il33[0] =~ /\|\|/)
	{
		@pipe_gene = split (/\|\|/, $il33[0]);	
		$z_size = @pipe_gene;
		@pipe_seq = split (/\|\|/, $il33[2]);	
		for ($z=0; $z < $z_size; $z++)
		{
			if ($il33[1] =~ /;/)
			{ #Even sequence has pipes 
				$pipeYesSemi=1;
				@pipe_site = split (/\;/, $il33[1]);
				$modfn_times = @pipe_site;
				@semi_seq = split (/\;/, $pipe_seq[$z]);
				for ($mt=0; $mt<$modfn_times; $mt++)
				{
					$geneQuery_win = "$pipe_gene[$z]\*$pipe_site[$mt]\*$semi_seq[$mt]";
					#print "-->$geneQuery_win<-- \n";
					if (exists $pspWindow{$geneQuery_win}) {
						push @short_pipeSemiArr, "Yes";
					}else{
						push @short_pipeSemiArr, "No";
					}
				}
				$short_pipeSemiArr_value = join ("\;", @short_pipeSemiArr);
				if ($short_pipeSemiArr_value =~ /^;/) {$short_pipeSemiArr_value =~ s/\;//;}
				push @long_pipeSemiArr, $short_pipeSemiArr_value; 
				$short_pipeSemiArr_value=""; @short_pipeSemiArr="";
				
			}
			else
			{
				$noSemi= 1;
				$query_win = "$pipe_gene[$z]\*$il33[1]\*$pipe_seq[$z]";
				if (exists $pspWindow{$query_win}) {
					push @short_pipeArr, "Yes";
				}else{
					push @short_pipeArr, "No";
				}
			}
		}
		if ($noSemi ==1)
		{
			$short_pipeArr_value = join ("||", @short_pipeArr);
			if ($short_pipeArr_value =~ /^\|\|/) {$short_pipeArr_value =~ s/\|\|//;}
			print OUT "$x\t$short_pipeArr_value\n"; $short_pipeArr_value=""; @short_pipeArr="";
		}
		elsif ($pipeYesSemi ==1 )
		{
			$long_pipeSemiArr_value = join ("||", @long_pipeSemiArr);
			if ($long_pipeSemiArr_value =~ /^\|\|/) {$long_pipeSemiArr_value =~ s/\|\|//;}
			print OUT "$x\t$long_pipeSemiArr_value\n"; $long_pipeSemiArr_value=""; @long_pipeSemiArr="";
		}
	}
	else 
	{
		if ($il33[1] =~ /;/)
		{
			@pipe_site = split (/\;/, $il33[1]);
			$modfn_times = @pipe_site;
			@seq_semi = split (/\;/, $il33[2]);
			$semi_size = @seq_semi;
			for ($mt=0; $mt<$modfn_times; $mt++)
			{
				$invar = "$il33[0]\*$pipe_site[$mt]\*$seq_semi[$mt]";
				if (exists $pspWindow{$invar})
				{
					#print "--->$invar<--YES--\n";
					push @short_semiArr, "Yes";
				}
				else
				{
					push @short_semiArr, "No";
				}
			}
			$short_semiArr_value = join ("\;", @short_semiArr);
			if ($short_semiArr_value =~ /^;/) {$short_semiArr_value =~ s/\;//;}
			print OUT "$x\t$short_semiArr_value\n"; $short_semiArr_value=""; @short_semiArr="";
		}
		else
		{
			$query_win = "$il33[0]\*$il33[1]\*$il33[2]";	
			#print "-->$query_win\n";
			if (exists $pspWindow{$query_win})
			{
				print OUT "$x\tYes\n";
			}
			else
			{
				print OUT "$x\tNo\n";
			}
		}
	}
}
#==> subset_IL33.txt <==
#Cdk2|| Cdk1     T14;Y15||T14;Y15        VEKIGEGtYGVVYKA;EKIGEGTyGVVYKAK||IEKIGEGtYGVVYKG;EKIGEGTyGVVYKGR


#==> subset_phosphorylation_site_dataset.txt <==
#1110035H17Rik	T11	mouse	PPPGSRstVAQSPPQ

#Acetylation
#Methylation
#O-GalNAc
#O-GlcNAc
#Phosphorylation
#Sumoylation
#Ubiquitination
