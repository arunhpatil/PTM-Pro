#!/usr/bin/perl
#Arun H Patil
#051714
$time_stamp = $ARGV[2];
$x = "$time_stamp"."_temp_pepSite_out.txt";
$not_mapped = "$time_stamp"."_temp_pepSite_Not_mapped.txt";
$outf = $ARGV[0];
$org = $ARGV[1];
open (NEW_OUT, ">>$not_mapped") or die $!;
open (FILE, "$org") or die $!;
while (<FILE>)
{
	chomp ($_);
	$_ =~ s/\r//g;
	if ($_ =~ /^\>/)
	{
		@headsplit = split (/\|/, $_);
		$headkey = $headsplit[1];
		$accession{$headkey} = "$headsplit[-1]";
	}
	else 
	{
		$fetch{$headkey} .= "$_";
	}
}
$out_file= "$outf";
open (OUT, ">$out_file") or die $!;
open (PHS, "$x") or die $!;
while ($x= <PHS>)
{
	chomp ($x);
	$x =~ s/\r//g;
	@lines = split (/\t/, $x);	
	$peptide = $lines[0];
	if ($peptide =~ /\./) #Based on PD 2.1, [R].sDIDSPPITAR.[N], the peptide should be split across the dots. 
	{
		@pure_peptide = split (/\./, $peptide);
		$peptide = $pure_peptide[1];
	}
	$peptide =~ tr /[a-z]/[A-Z]/;
	$site = $lines[-1];
	if ($site =~ /\;/)
	{
		@semi = split (/\;/,$site);
		for ($i=0; $i<=$#semi; $i++)
		{
			$semi[$i] =~ /(\d+)/;
			$sdigit = $1;
			$sdigit = $sdigit-1;
			$semi[$i] =~ /(\D+)/;
			$svalue = $1;
			$svalue =~ tr/[A-Z]/[a-z]/;
			#print "$digit\t$value\n";
			substr $peptide, $sdigit, 1, $svalue;
			$peptide =~ s/\s//g;
		}
	}
	$peptide =~ s/\s//g;
	$lines[0] = $peptide;
	$x = join ("\t", @lines);
	if ($site =~ /,/)
	{
		$site =~ s/\, /\;/;
	}
	$pga = $lines[1];
	$pga =~ s/\,//g;
	$pga =~ s/\s//g;
	$pga =~ s/\;/\|\|/g;
	$raw_peptide = $peptide;
	$peptide =~ tr/[a-z]/[A-Z]/;
	if ($pga =~ /\|\|/)
	{
		@pga_split = split (/\|\|/, $pga);
		#@pga_split = split (/\|\|\s/, $pga);
		#@pipe_site_split = split (/\;/,$site);
		#@pipe_site_split = split (/\|\|/,$site);
		$pipe_size = @pga_split;
		$flag=3;
		for($i=0; $i <$pipe_size; $i++)
		{
			#print "-->$i\t$pga_split[$i]<--\n";
			#print "-->$fetch{$pga_split[$i]}<--\n";
			if (exists $fetch{$pga_split[$i]})
			{
				$subject = $fetch{$pga_split[$i]};
				$total_strings = length($subject);
				$result = index($subject, $peptide);
				#$result=$result;
				#$req_site = $pipe_site_split[$i];
				if ($site =~ /\;/)
				{
					@sites = split (/\;/, $site);
					#@sites = split (/\;/, $req_site);
					foreach $sits (@sites)
					{
						$sits =~ /(\d+)/;
			                        $digit = $1;
       				                $new_site= $result+$digit;
       			                        $sits =~ /(\D+)/;
		       	                        $value = $1;
			                        $new = $value.$new_site;
						push @newsites, $new;
						$reverse = $new_site-$total_strings-8;
						$R = substr $subject, $reverse, 7;
						$F = substr $subject, $new_site, 7;
						$temp_value = $value;
						$temp_value =~ tr/[A-Z]/[a-z]/;
						$lr = length ($R);
						$lf = length ($F);
						if ($lr < 7 ) 
						{
							$R = "_".$R;
						}
						if ($lf < 7)
						{
							$F = $F."_";
						}
						push @news_seqs_semi, $R.$temp_value.$F;
					}
					$rec_seqs_semi = join (';', @news_seqs_semi);
					$rec = join(';',@newsites);
					push @news_seqs, $rec_seqs_semi;
					push @news_array, $rec;
					#print "---> $i\t$rec\n";
					@newsites = ();
					@news_seqs_semi = ();
					$flag =1;
				}
				else 
				{
					#$req_site =~ /(\d+)/;
					$site =~ /(\d+)/;
					$digit = $1;
					$new_site= $result+$digit;
					#$req_site =~ /(\D+)/;
					$site =~ /(\D+)/;
					$value = $1;
					$new = $value.$new_site;
					#print "$sspsite\t$new\n";
					#print "$subject\t$peptide\t$site\t$new\n";
					$reverse = $new_site-$total_strings-8;
					$R = substr $subject, $reverse, 7;
					$F = substr $subject, $new_site, 7;
					push @news_array, $new;
					$temp_value = $value;
					$temp_value =~ tr/[A-Z]/[a-z]/;
					$lr = length ($R);
					$lf = length ($F);
					if ($lr < 7 ) 
					{
						$R = "_".$R;
					}
					if ($lf < 7)
					{
						$F = $F."_";
					}
					push @news_seqs, $R.$temp_value.$F;
					#$flag =2;
				}
				@header = split (/\#/, $accession{$pga_split[$i]});
				for ($z=0; $z<=$#header; $z=$z+4)
				{
					push @NPS, $header[$z];
					push @gene,$header[$z+1];
					push @gid, $header[$z+2];
					push @dis, $header[$z+3];
				}
			}
			else {print NEW_OUT "$x\tNOT_MAPPED\n";}
		}
		$rec_news = join('||',@news_array);
		$rec_seqs1 = join('||',@news_seqs);
		@news_array=();
		@news_seqs = ();
		#print "$pga\t$site\t$rec_news\t$rec_seqs1\n";
		$NP_A = join ("\|\|", @NPS);
		$GeneSym = join ("\|\|", @gene);
		$GeneId = join ("\|\|", @gid);
		$discrep = join ("\|\|", @dis);
		$NP_A =~ s/\s//g;
		@header=@gene=@gid=@dis=@NPS=();
		#@header = split (/\#/, $accession{$pga});
		$rec_seqs1 =~ s/\s//g;
		#print "A: $site\n";
		#B: S17_Phospho;
		#B: M3_Oxidation; S4_Phospho;
		#A: S6_Phospho;
		
		@Asemi_ptmSites = split (/\;/, $site);
		foreach $As_site (@Asemi_ptmSites)
		{
			$As_site =~ s/\s//g;
			($asite, $amodi) = split (/\_/, $As_site);
			push @A_khaliSite, $asite;
			push @A_khaliModi, $amodi;
		}
		$A_khaliModiJoined = join (";", @A_khaliModi);
		$A_khaliSiteJoined = join (";", @A_khaliSite);
		$A_khaliModiJoined =~ s/^;//; $A_khaliSiteJoined =~ s/^;//;
		@A_khaliSite=@A_khaliModi="";
		#print "$A_khaliSiteJoined\t$A_khaliModiJoined\n";
		print OUT "$x\t$raw_peptide\t$site\t$A_khaliSiteJoined\t$A_khaliModiJoined\t$rec_news\t$pga\t$NP_A\t$GeneSym\t$GeneId\t$discrep\t$rec_seqs1\n";
		#print OUT "$x\t$raw_peptide\t$site\t$rec_news\t$pga\t$header[1]\t$header[2]\t$header[3]\t$rec_seqs1\n";
		#print OUT "$x\t$raw_peptide\t$site\t$rec_news\t$pga\t$rec_seqs1\n";
	}
	else
	{	
		#print "------> $pga\n";
		if (exists $fetch{$pga})
		{
			$subject = $fetch{$pga};
			$total_strings = length($subject);
			$result = index($subject, $peptide);
			$result=$result;
			if ($site =~ /\;/)
			{
				@sites = split (/\;/, $site);
				foreach $sits (@sites)
				{
					$sits =~ /(\d+)/;
	                                $digit = $1;
       		                        $new_site= $result+$digit;
       	                        	$sits =~ /(\D+)/;
                                	$value = $1;
	                                $new = $value.$new_site;
        	                        #print "$subject\t$peptide\t$site\t$result\t$new_site\t$digit\t$new\n";
                	                #print "$subject\t$peptide\t$pga\t$site\t$new\n";
					$reverse = $new_site-$total_strings-8;
					$R = substr $subject, $reverse, 7;
					$F = substr $subject, $new_site, 7;
					push @newsites, $new;
					$temp_value = $value;
					$temp_value =~ tr/[A-Z]/[a-z]/;
					$lr = length ($R);
					$lf = length ($F);
					if ($lr < 7 ) 
					{
						$R = "_".$R;
					}
					if ($lf < 7)
					{
						$F = $F."_";
					}
					push @seqs, $R.$temp_value.$F;
				}
				$rec = join(';',@newsites);
				$rec_seq = join(';',@seqs);
                	        #print "$subject\t$peptide\t$site\t@newsites\n";
                	        #print "$subject\t$peptide\t$pga\t$site\t$rec\n";
				@header = split (/\#/, $accession{$pga});
				$rec_seq =~ s/\s//g;
				#print "B: $site\n";
				@Bsemi_ptmSites = split (/\;/, $site);
				foreach $Bs_site (@Bsemi_ptmSites)
				{
					$Bs_site =~ s/\s//g;
					($bsite, $bmodi) = split (/\_/, $Bs_site);
					push @B_khaliSite, $bsite;
					push @B_khaliModi, $bmodi;
				}
				$B_khaliModiJoined = join (";", @B_khaliModi);
				$B_khaliSiteJoined = join (";", @B_khaliSite);
				$B_khaliModiJoined =~ s/^;//; $B_khaliSiteJoined =~ s/^;//;
				@B_khaliSite=@B_khaliModi="";
				#print "$B_khaliSiteJoined\t$B_khaliModiJoined\n";
				print OUT "$x\t$raw_peptide\t$site\t$B_khaliSiteJoined\t$B_khaliModiJoined\t$rec\t$pga\t$header[0]\t$header[1]\t$header[2]\t$header[3]\t$rec_seq\n";
				@newsites = ();
				@seqs = ();
			}
			else 
			{
				$site =~ /(\d+)/;
				$digit = $1;
				$new_site= $result+$digit;
				$site =~ /(\D+)/;
				$value = $1;
				$new = $value.$new_site;
				#print "$subject\t$peptide\t$site\t$result\t$new_site\t$digit\t$new\n";
				#print "$subject\t$peptide\t$pga\t$site\t$new\n";
				$reverse = $new_site-$total_strings-8;
				$R = substr $subject, $reverse, 7;
				$F = substr $subject, $new_site, 7;
				$temp_value = $value;
				$temp_value =~ tr/[A-Z]/[a-z]/;
				$lr = length ($R);
				$lf = length ($F);
				if ($lr < 7 ) 
				{
					$R = "_".$R;
				}
				if ($lf < 7)
				{
					$F = $F."_";
				}
				#print "--> $subject\t$peptide\t$R$value$F\n";
				@header = split (/\#/, $accession{$pga});
				#print "C: $site\n";
				@Csemi_ptmSites = split (/\;/, $site);
				foreach $Cs_site (@Csemi_ptmSites)
				{
					$Cs_site =~ s/\s//g;
					($csite, $cmodi) = split (/\_/, $Cs_site);
					push @C_khaliSite, $csite;
					push @C_khaliModi, $cmodi;
				}
				$C_khaliModiJoined = join (";", @C_khaliModi);
				$C_khaliSiteJoined = join (";", @C_khaliSite);
				$C_khaliModiJoined =~ s/^;//; $C_khaliSiteJoined =~ s/^;//;
				@C_khaliSite=@C_khaliModi="";
				#print "$C_khaliSiteJoined\t$C_khaliModiJoined\n";
				print OUT "$x\t$raw_peptide\t$site\t$C_khaliSiteJoined\t$C_khaliModiJoined\t$new\t$pga\t$header[0]\t$header[1]\t$header[2]\t$header[3]\t$R$temp_value$F\n";
			}
		}
		else {print NEW_OUT "$x\tNOT_MAPPED\n";}
	}
}
