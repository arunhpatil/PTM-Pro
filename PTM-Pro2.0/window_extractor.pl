#!/usr/bin/perl
$file = $ARGV[0];
$path = $ARGV[1];

open (FILE, "$file") or die $!;
#open (FILE, "/tmp/5cd7d7c463936/Ex_human_PTM_PSMs_prtSite_ptmWindow.txt") or die $!;
while (<FILE>)
{
	chomp ($_);
	@lines = split (/\t/,$_);
	#Step 1 to breakdown |pipes|
	if ($lines[-1] =~ /\|\|/)
	{
		@pipes = split (/\|\|/, $lines[-1]);
		foreach $nowForEach (@pipes)
		{
			if ($nowForEach !~ /\;/)
			{
				if ($nowForEach =~ /^_/)
				{
					$getLen_A = length ($nowForEach);
					$diff_A = 15 - $getLen_A;
					$prefix_A = "_" x $diff_A;
					$windowA = "$prefix_A"."$nowForEach";
					$temp_A = substr($windowA, 7,1);
					$file_A = "$lines[-8]"."_"."$temp_A";
					$filesList{$file_A}="1";
					open (A_OUT, ">>$path/$file_A\.txt") or die $!;
					$windowA =~ tr/[a-z]/[A-Z]/;
					print A_OUT "$windowA\n";
					#print "$windowA\t$lines[-8]\n";
				}
				elsif ($nowForEach =~ /_$/)
				{
					$getLen_B = length ($nowForEach);
					$diff_B = 15 - $getLen_B;
					$suffix_B = "_" x $diff_B;
					$windowB = "$nowForEach"."$suffix_B";
					$temp_B = substr($windowB, 7,1);
					$file_B = "$lines[-8]"."_"."$temp_B";
					$filesList{$file_B}="1";
					open (B_OUT, ">>$path/$file_B\.txt") or die $!;
					$windowB =~ tr/[a-z]/[A-Z]/;
					print B_OUT "$windowB\n";
					#print "$windowB\t$lines[-8]\n";
				}
				else
				{
					$temp_M = substr($nowForEach, 7,1);
					$file_M = "$lines[-8]"."_"."$temp_M";
					$filesList{$file_M}="1";
					open (M_OUT, ">>$path/$file_M\.txt") or die $!;
					$nowForEach =~ tr/[a-z]/[A-Z]/;
					print M_OUT "$nowForEach\n";
					#print "$nowForEach\t$lines[-8]\n";
				}
			}
			else
			{
				@semisplit = split (/\;/, $nowForEach);
				@semisplitModfn = split (/\;/, $lines[-8]);
				for ($i =0; $i<= $#semisplitModfn; $i++)
				{
					if ($semisplit[$i] =~ /^_/)
					{
						$getLen_C = length ($semisplit[$i]);
						$diff_C = 15 - $getLen_C;
						$prefix_C = "_" x $diff_C;
						$windowC = "$prefix_C"."$semisplit[$i]";
						$temp_C = substr($windowC, 7,1);
						$file_C = "$semisplitModfn[$i]"."_"."$temp_C";
						$filesList{$file_C}="1";
						open (C_OUT, ">>$path/$file_C\.txt") or die $!;
						$windowC =~ tr/[a-z]/[A-Z]/;
						print C_OUT "$windowC\n";
						#print "$windowC\t$semisplitModfn[$i]\n";
					}
					elsif ($semisplit[$i] =~ /_$/)
					{
						$getLen_D = length ($semisplit[$i]);
						$diff_D = 15 - $getLen_D;
						$suffix_D = "_" x $diff_D;
						$windowD = "$semisplit[$i]"."$suffix_D";
						$temp_D = substr($windowD, 7,1);
						$file_D = "$semisplitModfn[$i]"."_"."$temp_D";
						$filesList{$file_D}="1";
						open (D_OUT, ">>$path/$file_D\.txt") or die $!;
						$windowD =~ tr/[a-z]/[A-Z]/;
						print D_OUT "$windowD\n";
						#print "$windowD\t$semisplitModfn[$i]\n";
					}
					else
					{
						$temp_N = substr($semisplit[$i], 7,1);
						$file_N = "$semisplitModfn[$i]"."_"."$temp_N";
						$filesList{$file_N}="1";
						open (N_OUT, ">>$path/$file_N\.txt") or die $!;
						$semisplit[$i] =~ tr/[a-z]/[A-Z]/;
						print N_OUT "$semisplit[$i]\n";
						#print "$semisplit[$i]\t$semisplitModfn[$i]\n";
					}
				}
			}
		}
	}
	else
	{
		if ($lines[-1] !~ /\;/)
		{
			if ($lines[-1] =~ /^_/)
			{
				$getLen_E = length ($lines[-1]);
				$diff_E = 15 - $getLen_E;
				$prefix_E = "_" x $diff_E;
				$windowE = "$prefix_E"."$lines[-1]";
				$temp_E = substr($windowE, 7,1);
				$file_E = "$lines[-8]"."_"."$temp_E";
				$filesList{$file_E}="1";
				open (E_OUT, ">>$path/$file_E\.txt") or die $!;
				$windowE =~ tr/[a-z]/[A-Z]/;
				print E_OUT "$windowE\n";
				#print "$windowE\t$lines[-8]\n";
			}
			elsif ($lines[-1] =~ /_$/)
			{
				$getLen_F = length ($lines[-1]);
				$diff_F = 15 - $getLen_F;
				$suffix_F = "_" x $diff_F;
				$windowF = "$lines[-1]"."$suffix_F";
				$temp_F = substr($windowF, 7,1);
				$file_F = "$lines[-8]"."_"."$temp_F";
				$filesList{$file_F}="1";
				open (F_OUT, ">>$path/$file_F\.txt") or die $!;
				$windowF =~ tr/[a-z]/[A-Z]/;
				print F_OUT "$windowF\n";
				#print "$windowF\t$lines[-8]\n";
			}
			else
			{
				$temp_O = substr($lines[-1], 7,1);
				$file_O = "$lines[-8]"."_"."$temp_O";
				$filesList{$file_O}="1";
				open (O_OUT, ">>$path/$file_O\.txt") or die $!;
				$lines[-1] =~ tr/[a-z]/[A-Z]/;
				print O_OUT "$lines[-1]\n";
				#print "$lines[-1]\t$lines[-8]\n";
			}
		}
		else
		{
			@newSemiSplit = split (/\;/, $lines[-1]);
			@newSemiSplitModfn = split (/\;/, $lines[-8]);
			for ($j=0; $j<= $#newSemiSplitModfn; $j++)
			{
				if ($newSemiSplit[$j] =~ /^_/)
				{
					$getLen_G = length ($newSemiSplit[$j]);
					$diff_G = 15 - $getLen_G;
					$prefix_G = "_" x $diff_G;
					$windowG = "$prefix_G"."$newSemiSplit[$j]";
					$temp_G = substr($windowG, 7,1);
					$file_G = "$newSemiSplitModfn[$j]"."_"."$temp_G";
					$filesList{$file_G}="1";
					open (G_OUT, ">>$path/$file_G\.txt") or die $!;
					$windowG =~ tr/[a-z]/[A-Z]/;
					print G_OUT "$windowG\n";
					#print "$windowG\t$newSemiSplitModfn[$j]\n";
				}
				elsif ($newSemiSplit[$j] =~ /_$/)
				{
					$getLen_H = length ($newSemiSplit[$j]);
					$diff_H = 15 - $getLen_H;
					$suffix_H = "_" x $diff_H;
					$windowH = "$newSemiSplit[$j]"."$suffix_H";
					$temp_H = substr($windowH, 7,1);
					$file_H = "$newSemiSplitModfn[$j]"."_"."$temp_H";
					$filesList{$file_H}="1";
					open (H_OUT, ">>$path/$file_H\.txt") or die $!;
					$windowH =~ tr/[a-z]/[A-Z]/;
					print H_OUT "$windowH\n";
					#print "$windowH\t$newSemiSplitModfn[$j]\n";
				}
				else
				{
					$temp_P = substr($newSemiSplit[$j], 7,1);
					$file_P = "$newSemiSplitModfn[$j]"."_"."$temp_P";
					$filesList{$file_P}="1";
					open (P_PUT, ">>$path/$file_P\.txt") or die $!;
					$newSemiSplit[$j] =~ tr/[a-z]/[A-Z]/;
					print P_PUT "$newSemiSplit[$j]\n";
					#print "$newSemiSplit[$j]\t$newSemiSplitModfn[$j]\n";
				}
			}
		}
	}
}
open (MOMOOUT, ">$path/momo_exec.sh") or die $!;
while (($pani, $puri) = each %filesList)
{
#	print "$pani\n";
	print MOMOOUT "/home/ptmproinhousepro/meme/bin/momo motifx --width 15 -oc $path/$pani $path/$pani.txt\n";
}
