#!/usr/bin/perl
$input = $ARGV[0];
@under = split (/\_/, $input);
$output = "$under[0]"."_header_file_temp.txt";
open (FILE, "$input") or die $!;
while (<FILE>)
{
	chomp($_);
	@lines = split (/\t/, $_);
	$line = join ("\t", @lines);
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	$newline = "$line\t" . "PTM_Site (Peptide)\tPeptideSequence\tSite_modification\tPTM_Site (Peptide)\tPTM_Modification\tPTM_Site (Protein)\tProteinGroupAccession\tNP_Accession\tGeneSymbol\tGeneId\tProteinDescription\tPTM_window\tPhospho_sitePlusEvidence\n";
}
open (OUT, ">$output") or die $!;
#open (OUT, ">header_file_temp.txt") or die $!;
print OUT $newline;
