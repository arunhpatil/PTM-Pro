#!/bin/sh
#
# This file is part of the InHouse protocols distribution (http://inhouseprotocols.com).
# Copyright (c) 2018, Arun H. Patil. 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
clear
date=$(date +%Y%m%d%H%M%S)
rand=$(cat /dev/urandom | tr -cd [:alnum:] | head -c 4)
ts=$date$rand

echo "PTM-Pro: Post-translational modifications - profiling"
echo "Copyright (c) 2018, Arun H. Patil"
echo
ptmRsprob=$1
database=$2 
organism=$3 
#organism="human_subset_080315.txt" 
newdate=$(date +%Y%m%d)
filename=$4;
outfile=$(sed -e "s/.txt/_prtSite_ptmWindow.txt/g" <<< $filename)
newoutfile=$(sed -e "s/.txt/_prtSite_ptmWindow_$newdate\.txt/g" <<< $filename)
summary=$(sed -e "s/.txt/_Summary_$newdate\.txt/g" <<< $filename)
new_hf=$ts"_header_file.txt"
new_hf_out=$ts"_header_file_temp.txt"
head -1 $filename > $new_hf
perl header_parser.pl $new_hf 

echo "Filtering PTM-RS probability at $ptmRsprob"
perl PhosphoRS_probabilities.pl $filename $ptmRsprob $ts

echo "Profiling PTMs from $database database"
echo 
perl Fetch_phosphosite.pl $outfile $database $ts & PID=$!
echo "perl Fetch_phosphosite.pl $outfile $database $ts "
printf "["
while kill -0 $PID 2> /dev/null; do
    printf  "▓"
    sleep 1
done
printf "] "
echo
#perl found_in_phsopho_site_plus.pl $outfile $organism $newoutfile
#echo "perl found_in_phsopho_site_plus.pl $outfile $organism $newoutfile"
#perl found_in_phsopho_site_plus.pl $outfile $organism $newoutfile
perl psp.pl $outfile $organism $newoutfile
chmod 777 $new_hf_out $newoutfile
cat $newoutfile >> $new_hf_out
cp $new_hf_out $newoutfile

#chmod 777 $new_hf_out $outfile #--> HERE
#cat $outfile >> $new_hf_out #--> HERE
#cp $new_hf_out $newoutfile 
perl ptm_summary.pl $outfile $summary $database $ptmRsprob & HEAD=$!
echo "perl ptm_summary.pl $outfile $summary $database $ptmRsprob"
echo "Summarizing the PTMs!"
printf "["
while kill -0 $HEAD 2> /dev/null; do
    printf  "▓"
    sleep 1
done
printf "] "
echo
#echo "$new_hf"
#echo "$new_hf_out"
#tar="$file_names[0]".".tar";
temp_name=$(sed -e "s/.txt/_NOT_MAPPED\.txt/g" <<< $filename)
del_mapped_file=$ts"_temp_pepSite_out.txt";
notMAPPED=$ts"_temp_pepSite_Not_mapped.txt";
mv $notMAPPED $temp_name
sleep 2
#rm $del_mapped_file
#rm $new_hf
#rm $new_hf_out
#rm $outfile
echo "Execution Complete!!!"
echo
echo "Cite: Patil, A. H., Datta, K. K., et al., Dissecting Candida pathobiology: Post-translational modifications on the Candida tropicalis proteome. 2018. OMICS: A Journal of Integrative Biology."
echo
