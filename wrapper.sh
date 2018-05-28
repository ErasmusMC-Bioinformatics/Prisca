#!/bin/bash

inputFile=$1
outputFile=$2
outputDir=$3
min_freq=$4
min_cells=$5
merge_on="$6"

dir="$(cd "$(dirname "$0")" && pwd)"
mkdir $outputDir


Rscript --verbose $dir/RScript.r $inputFile $outputDir $outputFile $min_freq $min_cells "${merge_on}" 2>&1
cp $dir/jquery-1.11.0.min.js $outputDir
cp $dir/script.js $outputDir
cp $dir/style.css $outputDir
cp $dir/tabber.js $outputDir
mv "$outputFile" "$outputDir/log.html"

echo "<html><center><h1><a href='index.html'>Click here for the results</a></h1>Tip: Open it in a new tab (middle mouse button or right mouse button -> 'open in new tab' on the link above)</center></html>" > $outputFile

cd $outputDir

header="<html><head><script type='text/javascript' src='jquery-1.11.0.min.js'></script><script type='text/javascript' src='tabber.js'></script><script type='text/javascript' src='script.js'></script><link rel='stylesheet' type='text/css' href='style.css'></head><div id='hidden_div' style='display: none;'></div>"
singles=()
pairs_BM_PB=()
pairs_Left_Right=()
pairs_R_Dx=()
while read patient sample1 sample2 type
do
	echo "$patient"
	html="${patient}.html"
	echo "$header" > "$html"
	if [[ "$type" == *pair* ]] ; then
		if [[ "$sample1" == *_BM* ]] || [[ "$sample1" == *_PB* ]] ; then
			pairs_BM_PB+=( "$patient" )
		elif [[ "$sample1" == *_Left* ]] || [[ "$sample1" == *_Right* ]] ; then
			pairs_Left_Right+=( "$patient" )
		else
			pairs_R_Dx+=( "$patient" )
		fi
	else
		singles+=( "$patient" )
	fi
	oldLocus=""
	sample1="$(echo ${sample1} | tr -d '\r' | tr -d '\n')"
	sample2="$(echo ${sample2} | tr -d '\r' | tr -d '\n')"
	tail -n+2 "${patient}_freq.txt" | sed "s/>//" > tmp.txt
	echo "<div class='tabber'>" >> "$html"
	echo "<div class='tabbertab' title='Data frequency'>" >> "$html"
	echo "<table><tr><td style='vertical-align:top;'>" >> "$html"
	echo "<table border = 1 class='result_table summary_table' id='summary_table_${patient}_freq'>" >> "$html"
	echo "<thead><th>Ig/TCR gene rearrangement type</th><th>Proximal gene segment</th><th>Distal gene segment</th><th>Cut off value</th><th>Number of sequences ${patient}_Both</th><th>Number of sequences_$sample1</th><th>Read Count $sample1</th><th>Number of sequences_$sample2</th><th>Read Count $sample2</th><th>Sum number of sequences $patient</th><th>Percentage of sequences ${patient}_both</th></thead>" >> "$html"
	echo "<tbody>" >> "$html"
	scatterplot_tab="<div class='tabbertab' title='Scatter Plots Frequency'><table border='0'><tr>"
	while read locus j_segment v_segment cut_off_value both one read_count1 two read_count2 sum percent locusreadsum1 locusreadsum2
	do
		if [ "$locus" != "$oldLocus" ] ; then
			echo "<tr><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr>" >> "$html"
			echo "<tr><td><b>$locus</b></td>" >> "$html"
		else
			echo "<td></td>" >> "$html"
		fi
		echo "<td>$v_segment</td>" >> "$html"
		echo "<td>$j_segment</td>" >> "$html"
		echo "<td>>$cut_off_value</td>" >> "$html" 
		if [ "$both" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${sample2}_${locus}_${cut_off_value}.txt\", \"$patient\", \"freq\")'>$both</td>" >> "$html"
		else
			echo "<td>$both</td>" >> "$html"
		fi
		if [ "$one" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${locus}_${cut_off_value}.txt\", \"$patient\", \"freq\")'>$one</td>" >> "$html"
		else
			echo "<td>$one</td>" >> "$html"
		fi
		echo "<td>$read_count1</td>" >> "$html"
		if [ "$two" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample2}_${locus}_${cut_off_value}.txt\", \"$patient\", \"freq\")'>$two</td>" >> "$html"
		else
			echo "<td>$two</td>" >> "$html"
		fi
		echo "<td>$read_count2</td>" >> "$html"
		echo "<td>$sum</td>" >> "$html"
		echo "<td>${percent}&#37;</td>" >> "$html"
		echo "</tr>" >> "$html"
		oldLocus="$locus"
		if [ "${cut_off_value}" == "0" ] ; then
			scatterplot_tab="${scatterplot_tab}<td><img src='${patient}_${sample1}_${sample2}_freq_${locus}_scatter.png' /></td>"
		fi
	done < tmp.txt
	echo "</tbody></table>" >> "$html"
	echo "</td><td style='vertical-align:top;'><div id='result_div_${patient}_freq'></div></td></tr></table>" >> "$html"
	echo "</div>" >> "$html"
	echo "<div class='tabbertab' title='Graphs frequency'>" >> "$html"
	echo "<a href='${patient}_freq.png'><img src='${patient}_freq.png' width='1280' height='720' /></a><br />" >> "$html"
	echo "<a href='${patient}_freq_both.png'><img src='${patient}_freq_both.png' width='1280' height='720' /></a><br />" >> "$html"
	echo "<a href='${patient}_percent_freq.png'><img src='${patient}_percent_freq.png' width='1280' height='720' /></a></div>" >> "$html"
	echo "${scatterplot_tab}</tr></table></div>" >> "$html"
	
	tail -n+2 "${patient}_reads.txt" | sed "s/>//" > tmp.txt
	echo "<div class='tabbertab' title='Data reads'>" >> "$html"
	echo "<table><tr><td style='vertical-align:top;'>" >> "$html"
	echo "<table border = 1 class='result_table summary_table' id='summary_table_${patient}_reads'>" >> "$html"
	echo "<thead><th>Ig/TCR gene rearrangement type</th><th>Proximal gene segment</th><th>Distal gene segment</th><th>Cut off value</th><th>Number of sequences ${patient}_Both</th><th>Number of sequences_$sample1</th><th>Read Count $sample1</th><th>Number of sequences_$sample2</th><th>Read Count $sample2</th><th>Sum number of sequences $patient</th><th>Percentage of sequences ${patient}_both</th></thead>" >> "$html"
	echo "<tbody>" >> "$html"
	scatterplot_tab="<div class='tabbertab' title='Scatter Plots Reads'><table border='0'><tr>"
	while read locus j_segment v_segment cut_off_value both one read_count1 two read_count2 sum percent locusreadsum1 locusreadsum2
	do
		if [ "$locus" != "$oldLocus" ] ; then
			echo "<tr><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr>" >> "$html"
			echo "<tr><td><b>$locus</b></td>" >> "$html"
		else
			echo "<td></td>" >> "$html"
		fi
		echo "<td>$v_segment</td>" >> "$html"
		echo "<td>$j_segment</td>" >> "$html"
		echo "<td>>$cut_off_value</td>" >> "$html" 
		if [ "$both" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${sample2}_${locus}_${cut_off_value}.txt\", \"$patient\", \"reads\")'>$both</td>" >> "$html"
		else
			echo "<td>$both</td>" >> "$html"
		fi
		if [ "$one" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${locus}_${cut_off_value}.txt\", \"$patient\", \"reads\")'>$one</td>" >> "$html"
		else
			echo "<td>$one</td>" >> "$html"
		fi
		echo "<td>$read_count1</td>" >> "$html"
		if [ "$two" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample2}_${locus}_${cut_off_value}.txt\", \"$patient\", \"reads\")'>$two</td>" >> "$html"
		else
			echo "<td>$two</td>" >> "$html"
		fi
		echo "<td>$read_count2</td>" >> "$html"
		echo "<td>$sum</td>" >> "$html"
		echo "<td>${percent}&#37;</td>" >> "$html"
		echo "</tr>" >> "$html"
		oldLocus="$locus"
		if [ "${cut_off_value}" == "0" ] ; then
			scatterplot_tab="${scatterplot_tab}<td><img src='${patient}_${sample1}_${sample2}_reads_${locus}_scatter.png' /></td>"
		fi
	done < tmp.txt
	echo "</tbody></table>" >> "$html"
	echo "</td><td style='vertical-align:top;'><div id='result_div_${patient}_reads'></div></td></tr></table>" >> "$html"
	echo "</div>" >> "$html"
	echo "<div class='tabbertab' title='Graphs reads'>" >> "$html"
	echo "<a href='${patient}_reads.png'><img src='${patient}_reads.png' width='1280' height='720' /></a><br />" >> "$html"
	echo "<a href='${patient}_reads_both.png'><img src='${patient}_reads_both.png' width='1280' height='720' /></a><br />" >> "$html"
	echo "<a href='${patient}_percent_reads.png'><img src='${patient}_percent_reads.png' width='1280' height='720' /></a></div>" >> "$html"
	echo "${scatterplot_tab}</tr></table></div>" >> "$html"
	echo "</div>" >> "$html"
	echo "</div>" >> "$html"
	echo "</html>" >> "$html"
done < patients.txt

html="index.html"
echo "<html>" > $html
echo "<table>" >> "$html"
echo "<tr><td><b>Singles (<a href='singles_freq_scatterplot.png'>Frequency scatterplot</a>, <a href='singles_reads_scatterplot.png'>Reads scatterplot</a>):</b></td></tr>" >> "$html"
for patient in "${singles[@]}"
do
	echo "<tr><td><a href='${patient}.html'>$patient</a></td></tr>" >> "$html"
done
echo "<tr><td><b>Pairs (Left & Right):</b></td></tr>" >> "$html"
for patient in "${pairs_Left_Right[@]}"
do
	echo "<tr><td><a href='${patient}.html'>$patient</a></td></tr>" >> "$html"
done
echo "<tr><td><b>Pairs (BM & PB):</b></td></tr>" >> "$html"
for patient in "${pairs_BM_PB[@]}"
do
	echo "<tr><td><a href='${patient}.html'>$patient</a></td></tr>" >> "$html"
done
echo "<tr><td><b>Pairs (Dx & R):</b></td></tr>" >> "$html"
for patient in "${pairs_R_Dx[@]}"
do
	echo "<tr><td><a href='${patient}.html'>$patient</a></td></tr>" >> "$html"
done
echo "<tr><td><b>Triplets:</b></td></tr>" >> "$html"

while read sample1 sample2 sample3
do
	sample1="$(echo ${sample1} | tr -d '\r' | tr -d '\n')"
	sample2="$(echo ${sample2} | tr -d '\r' | tr -d '\n')"
	sample3="$(echo ${sample3} | tr -d '\r' | tr -d '\n')"
	patient="${sample1}_${sample2}_${sample3}"
	echo "$patient"
	html="${patient}.html"
	echo "<tr><td><a href='${patient}.html'>$patient</a></td></tr>" >> "index.html"
	echo "$header" > "$html"
	oldLocus=""
	tail -n+2 "${patient}_freq.txt" | sed "s/>//" > tmp.txt
	echo "<div class='tabber'>" >> "$html"
	echo "<div class='tabbertab' title='Data frequency'>" >> "$html"
	echo "<table><tr><td style='vertical-align:top;'>" >> "$html"
	echo "<table border = 1 class='result_table summary_table' id='summary_table_${patient}_freq'>" >> "$html"
	echo "<thead><th>Ig/TCR gene rearrangement type</th><th>Proximal gene segment</th><th>Distal gene segment</th><th>Cut off value</th><th>Number of sequences ${patient}_All</th><th>Number of sequences_$sample1</th><th>Number of sequences_$sample2</th><th>Number of sequences_$sample3</th><th>Number of sequences_${sample1}_${sample2}</th><th>Number of sequences_${sample1}_${sample3}</th><th>Number of sequences_${sample2}_${sample3}</th></thead>" >> "$html"
	echo "<tbody>" >> "$html"
	scatterplot_tab="<div class='tabbertab' title='Scatter Plots Frequency'><table border='0'><tr>"
	while read locus j_segment v_segment cut_off_value all one two three one_two one_three two_three 
	do
		if [ "$locus" != "$oldLocus" ] ; then
			echo "<tr><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr>" >> "$html"
			echo "<tr><td><b>$locus</b></td>" >> "$html"
		else
			echo "<td></td>" >> "$html"
		fi
		echo "<td>$v_segment</td>" >> "$html"
		echo "<td>$j_segment</td>" >> "$html"
		echo "<td>>$cut_off_value</td>" >> "$html" 
		if [ "$all" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${sample2}_${sample3}_${locus}_${cut_off_value}.txt\", \"$patient\", \"freq\")'>$all</td>" >> "$html"
		else
			echo "<td>$all</td>" >> "$html"
		fi
		if [ "$one" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${locus}_${cut_off_value}.txt\", \"$patient\", \"freq\")'>$one</td>" >> "$html"
		else
			echo "<td>$one</td>" >> "$html"
		fi		
		if [ "$two" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample2}_${locus}_${cut_off_value}.txt\", \"$patient\", \"freq\")'>$two</td>" >> "$html"
		else
			echo "<td>$two</td>" >> "$html"
		fi
		if [ "$three" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample3}_${locus}_${cut_off_value}.txt\", \"$patient\", \"freq\")'>$three</td>" >> "$html"
		else
			echo "<td>$three</td>" >> "$html"
		fi
		
		if [ "${one_two}" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${sample2}_${locus}_${cut_off_value}freq.txt\", \"$patient\", \"freq\")'>${one_two}</td>" >> "$html"
		else
			echo "<td>${one_two}</td>" >> "$html"
		fi
		if [ "${one_three}" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${sample3}_${locus}_${cut_off_value}freq.txt\", \"$patient\", \"freq\")'>${one_three}</td>" >> "$html"
		else
			echo "<td>${one_three}</td>" >> "$html"
		fi
		if [ "${two_three}" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample2}_${sample3}_${locus}_${cut_off_value}freq.txt\", \"$patient\", \"freq\")'>${two_three}</td>" >> "$html"
		else
			echo "<td>${two_three}</td>" >> "$html"
		fi
		
		echo "</tr>" >> "$html"
		oldLocus="$locus"
		if [ "${cut_off_value}" == "0" ] ; then
			scatterplot_tab="${scatterplot_tab}<td><img src='${sample1}_${sample2}_${sample3}_freq_${locus}_scatter.png' /></td>"
		fi
	done < tmp.txt
	echo "</tbody></table>" >> "$html"
	echo "</td><td style='vertical-align:top;'><div id='result_div_${patient}_freq'></div></td></tr></table>" >> "$html"
	echo "</div>" >> "$html"
	echo "<div class='tabbertab' title='Graphs frequency'>" >> "$html"
	echo "<a href='${patient}_freq_total_all.png'><img src='${patient}_freq_total_all.png' width='1280' height='720' /></a><br />" >> "$html"
	echo "<a href='${patient}_freq_indiv_all.png'><img src='${patient}_freq_indiv_all.png' width='1280' height='720' /></a><br /></div>" >> "$html"
	echo "${scatterplot_tab}</tr></table></div>" >> "$html"
	
	tail -n+2 "${patient}_reads.txt" | sed "s/>//" > tmp.txt
	echo "<div class='tabbertab' title='Data reads'>" >> "$html"
	echo "<table><tr><td style='vertical-align:top;'>" >> "$html"
	echo "<table border = 1 class='result_table summary_table' id='summary_table_${patient}_reads'>" >> "$html"
	echo "<thead><th>Ig/TCR gene rearrangement type</th><th>Proximal gene segment</th><th>Distal gene segment</th><th>Cut off value</th><th>Number of sequences ${patient}_All</th><th>Number of sequences_$sample1</th><th>Number of sequences_$sample2</th><th>Number of sequences_$sample3</th><th>Number of sequences_${sample1}_${sample2}</th><th>Number of sequences_${sample1}_${sample3}</th><th>Number of sequences_${sample2}_${sample3}</th></thead>" >> "$html"
	echo "<tbody>" >> "$html"
	scatterplot_tab="<div class='tabbertab' title='Scatter Plots Reads'><table border='0'><tr>"
	while read locus j_segment v_segment cut_off_value all one two three one_two one_three two_three 
	do
		if [ "$locus" != "$oldLocus" ] ; then
			echo "<tr><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr>" >> "$html"
			echo "<tr><td><b>$locus</b></td>" >> "$html"
		else
			echo "<td></td>" >> "$html"
		fi
		echo "<td>$v_segment</td>" >> "$html"
		echo "<td>$j_segment</td>" >> "$html"
		echo "<td>>$cut_off_value</td>" >> "$html" 
		if [ "$all" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${sample2}_${sample3}_${locus}_${cut_off_value}.txt\", \"$patient\", \"reads\")'>$all</td>" >> "$html"
		else
			echo "<td>$all</td>" >> "$html"
		fi
		if [ "$one" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${locus}_${cut_off_value}.txt\", \"$patient\", \"reads\")'>$one</td>" >> "$html"
		else
			echo "<td>$one</td>" >> "$html"
		fi
		if [ "$two" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample2}_${locus}_${cut_off_value}.txt\", \"$patient\", \"reads\")'>$two</td>" >> "$html"
		else
			echo "<td>$two</td>" >> "$html"
		fi
		if [ "$three" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample3}_${locus}_${cut_off_value}.txt\", \"$patient\", \"reads\")'>$three</td>" >> "$html"
		else
			echo "<td>$three</td>" >> "$html"
		fi
		
		if [ "${one_two}" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${sample2}_${locus}_${cut_off_value}reads.txt\", \"$patient\", \"reads\")'>${one_two}</td>" >> "$html"
		else
			echo "<td>${one_two}</td>" >> "$html"
		fi
		if [ "${one_three}" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample1}_${sample3}_${locus}_${cut_off_value}reads.txt\", \"$patient\", \"reads\")'>${one_three}</td>" >> "$html"
		else
			echo "<td>${one_three}</td>" >> "$html"
		fi
		if [ "${two_three}" != "0" ] && [ "$cut_off_value" != "0" ] ; then
			echo "<td data-patient='${patient}' style='cursor:pointer' onclick='javascript:loadfile(\"${sample2}_${sample3}_${locus}_${cut_off_value}reads.txt\", \"$patient\", \"reads\")'>${two_three}</td>" >> "$html"
		else
			echo "<td>${two_three}</td>" >> "$html"
		fi
		
		echo "</tr>" >> "$html"
		oldLocus="$locus"
		if [ "${cut_off_value}" == "0" ] ; then
			scatterplot_tab="${scatterplot_tab}<td><img src='${sample1}_${sample2}_${sample3}_reads_${locus}_scatter.png' /></td>"
		fi
	done < tmp.txt
	echo "</tbody></table>" >> "$html"
	echo "</td><td style='vertical-align:top;'><div id='result_div_${patient}_reads'></div></td></tr></table>" >> "$html"
	echo "</div>" >> "$html"
	echo "<div class='tabbertab' title='Graphs reads'>" >> "$html"
	echo "<a href='${patient}_reads_total_all.png'><img src='${patient}_reads_total_all.png' width='1280' height='720' /></a><br />" >> "$html"
	echo "<a href='${patient}_reads_indiv_all.png'><img src='${patient}_reads_indiv_all.png' width='1280' height='720' /></a><br /></div>" >> "$html"
	echo "${scatterplot_tab}</tr></table></div>" >> "$html"
	echo "</div>" >> "$html"
	echo "</div>" >> "$html"
	echo "</html>" >> "$html"
done < triplets.txt
rm tmp.txt


html="index.html"

echo "</table>" >> "$html"
echo "<a href='log.html'>log</a><br />" >> "$html"
echo "<a href='single_matches.html'>single_matches</a><br />" >> "$html"
echo "<a href='multiple_matches.html'>multiple_matches</a><br />" >> "$html"
echo "</html>" >> "$html"

cp "index.html" "$outputFile"

