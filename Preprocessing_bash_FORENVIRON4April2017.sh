#!/bin/bash

# Deals with multiple folders of mapping files and others.
# Trimming files to $trimlength bp using fastx.

# EDITED SPECIFICALLY FOR ENVIRON-- 99% REPSET

doMED=False
doOTUPICK=True

# ==================================================
# Sources and summary of contents

source VAR
	# print project_name
	# print fasta_files
	# print Demultiplexed
	# print PATHfiles
	# print mapping_fps --> mappingfile
	# print rep_set --> reftreefasta
	# print reference_tree --> reftreemap
	# print aligned_rep_set --> reftreealigned
	# print length_reads --> trimlength
	# print tree_building_method treemethod
	# print align_gap_filter --> gap
	# print entropy_filter --> entropy
	# print min_num_reads --> minimumcount
	# print MED_value MED
source mappingfilenames.txt
	# $mappingfilenames = comma delimited list of mappingfiles
	# 
# ==================================================
# Setting variables; making fasta files and mapping files

# Concatenate the all fna files to make a single file for trimming and clipping
echo "Concatenating fasta files..."
echo ""
cat ./FASTAFILES/*.fna > ./FASTAFILES/ALL_seqs.fna
for F in ./FASTAFILES/seqs*.fna; do
	rm $F
	done

# Merge the mapping files using a QIIME script
echo "Merging mapping files"
echo ""
merge_mapping_files.py -m $mappingfilenames -o ./MAPPINGFILES/Merged_mapping_file.txt
for F in ./MAPPINGFILES/mappingfile*.txt; do
	rm $F
	done


# ==================================================
# TRIMMING AND CLIPPING

mkdir Trimmed_Quality_Filtered_MSL
cd FASTAFILES

# Using the $trimlength inputted before, use fastx trimmer and clipper to trim the seq file from multiple split libraries and then filter out sequences less than the minimum sequence length.
echo "Trimming and clipping sequences using fastx..."

if [ ${doMED} == 'True' ]; then
 fastx_trimmer -l $length_reads -i ALL_seqs.fna | fastx_clipper -v -a NNNNNNNNNNNNN  -l $length_reads -o ../Trimmed_Quality_Filtered_MSL/seqs_trim_clip.fna >> ../LOG
fi

if [ ${doMED} == 'False' ]; then
	fastx_clipper -v -a NNNNNNNNNNNNN -l $length_reads -i ALL_seqs.fna -o ../Trimmed_Quality_Filtered_MSL/seqs_trim_clip.fna >> ../LOG
fi 

echo "Trimming/Clipping Complete"
echo "" >> ../LOG

# Counting sequences
echo "POST-FASTX SEQUENCE COUNT" >> ../LOG
echo "seqs_trim_clip.fna sequence count" | tee -a ../LOG
grep -c ">" ../Trimmed_Quality_Filtered_MSL/seqs_trim_clip.fna | tee -a ../LOG
echo "" >> ../LOG
# back to project home dir
cd ..

if [ ${doMED} == 'True' ]; then


	# ==================================================
	# FORMATTING

	# In home directory.
	# Now, we want to alter the seq.fna file so that the format fits with MED requirements.
	# MED likes it to be SampleXXX_ReadXXX, so this is what we'll do.
	# Note that this ONLY works for halifax data so far-- you'll likely have to change this for other formats. 
	# The halifax data gets us fasta headers look like this:
	# >JulyZosInnerH_1972165 M02352:23:000000000-ANBEB:1:2119:15710:25291 1:N:0:0 orig_bc=ACCTACTTGTCT new_bc=ACCTACTTGTCT bc_diffs=0
	# CCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGTGGTTAAAAAGCTCGTAGTTGGATCTCAACAGCCTTGAAGCGGTTAACTTTGTAGTTTTACTGCTTTATAGGCTGTGCTTTTGCTGGAGGTTTAGGGGTGCTCTTGATCGAGTGTCTCTGGATGCCAGCAGGTTTACTTTGAAAAAATTAGAGTGCTCAAAGCAGGCTAAAATGCCTGAATATTCGTGCATGGAATAATAGAATAGGA
	# where JulyZosInnerH is the sample name and 1972165 is the read number
	# In order to prep the sequence headers for MED, we want to get rid of everything after the first space, and add "Read" before the read number

	cd ./Trimmed_Quality_Filtered_MSL
	echo "Prepping fasta headers for MED..."
	cat seqs_trim_clip.fna | awk -F' ' '{ st = index($0," ");print $1}' | awk -F_ '{ print ($2!="" ? $1"_Read"$2 : $1) }' > ../seqs_MED.fna
	echo "done!"
	# TODO: the "Read" part of the above string is actually not mandatory, but this wasn't known at the time the script was written. Currently this has no effect on the output, so it's not a high priority. Should be fixed eventually.
	# the command above splits each line on " " and then returns the first element only, then splits that first element on _ and returns both but adds "_Read" between them. 
	# it doesn't modify the sequence lines since they don't satisfy any of the splitting criteria, and the insertion in the second part of the command is conditional on a second element existing

	#checking output fasta for correctly formatted headers and dumping output to LOG file for review
	echo "Recording sample and read counts in LOG file..."
	echo "SAMPLE AND READ COUNTS FOR MED" >> ../LOG
	o-get-sample-info-from-fasta ../seqs_MED.fna >> ../LOG
	echo "" >> ../LOG
	echo "done!"

	cd ..
	# Cleaning up and moving output to its own directory
	echo "Cleaning up workspace..."
	rm -rf ./Trimmed_Quality_Filtered_MSL
	mkdir ./MED
	mv seqs_MED.fna ./MED/seqs_MED.fna

	# In the home directory

	#***** TO DO: CLEAN UP EXTRA FILES IN HOME DIRECTORY*****

	# ==================================================

	# SECTION THREE: MED decompose and transposing
	# This section decomposes the fasta files using the $minimumentropy input
	# Also, transposes the resulting MATRIX-COUNT.txt file into a QIIME-friendly format
	# The last part changes the format of NOE-REPRESENTATIVES.fasta to be compatible with QIIME
	# (explained below)

	# ==================================================


	# Enter the minimum entropy. This is the only variable in this script that changes for MED (Minimum Entropy Decomposition)
	if [ $MED_value == None ] 
	then
		echo "Enter Minimum Substantive Abundance (for MED, see MED documentation for additional details) -- the default is total reads divided by 5000."
		read MED_value
		echo "minimum substantive abundance: $minimumentropy" >> LOG
		echo "" >> LOG
	else
		minimumentropy=$MED_value
	fi

	

	# Now, decompose fasta file. This is a single line.

	echo "MED is decomposing fasta file..."

	decompose ./MED/seqs_MED.fna -M $MED_value -o ./MED/decompose-$MED_value

	echo "MED run completed!"

	# After decomposition, we need to transpose the MATRIX-COUNT.txt file because it is in the wrong format.
	# I was too lazy to write my own script for transposing so I copied one from stackoverflow and tested it before including it below
	# This part of the script was copy and pasted from stackoverflow

	echo "Transposing MATRIX_COUNT.txt..."
	awk '
	{ 
		for (i=1; i<=NF; i++)  {
			a[NR,i] = $i
		}
	}
	NF>p { p = NF }
	END {    
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str"\t"a[i,j];
			}
			print str
		}
	}' ./MED/decompose-$MED_value/MATRIX-COUNT.txt > MATRIX-COUNT_tmp.txt

	# A few more alterations-- change the top left name from samples to #OTU ID, removing leading zeros from OTU IDs

	sed 's/samples/#OTU ID/g' MATRIX-COUNT_tmp.txt | sed 's/^0*//' > MATRIX-COUNT_transposed.txt

	rm MATRIX-COUNT_tmp.txt # Deletes intermediate file
	rm ./MED/seqs_MED.fna

	# Also, remove the pesky "|size:XXX" at the end of every OTU ID in the NODE_REPRESENTATIVES.fasta file.
	# This is necessary because:
	# 1. Adding metadata later on requires that the observation_metadata doesn't have the "|size:XXX"
	# 2. The tree cannot have "|size:XXX" otherwise the tree tips will not match up with the OTU table itself.

	echo "Modifying MED fasta headers for use with QIIME..."
	sed 's/|size:[0-9]*$//g' ./MED/decompose-$MED_value/NODE-REPRESENTATIVES.fasta | sed 's/^>0*/>/' > NODE_REP_FORDOWNSTREAM.fasta
	echo "done!"
	# ==================================================

	# SECTION FOUR: Making a tree and OTU Table
	# This last section goes back to QIIME to make a phylotree and the OTU table
	# You can specify which tree you want
	# It also records some of the intermediate results in the LOG
	# Finally, it adds observation metadata to the file.
	# The resulting OTU table and tree have been tested to be compatible for: alpha_rarefaction,beta_diversity_through_plots,summarize_taxa_through_plots 

	# ==================================================

	# Assign Taxonomy using previous filepaths of references database. Need to use the unaligned database.
	echo "QIIME STEPS BEGIN HERE"
	echo "Assigning taxonomy..."
	assign_taxonomy.py -i NODE_REP_FORDOWNSTREAM.fasta -o assign_taxonomy -r $rep_set -t $reference_tree

	# Align sequences. This step uses the aligned reference database you supplied.

	echo "Aligning sequences..."
	align_seqs.py -i NODE_REP_FORDOWNSTREAM.fasta -t $aligned_rep_set -o aligned_seqs

	# Now, writing the newly aligned sequences and failures to the log for easy access later.

	echo "" >> LOG
	echo "Aligned failures: " >> LOG
	less ./aligned_seqs/*failures.fasta >> LOG
	echo "" >> LOG

	# Filtering alignment

	echo "Filtering alignment..."
	filter_alignment.py -i ./aligned_seqs/*aligned.fasta -s -e $entropy_filter -g $align_gap_filter -o "filter_alignment_G${align_gap_filter}_E${entropy_filter}"

	# Make phylogenetic tree

	echo "Making phylogenetic tree..."
	make_phylogeny.py -i ./filter_alignment_G${align_gap_filter}_E${entropy_filter}/*.fasta -o tmp.tre -t $tree_building_method
	#remove "seq_" from trees generated with raxml
	sed 's/seq_//g' tmp.tre > makephylo_${tree_building_method}.tre
	rm tmp.tre

	# Make OTU table

	echo "Making OTU Table..."
	biom convert -i ./MATRIX-COUNT_transposed.txt -o OTU_Table.biom --to-json --table-type="OTU table"

	echo "####################" >> LOG
	echo "OTU Table Summary:" >> LOG
	echo "####################" >> LOG
	biom summarize-table -i OTU_Table.biom >> LOG
	echo "" >> LOG

	# Filter OTUs by removing singles and OTUs that have very few reads

	echo "Filtering OTUs..."
	filter_otus_from_otu_table.py -i OTU_Table.biom -o removed_singles_few_reads.biom --min_count $min_num_reads

	echo "" >> LOG
	echo "OTU Table Summary After Filtering:" >> LOG
	biom summarize-table -i removed_singles_few_reads.biom >> LOG
	echo "" >> LOG

	echo "Adding obs metadata..."

	biom add-metadata -i removed_singles_few_reads.biom -o OTU_Table_wtaxa.biom --observation-metadata-fp ./assign_taxonomy/*.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy

	biom summarize-table -i ./OTU_Table_wtaxa.biom | tee -a LOG

fi

# ==================================================

if [ ${doOTUPICK} == 'True' ]; then

	echo "pick_otus:otu_picking_method uclust_ref" > pickOTUParameters.txt
	echo "pick_otus:similarity 0.99" >> pickOTUParameters.txt
	
	echo "Picking OTUs..."
	pick_closed_reference_otus.py -i ./Trimmed_Quality_Filtered_MSL/seqs_trim_clip.fna -p pickOTUParameters.txt -r ${rep_set} -o PICKOTUS_uclust_ref_99 -s 

	
fi


echo "DONE!"
echo "DONE" >> LOG
