#!/bin/bash

# Deals with multiple folders of mapping files and others.
# Trimming files to $trimlength bp using fastx.
# ==================================================
# Sources and summary of contents

source VAR
	# print project_name
	# print fasta_files
	# print Demultiplexed
	# print PATH
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
	
# ==================================================
# Setting variables; making fasta files and mapping files

# Concatenate the all fna files to make a single file for trimming and clipping
cat ./FASTAFILES/*.fna > ./FASTAFILES/ALL_seqs.fna

# Merge the mapping files using a QIIME script
merge_mapping_files.py -m $mappingfilenames -o Merged_mapping_file.txt

validate_mapping_file.py -m Merged_mapping_file.txt -o validate_mapping_file

print "DONE!"

