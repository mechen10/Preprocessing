#!/bin/bash
#===============================================================
# PART ONE: DEMULTIPLEXING
# First, load in information from py. 
# source VAR
	# print project_name
	# print fasta_files
	# print Demultiplexed
	# print PATH
	# print mapping_fps
	# print rep_set
	# print reference_tree
	# print aligned_rep_set
	# print length_reads
	# print tree_building_method
	# print align_gap_filter
	# print entropy_filter
	# print min_num_reads
	# print MED_value
#===============================================================

source fasta_fps.txt
	# $fasta = fasta filepath
source map_fps.txt
	# $mappingfile = mapping filepath
source rep.txt
	# $rep = int
source sep.txt
	# $sep = 'True' or 'False'
echo "Copying files"
cp $fasta ./FASTAFILES/seqs$rep.fna
cp $mappingfile ./MAPPINGFILES/mappingfile$rep.txt
echo "./MAPPINGFILES/mappingfile$rep.txt" >> mappingfilenames.txt
echo "Done copying demultiplexed $rep files"
