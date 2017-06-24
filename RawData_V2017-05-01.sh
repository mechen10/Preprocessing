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
source fasta_fps.txt
	# $fasta = fasta filepath
source map_fps.txt
	# $mappingfile = mapping filepath
source rep.txt
	# $rep = int
source sep.txt
	# $sep = 'True' or 'False'

#======================================
# Run multiple_split_libraries.py using default settings
# split_libraries default is:
# 	max_bad_run_length:3
# 	min_per_read_length_fraction: 0.75
#	sequence_max_n:0
# etc
# We let quality filter threshold
# The sample ID indicator means that the stuff before this part will act as the new ID name.
echo "Running split_libraries QIIME script..."
if [ $sep == True ] 
then
	# touch parameters_mult_split_libraries.txt
	# echo "split_libraries_fastq:phred_quality_threshold 19" >> parameters_mult_split_libraries.txt
    multiple_split_libraries_fastq.py -i $fasta -o split_libraries_fastq$rep --demultiplexing_method sampleid_by_file # -p parameters_mult_split_libraries.txt
    # Output should yield histograms.txt; log; seqs.fna; split_library_log.txt
    echo "multiple_split_libraries script finished...for rep $rep"
    echo "Copying fna file"
    cp ./split_libraries_fastq$rep/seqs.fna ./FASTAFILES/seqs$rep.fna
    echo "Copying mappingfile"
    cp $mappingfile ./MAPPINGFILES/mappingfile$rep.txt
    echo "./MAPPINGFILES/mappingfile$rep.txt" >> mappingfilenames.txt
fi
###############
if [ $sep == False ]
then
	# I have removed quality filtering, as per Matt's recommendation. Will do this AFTER trimming/clipping
	echo $fasta
	ourindex=`ls $fasta | grep -v R1 | grep -v R2 | grep -v .gz | grep -v .txt` #get the name of the index file
	echo $ourindex
	ourRun=`ls $fasta | grep R1 | grep -v gz`
	echo $ourRun
    echo "preparing sequences with QIIME..."	
	split_libraries_fastq.py -i $fasta$ourRun -o ./split_libraries_fastq$rep --barcode_read_fps $fasta$ourindex --mapping_fps $mappingfile --barcode_type 12 # --phred_quality_threshold 19 
    echo "split_libraries_fastq script finished rep $rep..."
    echo "Copying fna file"
    cp ./split_libraries_fastq$rep/seqs.fna ./FASTAFILES/seqs$rep.fna
	echo "Copying mappingfile"
    cp $mappingfile ./MAPPINGFILES/mappingfile$rep.txt
    echo "./MAPPINGFILES/mappingfile$rep.txt" >> mappingfilenames.txt
fi
echo "Done copying non-Demultiplexed files"