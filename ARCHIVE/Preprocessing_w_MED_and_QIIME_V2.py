#!/usr/bin/env python

## EDITED: Melissa-- to make it compatible with multiple things.
# Version 1:
# Adapted from original preprocessing script; re-done in python

# Input:
# Folder containing sequence files in .fastq or .fastq.gz format, path of mapping file, paths of database files, parameters for MED and QIIME to process sequence file
# Output:
# filtered OTU table with taxonomy added, representative sequences for each OTU, phylogenetic tree of sequences, various log files, other intermediate files from MED and QIIME


# Requirements: #links good as of May 2016
# QIIME or macqiime is installed http://qiime.org/
# fastx toolkit is installed http://hannonlab.cshl.edu/fastx_toolkit/
# MED is installed http://merenlab.org/2014/08/16/installing-the-oligotyping-pipeline/
# folder with raw fastq files must not contain any non-fastq files

#=========================================================

import argparse
import os, os.path
import sys
import subprocess

#=========================================================

parser = argparse.ArgumentParser(
	description="Pre-process data to be ready for analysis.")
parser.add_argument(
	'-o',
	'--project_name',
	help = "Project name; will be name of output folder",
	required = False,
	default = 'Preprocessing_Output')
parser.add_argument(
	'-f',
	'--fasta_files',
	help = 'Full filepath to folders with fasta files for pre-processing. If multiple, make sure they are comma-separated. Files can be raw seq from Dalhousie, or they can be demultiplexed reads',
	required = True)
parser.add_argument(
	'-D',
	'--Demultiplexed',
	help = "Set as 'True' if the input files are actually demultiplexed already and just need to be concatenated. Can be comma-separated, if multiple files.",
	required = True,)
parser.add_argument(
	'-P',
	'--PATH',
	help = 'Full path to extra shell scripts',
	required = True)
parser.add_argument(
	'--SEP_files',
	help = 'Raw fasta files can either be a single file (usually 16s) or individual files (usually 18s). If there are separate files for each sample, mark as True. If not, leave blank or mark as False. If included, should have one for each set of files.',
	required = False,
	default = 'False')
parser.add_argument(
	'-m',
	'--mapping_fps',
	help = 'Full filepath to mappingfile(s). MUST BE IN CORRESPONDING ORDER TO FASTA FILES, IF MULTIPLE. Mapping file MUST be QIIME-compatible and contain all desired samples. To ensure QIIME-compatibility, please see validate_mapping_file.py in QIIME',
	required = True)
parser.add_argument(
	'-r',
	'--rep_set',
	help = 'Full filepath to reference database rep set (Fasta file)',
	required = True)
parser.add_argument(
	'-t',
	'--reference_tree',
	help = 'Map of reference tree (Text file)',
	required = True)
parser.add_argument(
	'-a',
	'--aligned_rep_set',
	help = "Full filepath to aligned rep set (Fasta file)",
	required = True)
parser.add_argument(
	'-l',
	'--length_reads',
	help = "Trimming length for reads [Default: 250]",
	required = False,
	type = int,
	default = 250)
parser.add_argument(
	'-b',
	'--tree_building_method',
	help = "Tree Building Method. Choices: clustalw, raxml_v730, muscle, fasttree, clearcut. [Default: fasttree]",
	required = False,
	choices = ["clustalw","raxml_v730","muscle","fasttree","clearcut"],
	default = 'fasttree')
parser.add_argument(
	'-g',
	'--align_gap_filter',
	help = "Alignment Gap Filtering Parameter. A value from 0.0 to 1.0, representing the percentage of gaps that are required to filter out an alignment position. Please see filter_alignment.py documentation for more details. [Default: 0.9]",
	required = False,
	type = float,
	default = 0.9)
parser.add_argument(
	'-e',
	'--entropy_filter',
	help = "Alignment Entropy Filtering Paramtere. A value from 0.0 to 1.0, representing the N/1.0 percentage of positions to be removed based on entropy. Please see filter_alignment.py documentation for more details. [Default: 0.05]",
	required = False,
	type = float,
	default = 0.05)
parser.add_argument(
	'-n',
	'--min_num_reads',
	help = "Minimum number of reads to include in final OTU table [Default: 2]",
	required = False,
	type = int,
	default = 2)
parser.add_argument(
	'-M',
	'--MED_value',
	help = "MED value; if given, program will jump directly into calculating MEDs. If not given, then it will prompt you to give MED value later. [Default: None]",
	required = False,
	default = "None")

#=========================================================
# SET ARGS

args = parser.parse_args()

project_name = args.project_name
fasta_files = args.fasta_files.split(",")
Demultiplexed = args.Demultiplexed.split(",")
PATH = args.PATH
SEP_files = args.SEP_files.split(",")
mapping_fps = args.mapping_fps.split(",")
rep_set = args.rep_set
reference_tree = args.reference_tree
aligned_rep_set = args.aligned_rep_set
length_reads = args.length_reads
tree_building_method = args.tree_building_method
align_gap_filter = args.align_gap_filter
entropy_filter = args.entropy_filter
min_num_reads = args.min_num_reads
MED_value = args.MED_value


#==============
# To be deleted; checking
#TBD

# print project_name
# print fasta_files
# print Demultiplexed
# print PATH
# print SEP_files
# print mapping_fps
# print rep_set
# print reference_tree
# print aligned_rep_set
# print length_reads
# print tree_building_method
# print align_gap_filter
# print min_num_reads
# print MED_value

#=========================================================
# CHECK TO MAKE SURE ALL FILES MAKE SENSE

# Project Name can be any string
# Fasta Files must be folder

n = 0
for i in fasta_files:
	if os.path.isfile(i) and Demultiplexed[n]:
		print "File and Demultiplexed status matches"
		n += 1
	elif os.path.isdir(i) and Demultiplexed[n]:
		print "File and Demultiplexed status matches"
		n += 1
	else:
		sys.exit("ERROR: Raw fasta file path is not a directory OR demultiplexed file paths are not files.")

# Demultiplexed is dealt with above
# Make sure that there are same number of fasta files, mapping files, and 'types'
if len(SEP_files) == 1 and SEP_files[0] == 'False' and len(fasta_files) != 1: # If you put nothing under SEP, it will default put 'False'-- just checking to make sure you have the appropriate number of fasta files
	if len(fasta_files) == len(mapping_fps):
		pass
	else:
		sys.exit("ERROR: Number of fasta files and mapping files to not match")
else:
	if len(fasta_files) == len(SEP_files):
		pass
	else:
		sys.exit("ERROR: Number of SEP_files does not match number of fasta files")
# Mapping files must be valid, but we do not check for this.
# Rep set, ref tree, aligned rep are not checked
# Length reads should be number
# Length, tree building method, gap, entropy, min_num_reads and MED value all not checked

#=========================================================
# WRITE VARIABLES	
os.system("mkdir " + project_name)
os.system("cd " + project_name)

var = open("VAR", 'a')
#---
var.write('project_name=' + str(project_name) + "\n")
#---
n = 0
fasta_files_string = ''
for i in range(len(fasta_files)):
	if n == len(fasta_files):
		fasta_files_string += str(i)
	else: 
		fasta_files_string += str(i) + ','
	n += 1
var.write('fasta_files=' + str(fasta_files_string) + "\n")
#---
Demultiplexed_string = ''
n = 0
for i in range(len(Demultiplexed)):
	if n == len(Demultiplexed):
		Demultiplexed_string += str(i)
	else: 
		Demultiplexed_string += str(i) + ','
	n += 1
var.write('Demultiplexed=' + str(Demultiplexed_string) + "\n")
#---
var.write('PATHfiles=' + str(PATH) + "\n")
#---
mapping_fps_string = ''
n = 0
for i in range(len(mapping_fps)):
	if n == len(mapping_fps):
		mapping_fps_string += str(i)
	else: 
		mapping_fps_string += str(i) + ','
	n += 1
var.write('mapping_fps=' + str(mapping_fps_string) + "\n")
#---
var.write('rep_set=' + str(rep_set) + "\n")
#---
var.write('reference_tree=' + str(reference_tree) + "\n")
#---
var.write('aligned_rep_set=' + str(aligned_rep_set) + "\n")
#---
var.write('length_reads=' + str(length_reads) + "\n")
#---
var.write('tree_building_method=' + str(tree_building_method) + "\n")
#---
var.write('align_gap_filter=' + str(align_gap_filter) + "\n")
#---
var.write('entropy_filter=' + str(entropy_filter) + "\n")
#---
var.write('min_num_reads=' + str(min_num_reads) + "\n")
#---
var.write('MED_value=' + str(MED_value) + "\n")
#---
var.close()

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

#=========================================================

#=========================================================
# FIGURE OUT WHAT TYPE OF DATA INPUT
os.chdir("./"+project_name
os.system("mkdir FASTAFILES")
os.system("mkdir MAPPINGFILES")

os.system("touch mappingfilenames.txt") # for storying mapping file names; used to call names later.

pos = 0
for i in fasta_files:
	TF_demulti = Demultiplexed[pos]
	if len(mapping_fps) == 1:
		mapFPS = mapping_fps[0]
	else:
		mapFPS = mapping_fps[pos]
	separateFiles = SEP_files[pos]
	open("map_fps.txt", 'w').write('mappingfile='+str(mapFPS)) # write mapping file; will be over-written each loop
	open("rep.txt", 'w').write('rep='+str(pos+1))
	open("sep.txt", 'w').write('sep='+str(separateFiles))
	print TF_demulti
	if TF_demulti == 'False': # Raw data
		open("fasta_fps.txt", 'w').write('fasta='+str(i)) # write fp for fasta file; will be over-written each loop
		pwdRawData = PATH.strip()+'/RawData.sh'
		print 'Executing RawData script' #TBD
		subprocess.call([str(pwdRawData)])
	else: # Demultiplexed already
		open("fasta_fps.txt", 'w').write('fasta='+str(i)) # write fp for fasta file; will be over-written each loop
		pwdDemulti = PATH.strip()+'/DemultiplexedData.sh'
		print 'Executing demultiplexed script' #TBD
		subprocess.call([str(pwdDemulti)])
	pos += 1

totalFiles = pos # Save the total number of fasta files we have


#=========================================================
# Now, we should have several folders that have the demultiplexed data and new metadata files
# Next script will concatenate them and do some trimming/filtering etc

# Make a file that is a comma separated list of mapping file names
mappingfilenames = open("mappingfilenames.txt",'r').read().strip()
mappingfilenames = mappingfilenames.replace("\n",",")
open("mappingfilenames.txt",'w').write('mappingfilenames='+mappingfilenames)

# Call shell script that will do the rest

pwdPartTwo = PATH.strip()+'/Processing_bash.sh'
subprocess.call([str(pwdPartTwo)])

#=========================================================

print "DONE"

