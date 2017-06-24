#!/bin/bash

# Input:
# Demultiplexed 18s data from Halifax
# All fastq files should be in a single folder

# Requirements:
# macqiime is installed 
# fastx toolkit is installed
# MED is installed




#------------------------


# SECTION ONE: Filepaths and parameters
# This section prompts the user to specific a bunch of filepaths and parameters
# It also records their specifications in the log so they know exactly what they did
# Last thing it does is count the seqs in each file prior to processing.


#------------------------
# Making Project directory. script will prompt user to type in a folder name and press enter. Then, we move inside the folder.
# All files will be created inside this directory. This means I don't need complete filepaths in order to do all commands.
# First, there will be a series of questions to ask. These answers will specify the parameters you want.

echo "Please enter project name"
	read projectname
	mkdir "$projectname"
	cd "$projectname"
	
# Now, we are in the project directory.

# Copying raw data from folder and putting it in $projectname folder.
# I wanted to copy it because then you can't mess up your raw data. Just in case!
# cat is used so I'm not tampering with original file.

echo "Enter complete pathway to raw fastq files. Files must be in a single folder."
	read fastqlocation
	
# You should end up with a folder called "raw_fastq" inside your project directory.
# Next, we will save other complete file paths for the database as variables to be used later.
# I used the SILVA111 database for this, but other databases can be used as well.
# I'm also copying these variables into the LOG file so we know what we used in the past if we need to go back.

	
echo "Enter complete pathway to database rep set (fasta)"
	read reftreefasta
echo "Enter complete pathway to map of ref tree (txt)"
	read reftreemap
echo "Enter complete pathway to aligned rep set (fasta)"
	read reftreealigned
	
echo " " >> LOG
echo "Input information: complete pathways to reference database" >> LOG
echo " " >> LOG

echo "Reftreefasta

$reftreefasta


Reftreemap

$reftreemap


Reftreealigned

$reftreealigned
" >> LOG

# Now, we enter some parameters and variables. There is a prompt first, and then the user may type in their conplete filepath.


# Enter quality threshold. This will be used in fastq_quality_filter. All sequences below this filter level will be discarded.

echo "Enter quality threshold (usually 20) for filtering fasta files"
	read qthreshold
	
# Enter the trimming length for fastx_clipper and fastx_trimmer.
echo "Enter trimming length for reads"
	read trimlength

# Enter Tree building method. Fasttree is fast and the default; raxml seems to be more popular. Although, I've read somewhere they're about the same?
echo "Enter tree building method. Method for tree building. Valid choices are: clustalw,raxml_v730, muscle, fasttree, clearcut"
	read treemethod
	
# Enter the minimum count fraction. This is for filtering the OTU table.
echo "Enter minimum count fraction to include in final OTU table"
	read minimumcountfraction
	
echo "
Quality threshold: $qthreshold
trimlength: $trimlength
Tree method: $treemethod
Minimumcountfraction: $minimumcountfraction
" >> LOG
	
	
echo "Copying raw fast files..."
		cp -r "$fastqlocation" ./raw_fastq
		
echo "Doing preamble stuff..."

# Printing preamble for LOG. The LOG file will contain intermediate summaries for use of troubleshooting.
echo " " >> LOG
echo "Log and script summaries for $projectname" >> LOG
echo " " >> LOG
	
# Note that all other parameters are default.


# Counting sequences for each sample prior to processing. 
# First, we go into the raw-fastq file. Then, for 'file' that ends in fastq, 
# Print the name of the file and then print the count of ">" in that file (aka number of sequence of reads in file)


echo "RAW SEQUENCE COUNTS" >> LOG
cd raw_fastq
for f in *.fastq; do
	echo "$f" | tee -a ../LOG | grep -c ">" "$f" >> ../LOG
done
echo " " >> LOG

cd ..

# We are now in $projectname directory again


#------------------------

# SECTION TWO: Pre-MED quality filtering and trimming
# This section does the following:
# Quality filters all sequences according to $qthreshold input
# Combines seq files into a single fasta file with proper names (multiple split libraries)
# Trims and clips all sequences to $trimlength input
# Finally, it counts all sequences again and records it in the LOG
# At the very end, it alters the format of the seq.fasta files to the MED-friendly format


#------------------------



# Quality filtering of raw fastq files using parameters of $qthreshold.
# This section goes through all files separately and filters through quality threshold. 
# The reason we do this separately is because quality_filtered_fastq only operates on fastq files.
# After combining all files with multiple split libraries, the format becomes a regular FNA file and quality_filtered_fastq doesn't like that.


echo "Starting Quality Filtering..."

mkdir Quality_Filtered_Fastq

cd raw_fastq
for f in *.fastq; do
	fastq_quality_filter -i "$f" -q "$qthreshold" -p 99 -o ../Quality_Filtered_Fastq/"$f"
done

cd ..
# in $projectname directory again

# Now, recounting number of sequences after quality filtering

echo " " >> LOG
echo "POST-QUALITY FILTER SEQUENCE COUNT" >> LOG
echo "q=$qthreshold" >> LOG
echo " " >> LOG

# This prints the filename of each sample, then uses that name to count the reads in each file.

cd Quality_Filtered_Fastq
for f in *.fastq; do
	echo "$f" | tee -a ../LOG | grep -c ">" "$f" >> ../LOG
done
cd ..
# in $projectname directory again

echo "Quality Filtering Complete"

#------------------------

# Run multiple_split_libraries.py using default settings
# split_libraries default is:
# 	max_bad_run_length:3
# 	min_per_read_length_fraction:0.75
#	sequence_max_n:0
#	phred_quality_threshold:3
# etc
# The sample ID indicator means that the stuff before this part will act as the new ID name.


multiple_split_libraries_fastq.py -i ./Quality_Filtered_Fastq -o multiple_split_libraries_fastq --sampleid_indicator ".fastq"

# Output should yield histograms.txt; log; seqs.fna; split_library_log.txt

echo "Multiple_split_libraries complete."

# Still in 'home' directory

# Trimming files to $trimlength bp using fastx.

mkdir Trimmed_Quality_Filtered_MSL

cd multiple_split_libraries_fastq

# Using the $trimlength inputted before, use fastx trimmer to trip the seq file from multiple split libraries.

fastx_trimmer -i seqs.fna -f 1 -l "$trimlength" -o ../Trimmed_Quality_Filtered_MSL/seqs_trim.fna

cd ..
# in $projectname directory again

#------------------------

# Counting sequences after trimming

echo " " >> LOG
echo "POST-TRIMMING SEQUENCE COUNT" >> LOG
echo " " >> LOG

cd Trimmed_Quality_Filtered_MSL

echo "seqs_trim.fna" >> LOG
grep -c ">" seqs_trim.fna >> ../LOG

echo "Trimming Complete"

cd ..
# Back in home directory

# Filter out sequences with less than $trimlength. Still in Trimmed_Quality_Filtered_MSL

fastx_clipper -i ./Trimmed_Quality_Filtered_MSL/seqs_trim.fna -l "$trimlength" -o ./Trimmed_Quality_Filtered_MSL/seqs_trim_clip.fna

# Counting sequences after clipping

# Writing title

echo " " >> LOG
echo "POST-CLIPPING SEQUENCE COUNT" >> LOG
echo " " >> LOG

# Counting sequences
echo "seqs_trim.fna" >> LOG
grep -c ">" ./Trimmed_Quality_Filtered_MSL/seqs_trim_clip.fna | tee -a ../LOG

# in $projectname directory again

echo "Clipping Complete"

#------------------------

# In home directory.
# Now, we want to alter the seq.fna file so that the format fits with MED requirements.
# MED likes it to be SampleXXX_ReadXXX, so this is what we'll do.
# Note that this ONLY works for halifax data so far-- you'll likely have to change this for other formats. 
# The halifax data has the reads look like this:
# >E-s-PM-1-Jul-a_S265_L001_R1_001_7 M02352:69:000000000-AKVT6:1:1101:21292:5274 1:N:0:265 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
# CCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGTGGTTAAAAAGCTCGTAGTTGGATCTCAACAGCCTTGAAGCGGTTAACTTTGTAGTTTTACTGCTTTATAGGCTGTGCTTTTGCTGGAGGTTTAGGGGTGCTCTTGATCGAGTGTCTCTGGATGCCAGCAGGTTTACTTTGAAAAAATTAGAGTGCTCAAAGCAGGCTAAAATGCCTGAATATTCGTGCATGGAATAATAGAATAGGA
# So, I want to get rid of everything that is between "_S" and "_R" (inclusive) and replace it all with "_Read"
# Then, the next line gets rid of everything after _001 (although I'm not sure this works for ALL halifax data???


cd ./Trimmed_Quality_Filtered_MSL

sed 's/\_\S.*\_\R/\_\Read/g' seqs_trim_clip.fna > seqs_MED1.fna
# The above gets rid of all the spaces and replaces R with Read

sed 's/_001.*//g' seqs_MED1.fna > seqs_MED.fna
# This gets rid of all the post-junk. This part may not work with other formats.

# Putting output in its own MED file
rm seqs_MED1.fna
mkdir MED
mv seqs_MED.fna ./MED/seqs_MED.fna

cd ..
# In the home directory


#------------------------

# SECTION THREE: MED decompose and transposing
# This section decomposes the fasta files using the $minimumentropy input
# Also, transposes the resulting MATRIX-COUNT.txt file into a QIIME-friendly format
# The last part changes the format of NOE-REPRESENTATIVES.fasta to be compatible with QIIME
# (explained below)



#------------------------

	
# Enter the minimum entropy. This is the only variable in this script that changes for MED (Minimum Entropy Decomposition)
echo "Enter Minimum Entropy (MED)-- the default is total reads divided by 5000."
	read minimumentropy
	
echo "

minimum entropy: $minimumentropy

"

# Now, decompose fasta file. This is a single line.

echo "Decomposing fasta file..."

decompose ./Trimmed_Quality_Filtered_MSL/MED/seqs_MED.fna -M "$minimumentropy" --gen-html -o decompose


# After decomposition, we need to transpose the MATRIX-COUNT.txt file because it is in the wrong format.
# I was too lazy to write my own script for transposing so I copied one from stackoverflow and tested it before including it below
# This part of the script was copy and pasted from stackoverflow

echo "Transposing..."
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
            str=str" "a[i,j];
        }
        print str
    }
}' ./decompose/MATRIX-COUNT.txt > MATRIX-COUNT_trans.txt

# Now, we need to replace the spaces with tabs, since the transposing above puts spaces instead of tabs.
echo "Replacing spaces with tabs"
awk -v OFS="\t" '$1=$1' MATRIX-COUNT_trans.txt > MATRIX-COUNT_transposed_temp.txt
echo "Deleting Matrix-count_trans.txt"
rm MATRIX-COUNT_trans.txt	# Deletes intermediate file


# A few more alterations-- changes the top left name from samples to #OTU ID

sed 's/samples/#OTU ID/g' MATRIX-COUNT_transposed_temp.txt > MATRIX-COUNT_transposed.txt

rm MATRIX-COUNT_transposed_temp.txt 	# Deletes intermediate file

# Also, remove the pesky "|size:XXX" at the end of every OTU ID.
# This is necessary because of two reasons:
# 1. Adding metadata requires that the observation_metadata doesn't have the "|size:XXX"
# 2. The tree cannot have "|size:XXX" otherwise the tree tips will not match up.

sed 's/|size:[0-9]*$//g' ./decompose/NODE-REPRESENTATIVES.fasta > NODE_REP_FORDOWNSTREAM.fasta


#------------------------

# SECTION FOUR: Making a tree and OTU Table
# This last section goes back to QIIME to make a phylotree and the OTU table
# You can specify which tree you want
# It also records some of the intermediate results in the LOG
# Finally, it adds observation metadata to the file.
# The resulting OTU table and tree have been tested to be compatible for: alpha_rarefaction,beta_diversity_through_plots,summarize_taxa_through_plots 


#------------------------


# Assign Taxonomy using previous filepaths of references database. Need to use the unaligned database.

echo "Assigning taxonomy..."
assign_taxonomy.py -i NODE_REP_FORDOWNSTREAM.fasta -o assign_taxonomy -r "$reftreefasta" -t "$reftreemap"

# Align sequences. This step uses the aligned reference database you supplied.

echo "Aligning sequences..."
align_seqs.py -i NODE_REP_FORDOWNSTREAM.fasta -t "$reftreealigned" -o aligned_seqs


# Now, writing the newly aligned sequences and failures to the log for easy access later.

echo " "
echo "Aligned Sequences" >> LOG
echo " "
less ./aligned_seqs/*aligned.fasta >> LOG
echo " "
echo "Aligned failures" >> LOG
echo " "
less ./aligned_seqs/*failures.fasta >> LOG
echo " "


# Filtering alignment

echo "Filtering alignment..."
filter_alignment.py -i ./aligned_seqs/*aligned.fasta -s -e 0.10 -o filter_alignment

# Make phylogeny

echo "Making phylogeny..."
make_phylogeny.py -i ./filter_alignment/*.fasta -o 18s_makephylo_fastree.tre -t "$treemethod"


# Make OTU table

echo "Making OTU Table..."
biom convert -i ./MATRIX-COUNT_transposed.txt -o 18s_OTU_Table.biom --to-json --table-type="OTU table"


echo "OTU Table Summary" >> LOG
echo " " >> LOG
biom summarize-table -i 18s_OTU_Table.biom >> LOG

# Filter OTUs by removing singles and OTUs that have very few reads

echo "Filtering OTUs..."
filter_otus_from_otu_table.py -i 18s_OTU_Table.biom -o removed_singles_few_reads.biom --min_count_fraction "$minimumcountfraction"

echo " " >> LOG
echo "OTU Table Summary After Filtering" >> LOG
echo " " >> LOG

biom summarize-table -i removed_singles_few_reads.biom >> LOG

echo "Adding obs metadata..."

biom add-metadata -i removed_singles_few_reads.biom -o 18s_OTU_Table_wtaxa.biom --observation-metadata-fp ./assign_taxonomy/*.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy

biom summarize-table -i ./18s_OTU_Table_wtaxa.biom | tee -a LOG

echo "DONE"
