#!/bin/bash/environ python



# This is the master script controlling all the workflows post-preprocessing.
# Here, you can change the pathways of the inputs to other workflows. Alternatively you can make manual source files that will be sourced from other workflows. 
# Note that BEFORE running this script, you should already know which OTUs you want to remove due to contaminations. Do this by converting the biom table into text format and comparing negatives
# (PS I should make a script to do this later)


# Make choice of running different things



WORKFLOWFOLDER=
RunRScripts=


doFilter=true
doAlphaBetaTaxa=true
doRscripts=true

############### VARIABLE FP ###############
# Raw OTU table from pre-processing script
OTUTABLE=
# Original merged MF; output form pre-processing script
MF=
# Minimum reads in sample to keep sample
minthresh=
# Minimum count of OTU in sample to keep OTU
mincount=
# Mapping file column:type to get rid of (usually negatives)
getrid=
# Filepath to OTU contaminants to remove; list in txt format
OTUCONTAMFP=
# Phylogenetic tree
treeFP=
# Max rarefying depth when doing alpha_rarefaction
maxE=
# Rarefying depth
rareValue=
# Alpha and Beta diversity metrics
alphametrics=chao1,PD_whole_tree,observed_otus
betametrics=bray_curtis,unweighted_unifrax,weighted_unifrac



############### FILTERING MITO, CHLP, EUK, CONTAM ###############


# Get rid of CHLP, MITO, CONTAM, EUK

if $doFilter ; then
	
	mkdir OTU_MP_filt
	mkdir OTU_MP_filt/intermediate_files
	
	cp $OTUCONTAMFP ./OTU_MP_filt/intermediate_files/TODELETE_contaminants.txt

	biom convert -i $OTUTABLE -o OTUTABLE_original.txt --to-tsv --header-key taxonomy

	grep -v "Chloroplast" ./OTUTABLE_original.txt | grep -v "Mitochondria" | grep -v "Eukaryota" >> nochlpmito.txt

	filter_otus_from_otu_table.py -i $OTUTABLE -e ./nochlpmito.txt --negate_ids_to_exclude -o TEMPBIOM.biom
	filter_otus_from_otu_table.py -i TEMPBIOM.biom -e ./OTU_MP_filt/intermediate_files/TODELETE_contaminants.txt -n $mincount -o TEMPBIOM2.biom
	filter_samples_from_otu_table.py -i TEMPBIOM2.biom -o ./OTU_MP_filt/OTU_Table_nochlpmito_m${minthresh}.biom -m $MF -s $getrid -n $minthresh --output_mapping_fp ./OTU_MP_filt/MF_nochlpmito_m${minthresh}.txt

	rm TEMPBIOM.biom
	rm TEMPBIOM2.biom
	
	mv OTUTABLE_original.txt ./OTU_MP_filt/intermediate_files/OTUTABLE_original.txt
	mv nochlpmito.txt ./OTU_MP_filt/intermediate_files/nochlpmito.txt

	echo "DONE filtering chlp, mito, contam, euk"

fi

############### ALPHA BETA TAXA ####################

if $doAlphaBetaTaxa ; then

	mkdir ANALYSIS_ALPHABETATAXA
	cd ANALYSIS_ALPHABETATAXA

	mkdir OTU_Tables_and_MP


	########
	# ALPHA

	echo "alpha_diversity:metrics "$alphametrics > alphaparams.txt

	alpha_rarefaction.py -i ../OTU_MP_filt/OTU_Table_nochlpmito_m${minthresh}.biom -m ../OTU_MP_filt/MF_nochlpmito_m${minthresh}.txt -p alphaparams.txt -t $treeFP -e $maxE -o alpha_rare

	multiple_rarefactions_even_depth.py -i ../OTU_MP_filt/OTU_Table_nochlpmito_m${minthresh}.biom -d $rareValue -o ./OTU_Tables_and_MP/Rarefied_Tables_r${rareValue}/

	alpha_diversity.py -i ./OTU_Tables_and_MP/Rarefied_Tables_r${rareValue}/ -m $alphametrics -t $treeFP -o alpha_div/

	collate_alpha.py -i alpha_div/ -o collated_alpha/
	
	# Get filepaths of files
	echo "tocollate=" > tocollateFP.txt
	first=true
	for F in ./collated_alpha/*; do
		if $first ; then
			echo $F >> tocollateFP.txt
			first=false
		else
			echo ,$F >> tocollateFP.txt
		fi
	done
	tr -d '\n' < tocollateFP.txt > tocollateFP.txt
	
	source tocollateFP.txt
	rm tocollateFP.txt
	
	add_alpha_to_mapping_file.py -i $tocollate -m ../OTU_MP_filt/MF_nochlpmito_m${minthresh}.txt --depth $rareValue --collated_input -o ./OTU_Tables_and_MP/MF_withalpha.txt


	########
	# BETA

	echo "beta_diversity:metrics "$betametrics > betaparams.txt

	beta_diversity_through_plots.py -i ../OTU_MP_filt/OTU_Table_nochlpmito_m${minthresh}.biom -m ../OTU_MP_filt/MF_nochlpmito_m${minthresh}.txt -o beta_div -e $rareValue -p betaparams.txt -t $treeFP

	jackknifed_beta_diversity.py -i ../OTU_MP_filt/OTU_Table_nochlpmito_m${minthresh}.biom -m ../OTU_MP_filt/MF_nochlpmito_m${minthresh}.txt -o jackknifed_beta -e $rareValue -t $treeFP -p betaparams.txt 


	########
	# TAXASUMMARIES

	summarize_taxa_through_plots.py -i ../OTU_MP_filt/OTU_Table_nochlpmito_m${minthresh}.biom -m ./OTU_MP_filt/MF_nochlpmito_m${minthresh}.txt -o summarize_taxa -s	
	
fi

	############### RUN ADDITIONAL R SCRIPTS ####################

	if $doRscripts ; then

		$RunRScripts

	fi

