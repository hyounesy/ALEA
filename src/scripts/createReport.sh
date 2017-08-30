#!/bin/bash

pushd `dirname $0` > /dev/null
AL_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

# MANUAL INSTALLATION
#source $AL_DIR_TOOLS/alea.config
# DOCKER INSTALLATION
source /alea-data/alea.config

##############################################################################
#############   Module 5: statistical analysis
##############################################################################

PARAM_VALID=1
PARAM_SINGLE_READS=1
if [ $AL_USE_CONCATENATED_GENOME = 1 ]; then
    if [ $# -eq 4 ]; then
        PARAM_BAM_PREFIX=$1
	PARAM_STRAIN1=$2
	PARAM_STRAIN2=$3
        PARAM_INTERVALS=$4
    else
        PARAM_VALID=0
    fi
else
    PARAM_VALID=0
fi
VAR_OUTPUT_BASENAME=$(basename $PARAM_BAM_PREFIX)
VAR_OUTPUT_DIR=$(dirname $PARAM_BAM_PREFIX)
VAR_OUTPUT_PREFIX1="$VAR_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"
VAR_OUTPUT_PREFIX2="$VAR_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"
VAR_OUTPUT_PREFIX3="$VAR_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_total_F"$AL_BAM2WIG_PARAM_FILTERING_FLAG"_q"$AL_BAM2WIG_PARAM_MIN_QUALITY"

if [ $PARAM_VALID = 0 ]; then
    echo "
Usage:
    alea createReport <alignReadsDir/bedGraphPrefix> <strain1> <strain2> <interval file>
Options:
    bedGraphPrefix       bedGraphPrefix used in createTracks
		    can be path to file (e.g. alignReadsDir/bam_prefix)
		    do not include _total.bedGraph in bedGraphPrefix

    strain1	    name of strain 1 (e.g. hap1)
    strain2	    name of strain 2 (e.g. hap2)

    intervals       General Feature Format (GFF) used for counting reads
		    can be user-define or chosen from the list in the folder alea-data/gff/
		    if using custom gff, please match supplied format
Output:

    analysis.tsv    a table with allelic coverage, total RPKM, parental ratios
		    and other statistics for each interval in the GFF file.

Examples:
    (AL_USE_BWA=1)
    alea createReport H3K4me3Liver C57BL6J CAST_EiJ ./gff/transcriptionStartSites_mm10.gff
    (AL_USE_STAR=1)
    alea createReport F1hybridLiver C57BL6J CAST_EiJ ./gff/RefseqGenes_mm10.gff

"
exit 1
fi



# choose from count, mean, max, min, etc.
OPERATION="mean"

function groomToBed5 {
	local INPUT_FILE=$1
	local FILE=$(basename "$INPUT_FILE")
	local FILE="${FILE%.*}"

	echo "converting $INPUT_FILE --> $FILE.bed5"
	if [[ $INPUT_FILE == *.*raph ]] ; then
		awk -v name="$FILE" '{OFS="\t";FS="\t"}{print $1, $2, $3, name, $4}' "$INPUT_FILE" > ${INPUT_FILE//.bedGraph/_unsort.bed5}

	elif [[ $INPUT_FILE == *.bed ]] ; then
		awk -v name="$FILE" '{OFS="\t";FS="\t"}{print $1, $2, $3, name, 1}' "$INPUT_FILE" > ${INPUT_FILE//.bed/_unsort.bed5}

	elif [[ $INPUT_FILE == *.bw ]] ; then
		bigWigToBedGraph "$INPUT_FILE" ${INPUT_FILE//.bw/.bedGraph}
		awk -v name="$FILE" '{OFS="\t";FS="\t"}{print $1, $2, $3, name, $4}' ${INPUT_FILE//.bw/.bedGraph} > ${INPUT_FILE//.bw/_unsort.bed5}
		rm ${INPUT_FILE//.bw/.bedGraph}

	elif [[ $INPUT_FILE == *.bed5 ]] ; then
		echo "file $INPUT_FILE already in bed5 format, sorting just in-case!"
		mv $INPUT_FILE ${INPUT_FILE//.bed5/_unsort.bed5}

	else
		echo "file $INPUT_FILE format not recognized. Please submit .bedGraph, .bed, .bw or .bed5 files"
	fi

	sort-bed "${INPUT_FILE%.*}"_unsort.bed5 > "${INPUT_FILE%.*}".bed5
	rm "${INPUT_FILE%.*}"_unsort.bed5
	echo "the bed5 table should have columns:"
	echo "#chr	start	end	name	score"
	echo "you file has columns:"
	head -n 2  "${INPUT_FILE%.*}".bed5
}

function bedmap_coverage {
	# THIS ADDS A COLUMN TO YOUR COORDINATES FILE. DUPLICATE IT TO PRESERVE ORIGINAL
	printProgress "[Calculating coverage] $COVERAGE over $COORDINATES"
	local COORDINATES=$1
	local COVERAGE=$2
	# keep coverage file name (column 4) for data integration
	local COLUMN_FILENAME=$(head -n 10 $COVERAGE | cut -f4 | tail -n 1)
	local COLUMN_NAME="$COLUMN_FILENAME"_"$OPERATION"
	echo "$COLUMN_NAME" > tmp

	# calculate coverage
	echo "calculating coverage over $COORDINATES using $COVERAGE and $OPERATION"
	bedmap --header --delim "\t" --$OPERATION $COORDINATES $COVERAGE >> tmp
	#for --mean & sum  (use COUNT), there are no values "0", instead you must sed NAN -> 0 #

	# add column with coverage values
	paste $COORDINATES tmp > tmp2
	sed 's/NAN/NaN/g' tmp2 > $COORDINATES
	rm tmp tmp2
}

# preserve original user-input interval file
DATE=$(date '+%y-%m-%d')
cp "$PARAM_INTERVALS" "$PARAM_INTERVALS"_"$DATE"

#ChIP
if [ $AL_USE_BWA = 1 ] || [ $AL_USE_BOWTIE2 = 1 ]; then
	groomToBed5 $VAR_OUTPUT_PREFIX1.bedGraph
	groomToBed5 $VAR_OUTPUT_PREFIX2.bedGraph
	groomToBed5 $VAR_OUTPUT_PREFIX3.bedGraph
	bedmap_coverage "$PARAM_INTERVALS"_"$DATE" $VAR_OUTPUT_PREFIX1.bed5
	bedmap_coverage "$PARAM_INTERVALS"_"$DATE" $VAR_OUTPUT_PREFIX2.bed5
	bedmap_coverage "$PARAM_INTERVALS"_"$DATE" $VAR_OUTPUT_PREFIX3.bed5
	printProgress "[Calculating coverage] succesfully finished"

#RNA
elif [ $AL_USE_STAR = 1 ] || [ $AL_USE_TOPHAT2 = 1 ]; then
        groomToBed5 $VAR_OUTPUT_PREFIX1.bedGraph
        groomToBed5 $VAR_OUTPUT_PREFIX2.bedGraph
        groomToBed5 $VAR_OUTPUT_PREFIX3.bedGraph
        bedmap_coverage "$PARAM_INTERVALS"_"$DATE" $VAR_OUTPUT_PREFIX1.bed5
        bedmap_coverage "$PARAM_INTERVALS"_"$DATE" $VAR_OUTPUT_PREFIX2.bed5
        bedmap_coverage "$PARAM_INTERVALS"_"$DATE" $VAR_OUTPUT_PREFIX3.bed5
	printProgress "[Merging exons into transcripts] started"
	sed -i 's/NaN/0/g' "$PARAM_INTERVALS"_"$DATE"
	grep -e "#" "$PARAM_INTERVALS"_"$DATE" > "$PARAM_INTERVALS"_"$DATE"_merged
	#bedtools groupby is BROKEN in v2.26.0 . hoping itll get fixed.... this code works correctly for v2.22.0
        $AL_BIN_BEDTOOLS groupby -full -i "$PARAM_INTERVALS"_"$DATE" -g 6 -c 7,8,9 -o mean,mean,mean >> "$PARAM_INTERVALS"_"$DATE"_merged
        printProgress "[Calculating coverage] succesfully finished"
#Methylation
#$AL_USE_BISMARK=1
else
	groomToBed5 "$VAR_OUTPUT_PREFIX1"_methyl.bedGraph
	groomToBed5 "$VAR_OUTPUT_PREFIX2"_methyl.bedGraph
	groomToBed5 "$VAR_OUTPUT_PREFIX3"_methyl.bedGraph
	bedmap_coverage "$PARAM_INTERVALS"_"$DATE" "$VAR_OUTPUT_PREFIX1"_methyl.bed5
	bedmap_coverage "$PARAM_INTERVALS"_"$DATE" "$VAR_OUTPUT_PREFIX2"_methyl.bed5
	bedmap_coverage "$PARAM_INTERVALS"_"$DATE" "$VAR_OUTPUT_PREFIX3"_methyl.bed5
	OPERATION="count"
        bedmap_coverage "$PARAM_INTERVALS"_"$DATE" "$VAR_OUTPUT_PREFIX1"_methyl.bed5
        bedmap_coverage "$PARAM_INTERVALS"_"$DATE" "$VAR_OUTPUT_PREFIX2"_methyl.bed5
        bedmap_coverage "$PARAM_INTERVALS"_"$DATE" "$VAR_OUTPUT_PREFIX3"_methyl.bed5
	printProgress "[Calculating coverage] succesfully finished"
fi
