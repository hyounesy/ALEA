#-----------------------------------------------------------#
#							MODULE 5						#
#-----------------------------------------------------------#	
# generates .tsv data files containing allelic coverage and #
# total RPKM values over user-define intervals. These int-  #
# ervals can be gene exons, gene promoters, known enhancers # 
# CpG islands, etc.											#
# Requirements: indexed reference genome, interval file in  #
# gff formats (General Feature Format).						#
#-----------------------------------------------------------#

#!/bin/bash

pushd `dirname $0` > /dev/null
AL_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

source $AL_DIR_TOOLS/alea.config

##############################################################################
#############   Module 5: statistical analysis
##############################################################################

PARAM_VALID=1
PARAM_SINGLE_READS=1
if [ $AL_USE_CONCATENATED_GENOME = 1 ]; then
    if [ "$1" = "-s" ]; then
        if [ $# -eq 6 ]; then
            PARAM_FASTQ_FILE=$2
            PARAM_GENOME=$3
			PARAM_STRAIN1=$4
			PARAM_STRAIN2=$5
            PARAM_INTERVALS=$6
            PARAM_FLAG=$7
            PARAM_MINMAPQ=$8
			PARAM_PREFIX="${PARAM_FASTQ_FILE##*/}"
        else
            PARAM_VALID=0
        fi
    elif [ "$1" = "-p" ]; then
        if [ $# -eq 7 ]; then
            PARAM_SINGLE_READS=0
            PARAM_FASTQ_FILE1=$2
            PARAM_FASTQ_FILE2=$3
            PARAM_GENOME=$4
			PARAM_STRAIN1=$5
			PARAM_STRAIN2=$6
            PARAM_INTERVALS=$7
            PARAM_FLAG=$8
            PARAM_MINMAPQ=$9
			PARAM_PREFIX="${PARAM_FASTQ_FILE_1##*/}"
        else
            PARAM_VALID=0
        fi
    else
        PARAM_VALID=0
    fi
fi

if [ $PARAM_VALID = 0 ]; then
    echo "
Usage:
    using concatenated genome method (AL_USE_CONCATENATED_GENOME=1):
        alea alignReads <-s/-p> <input_reads_1 [input_reads_2]> <genome_concat> <strain1 strain2> <outputPrefix>

Options:
    -s              to align single-end reads (requires one input file)
    -p              to align paired-end reads (requires two input files)
         
    input_reads_1   the 1st input reads file in fastq.
					will be used as basename for all subsequent files in this module
                    
    input_reads_2   (paired end) the 2nd input reads file in fastq.
					will be used as basename for all subsequent files in this module
                    
    genome   		path to the indexed REFERENCE genome.
                    for BWA, specifiy path to the fasta
                    for Bowtie2 and Tophat2, specify basename of index file
                    for Bismark, specify genome folder
					for STAR, specify genome folder

	intervals		General Feature Format (GFF) intervals over which to measure
					allelic coverage, total RPKM, and other statistics

	-F 				FLAG to filter-out reads (an F value of 1028 will remove reads
					that are unmapped OR are PCR duplicates/optical duplicates).
	-q				Minimum mapping quality (aligner-dependent, see examples). 


Output:
        
    outputPrefix_total.bam				all reads aligned to the reference genome (sorted bam)

    outputPrefix_total_FXXXX_qXX.bam	reads aligned to the reference genome that pass FLAG and MAPQ filters

	outputPrefix_analysis.tsv			a table with allelic coverage, total RPKM, parental ratios
										and other statistics for each interval in the GFF file.

Examples:
    (AL_USE_BWA=1)
    alea createReport -s H3K36me3.fastq reference_genome.fasta transcriptionStarSites.gff -F 1540 -q 10
    alea createReport -p H3K36me3_1.fastq H3K36me3_2.fastq reference_genome.fasta transcriptionStarSites.gff -F 1540 -q 10
    (AL_USE_STAR=1)
    alea createReport -p totalRNA_1.fastq totalRNA_2.fastq reference_genome.fasta exons.gff -F 1540 -q 2

"
exit 1
fi


if [ $AL_USE_CONCATENATED_GENOME = 1 ]; then
    
    if [ $PARAM_SINGLE_READS = 1 ]; then
        #align single-end reads to reference genome
        aleaCheckFileExists "$PARAM_FASTQ_FILE"

        if [ $AL_USE_BWA = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME"
    		$AL_BIN_BWA aln "$PARAM_GENOME" "$PARAM_FASTQ_FILE" > "$PARAM_PREFIX".sai
    		$AL_BIN_BWA samse "$PARAM_GENOME" "$PARAM_PREFIX".sai $PARAM_FASTQ_FILE > "$PARAM_PREFIX"_total.sam
			$AL_BIN_SAMTOOLS view -bS "$PARAM_PREFIX"_total.sam > "$PARAM_PREFIX"_total.bam
			$AL_BIN_SAMTOOLS index "$PARAM_PREFIX"_total.bam


       	elif [ $AL_USE_BISMARK = 1 ]; then
            aleaCheckDirExists "$PARAM_GENOME"/Bisulfite_Genome
            $AL_BIN_BISMARK $AL_BISMARK_ALN_PARAMS --basename "${PARAM_BAM_PREFIX##*/}"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2" -o "${PARAM_BAM_PREFIX%/*}" "$PARAM_GENOME" $PARAM_FASTQ_FILE
        
       	elif [ $AL_USE_STAR = 1 ]; then
            aleaCheckDirExists "$PARAM_GENOME"
            $AL_BIN_STAR --runMode alignReads --genomeDir "$PARAM_GENOME" "$AL_STAR_ALN_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE"
            mv Aligned.out.sam "$PARAM_PREFIX"_total.sam
			$AL_BIN_SAMTOOLS view -bS "$PARAM_PREFIX"_total.sam > "$PARAM_PREFIX"_total.bam
			$AL_BIN_SAMTOOLS index "$PARAM_PREFIX"_total.bam
            mv Log.out "$PARAM_PREFIX"_total_STAR_RunParameters.tsv
            mv Log.final.out "$PARAM_PREFIX"_total_STAR_AlignmentSummary.tsv
            
        elif [ $AL_USE_TOPHAT2 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME".1.bt2l
            $AL_BIN_TOPHAT2 --no-sort-bam --no-convert-bam "$PARAM_GENOME" $PARAM_FASTQ_FILE > "$PARAM_PREFIX"_total.sam
	    	#mv ./tophat_out/accepted_hits.sam ./"$PARAM_PREFIX".sam
			$AL_BIN_SAMTOOLS view -bS "$PARAM_PREFIX"_total.sam > "$PARAM_PREFIX"_total.bam
			$AL_BIN_SAMTOOLS index "$PARAM_PREFIX"_total.bam

		elif [ $AL_USE_BOWTIE2 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME".1.bt2l
            $AL_BIN_BOWTIE2 -x "$PARAM_GENOME" $PARAM_FASTQ_FILE > "$PARAM_PREFIX"_total.sam
			$AL_BIN_SAMTOOLS view -bS "$PARAM_PREFIX"_total.sam > "$PARAM_PREFIX"_total.bam
			$AL_BIN_SAMTOOLS index "$PARAM_PREFIX"_total.bam
        fi

    else #[ $PARAM_SINGLE_READS = 0 ]
        #align paired-end reads to reference genome
        
        aleaCheckFileExists $PARAM_FASTQ_FILE1
        aleaCheckFileExists $PARAM_FASTQ_FILE2

        if [ $AL_USE_BWA = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME"
    		$AL_BIN_BWA aln "$PARAM_GENOME" "$PARAM_FASTQ_FILE1" > "$PARAM_PREFIX_1".sai
			$AL_BIN_BWA aln "$PARAM_GENOME" "$PARAM_FASTQ_FILE2" > "$PARAM_PREFIX_1".sai
    		$AL_BIN_BWA sampe "$PARAM_GENOME" "$PARAM_PREFIX_1".sai "$PARAM_PREFIX_2".sai "$PARAM_FASTQ_FILE1" "$PARAM_FASTQ_FILE2" > "$PARAM_PREFIX"_total.sam
			$AL_BIN_SAMTOOLS view -bS "$PARAM_PREFIX"_total.sam > "$PARAM_PREFIX"_total.bam
			$AL_BIN_SAMTOOLS index "$PARAM_PREFIX"_total.bam


       	elif [ $AL_USE_BISMARK = 1 ]; then
            aleaCheckDirExists "$PARAM_GENOME"/Bisulfite_Genome
            $AL_BIN_BISMARK $AL_BISMARK_ALN_PARAMS --basename "${PARAM_BAM_PREFIX##*/}"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2" -o "${PARAM_BAM_PREFIX%/*}" "$PARAM_GENOME" $PARAM_FASTQ_FILE
        
       	elif [ $AL_USE_STAR = 1 ]; then
            aleaCheckDirExists "$PARAM_GENOME"
            $AL_BIN_STAR --runMode alignReads --genomeDir "$PARAM_GENOME" "$AL_STAR_ALN_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE1" "$PARAM_FASTQ_FILE2"
            mv Aligned.out.sam "$PARAM_PREFIX"_total.sam
			$AL_BIN_SAMTOOLS view -bS "$PARAM_PREFIX"_total.sam > "$PARAM_PREFIX"_total.bam
			$AL_BIN_SAMTOOLS index "$PARAM_PREFIX"_total.bam
            mv Log.out "$PARAM_PREFIX"_total_STAR_RunParameters.tsv
            mv Log.final.out "$PARAM_PREFIX"_total_STAR_AlignmentSummary.tsv
            
        elif [ $AL_USE_TOPHAT2 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME".1.bt2l
            $AL_BIN_TOPHAT2 --no-sort-bam --no-convert-bam "$PARAM_GENOME" $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2 > "$PARAM_PREFIX"_total.sam
			#mv ./tophat_out/accepted_hits.sam ./"$PARAM_PREFIX".sam
			$AL_BIN_SAMTOOLS view -bS "$PARAM_PREFIX"_total.sam > "$PARAM_PREFIX"_total.bam
			$AL_BIN_SAMTOOLS index "$PARAM_PREFIX"_total.bam

		elif [ $AL_USE_BOWTIE2 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME".1.bt2l
            $AL_BIN_BOWTIE2 -x "$PARAM_GENOME" $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2 > "$PARAM_PREFIX"_total.sam
			$AL_BIN_SAMTOOLS view -bS "$PARAM_PREFIX"_total.sam > "$PARAM_PREFIX"_total.bam
			$AL_BIN_SAMTOOLS index "$PARAM_PREFIX"_total.bam
        fi
    fi
fi






function BAM2WIGbedtools {
	PARAM_PREFIX=$1
	PARAM_CHROMOSOME_SIZES=""

	$AL_BIN_SAMTOOLS view -bh -F "$PARAM_FLAG" -q "$PARAM_MINMAPQ" "$PARAM_PREFIX".bam > "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".bam
    $AL_BIN_BEDTOOLS genomecov -ibam "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".bam -bg -split -scale "$RPM_SCALING_FACTOR" > "$PARAM_PREFIX"l_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".bedGraph
    $AL_BIN_BEDGRAPH_TO_BW "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".bedGraph "$PARAM_CHROMOSOME_SIZES" "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".bw #output this file for viz

    awk '
        BEGIN {
                print "track type=wiggle_0"
        }
        NF == 4 {
                print "fixedStep chrom="$1" start="$2+1" step=1 span=1"
                for(i = 0; i < $3-$2; i++) {
                        print $4
                }
        }' "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".bedGraph > "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".wig
    bgzip -c "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".wig > "$PARAM_PREFIX".wig.gz
}

BAM2WIGbedtools "$PARAM_PREFIX"_total
BAM2WIGbedtools "$PARAM_BAM_PREFIX"_"$STRAIN1"
BAM2WIGbedtools "$PARAM_BAM_PREFIX"_"$STRAIN2"





function TrackHubGenerate {
	local PARAM_INPUT_PREFIX=$1
	local PARAM_CHROMOSOME_SIZES=""
	PARAM_STRAIN1=""
	PARAM_STRAIN2=""
	local PARAM_GENOME=""
	

 	RPM_SCALING_FACTOR_tmp=`samtools view -c "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".bam`
    RPM_SCALING_FACTOR=$(echo "scale=25;1000000/$RPM_SCALING_FACTOR_tmp" | bc)
    echo ""$PARAM_PREFIX" scaling factor: "$RPM_SCALING_FACTOR"" >> "$AL_LOG"/log.txt

	awk -F "\t" -v RPM_SCALE="$RPM_SCALING_FACTOR" 'BEGIN {OFS="\t"; print $0} {
		$4 *= RPM_SCALE
		print $0;}' "$PARAM_PREFIX"_"$STRAIN1".bedGraph > "$PARAM_PREFIX"_"$STRAIN1"_tmp1.bedGraph
	grep -v "type" "$PARAM_PREFIX"_"$STRAIN1"_tmp1.bedGraph > "$PARAM_PREFIX"_"$STRAIN1"_tmp2.bedGraph
	sort -k1,1 -k2,2n "$PARAM_PREFIX"_"$STRAIN1"_tmp2.bedGraph > "$PARAM_PREFIX"_"$STRAIN1"_tmp3.bedGraph
	$AL_BIN_BEDGRAPH_TO_BW "$PARAM_PREFIX"_"$STRAIN1"_tmp3.bedGraph "$PARAM_CHROMOSOME_SIZES" ./Track_Hub/"$PARAM_GENOME"/"$PARAM_PREFIX"_"$STRAIN1".bw

	awk -F "\t" -v RPM_SCALE="$RPM_SCALING_FACTOR" 'BEGIN {OFS="\t"; print $0} {
                $4 *= RPM_SCALE
                print $0;}' "$PARAM_PREFIX"_"$STRAIN2".bedGraph > "$PARAM_PREFIX"_"$STRAIN2"_tmp1.bedGraph
    grep -v "type" "$PARAM_PREFIX"_"$STRAIN2"_tmp1.bedGraph > "$PARAM_PREFIX"_"$STRAIN2"_tmp2.bedGraph
    sort -k1,1 -k2,2n "$PARAM_PREFIX"_"$STRAIN2"_tmp2.bedGraph > "$PARAM_PREFIX"_"$STRAIN2"_tmp3.bedGraph
    $AL_BIN_BEDGRAPH_TO_BW "$PARAM_PREFIX"_"$STRAIN2"_tmp3.bedGraph "$PARAM_CHROMOSOME_SIZES" ./Track_Hub/"$PARAM_GENOME"/"$PARAM_PREFIX"_"$STRAIN2".bw

	awk -F "\t" -v RPM_SCALE="$RPM_SCALING_FACTOR" 'BEGIN {OFS="\t"; print $0} {
                $4 *= RPM_SCALE
                print $0;}' "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ".bedGraph > "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ"_tmp1.bedGraph
    grep -v "type" "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ"_tmp1.bedGraph > "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ"_tmp2.bedGraph
    sort -k1,1 -k2,2n "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ"_tmp2.bedGraph > "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ"_tmp3.bedGraph
    $AL_BIN_BEDGRAPH_TO_BW "$PARAM_PREFIX"_F"$PARAM_FLAG"_q"$PARAM_MINMAPQ"_tmp3.bedGraph "$PARAM_CHROMOSOME_SIZES" ./Track_Hub/"$PARAM_GENOME"/"$PARAM_PREFIX"_total.bw


	printf "genome "$PARAM_GENOME"\ntrackDb "$PARAM_GENOME"/trackDb.txt" > ./Track_Hub/genomes.txt
	printf "hub <name>\nshortLabel <short name>\nlongLabel <Hub to dispaly data at UCSC>\ngenomesFile genomes.txt\nemail <email>" > ./Track_Hub/hub.txt
	printf "track %s\ncontainer multiWig\nshortLabel %s\nlongLabel %s\ntype bigWig\nvisibility full\nmaxHeightPixels 150:60:32\nconfigurable on\nautoScale on\naggregate transparentOverlay\nshowSubtrackColorOnUi on\npriority 1.0\n\n" "$PARAM_PREFIX"_total "$PARAM_PREFIX"_total "$PARAM_PREFIX"_total.bw | tee -a ./Track_Hub/"$PARAM_GENOME"/trackDb.txt
	printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 215,215,215\n\taltColor 225,225,225\n\tautoScale on\n\n" "$PARAM_PREFIX"_total "$PARAM_PREFIX"_total "$PARAM_PREFIX"_total "$PARAM_PREFIX"_total "$PARAM_PREFIX"_total.bw | tee -a ./Track_Hub/"$PARAM_GENOME"/trackDb.txt
	printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 215,215,215\n\taltColor 225,225,225\n\tautoScale on\n\n" "$PARAM_PREFIX"_"$STRAIN1" "$PARAM_PREFIX"_"$STRAIN1" "$PARAM_PREFIX"_"$STRAIN1" "$PARAM_PREFIX"_"$STRAIN1" "$PARAM_PREFIX"_"$STRAIN1".bw | tee -a ./Track_Hub/"$PARAM_GENOME"/trackDb.txt
	printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 215,215,215\n\taltColor 225,225,225\n\tautoScale on\n\n" "$PARAM_PREFIX"_"$STRAIN2" "$PARAM_PREFIX"_"$STRAIN2" "$PARAM_PREFIX"_"$STRAIN2" "$PARAM_PREFIX"_"$STRAIN2" "$PARAM_PREFIX"_"$STRAIN2".bw | tee -a ./Track_Hub/"$PARAM_GENOME"/trackDb.txt
	$AL_BIN_HUBCHECK ./Track_Hub/hub.txt
}





function generateRPKM {
	$AL_BIN_PYTHON $AL_BIN_RPKMCOUNT -i "$PARAM_PREFIX".bam -r "$PARAM_INTERVALS" -o "$PARAM_PREFIX" -q "$PARAM_MINMAPQ" -e
}

function generateRPKMcustom {
#input gff file
#input bedGraph files from 1) bedtools genomecov (total.bedgraph) and 2) alea project (strain1.bedgraph strain2.bedgraph)
#convert to smart format: unique "genes" i.e. no isoforms (done in Python) BY EXON
#bedtools intersect the unique exon and bedgraph files
#using geneName value to merge (SUM) counts over all exons from the same gene
#divide the total.bedgraph column by the calculated mRNA size 
#report the strain1.bedgraph and strain2.bedgraph columns as COVERAGE
#calculate final column: parental_skew / parental_ratio / bias (strain2/(strain1+2)
}

	





