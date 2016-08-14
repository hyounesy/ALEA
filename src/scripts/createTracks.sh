#!/bin/bash

pushd `dirname $0` > /dev/null
AL_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

source $AL_DIR_TOOLS/alea.config


if [ $AL_USE_BISMARK = 0 ]; then

if test $# -ne 8
then
    echo "
Usage:   
         alea createTracks <-s/-p> bamPrefix strain1 strain2 genome1.refmap genome2.refmap chrom.sizes outputDIR
         
Options:
         -s to create tracks for the single-end aligned reads
         -p to create tracks for the paired-end aligned reads

         bamPrefix      prefix used for the output of alignReads command
         strain1        name of strain1 (e.g. hap1)
         strain2        name of strain2 (e.g. hap2)
         genome1.refmap path to the refmap file created for insilico genome 1
         genome1.refmap path to the refmap file created for insilico genome 2
         chrom.sizes    path to the chromosome size file (required for creating .bw)
         outputDIR      output directory (where to create track files)
         
Output:
         outputDIR/outputPrefix_strain1.bedGraph
         outputDIR/outputPrefix_strain1.bw        read profiles for strain1 projected to reference genome
         
         outputDIR/outputPrefix_strain2.bedGraph 
         outputDIR/outputPrefix_strain2.bw        read profiles for strain2 projected to reference genome
         
         outputDIR/outputPrefix_strain1.wig.gz
         outputDIR/outputPrefix_strain2.wig.gz    unprojected read profiles for strain1 and strain2
"
exit 1
fi

else #[ $AL_USE_BISMARK = 1 ]

if test $# -ne 9
then
    echo "
Usage:   
         alea createTracks <-s/-p> bamPrefix strain1 strain2 genome1.refmap genome2.refmap chrom.sizes outputDIR
         
Options:
         -s to create tracks for the single-end aligned reads
         -p to create tracks for the paired-end aligned reads

         bamPrefix      prefix used for the output of alignReads command
         strain1        name of strain1 (e.g. hap1)
         strain2        name of strain2 (e.g. hap2)
         genome1.refmap path to the refmap file created for insilico genome 1
         genome1.refmap path to the refmap file created for insilico genome 2
         chrom.sizes    path to the chromosome size file (required for creating .bw)
         referenceDIR   path to the genome folder (required for calling methylation)
         outputDIR      output directory (where to create track files)
         
Output:
         outputDIR/outputPrefix_strain1.bedGraph
         outputDIR/outputPrefix_strain1.bw        read profiles for strain1 projected to reference genome
         
         outputDIR/outputPrefix_strain2.bedGraph 
         outputDIR/outputPrefix_strain2.bw        read profiles for strain2 projected to reference genome
         
         outputDIR/outputPrefix_strain1.wig.gz
         outputDIR/outputPrefix_strain2.wig.gz    unprojected read profiles for strain1 and strain2
"
exit 1
fi

fi

##############################################################################
#############   Module 4: projection to reference genome
##############################################################################

### Converts filtered bam files to wig
function convertBam2WigSE {
    local PARAM_INPUT_PREFIX=$1
    local PARAM_OUTPUT_DIR=$2

    local VAR_q=$AL_BAM2WIG_PARAM_MIN_QUALITY     # min read quality [0]
    local VAR_F=$AL_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag [0]
    local VAR_x=$AL_BAM2WIG_PARAM_SE_EXTENSION    # average fragment length used for fixed length of the read extension [0]. used for ChIP-seq (SET) only
    local VAR_INPUT_BASENAME=`basename $PARAM_INPUT_PREFIX`
    
    aleaCheckFileExists "$PARAM_INPUT_PREFIX".bam

    # Create a wig profile from the bam file
    $AL_BIN_BAM2WIG \
        -samtools $AL_BIN_SAMTOOLS \
        -bamFile "$PARAM_INPUT_PREFIX".bam \
        -out $PARAM_OUTPUT_DIR/ \
        -q $VAR_q \
        -F $VAR_F \
        -cs \
        -x $VAR_x
    
    mv $PARAM_OUTPUT_DIR/$VAR_INPUT_BASENAME.q"$VAR_q".F"$VAR_F".SET_"$VAR_x".wig.gz $PARAM_OUTPUT_DIR/$VAR_INPUT_BASENAME.wig.gz
}

function convertBam2WigPE {
    local PARAM_INPUT_PREFIX=$1
    local PARAM_OUTPUT_DIR=$2

    local VAR_q=$AL_BAM2WIG_PARAM_MIN_QUALITY     # min read quality [0]
    local VAR_F=$AL_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag [0]
    local VAR_x=$AL_BAM2WIG_PARAM_SE_EXTENSION    # average fragment length used for fixed length of the read extension [0]. used for ChIP-seq (SET) only
    local VAR_INPUT_BASENAME=`basename $PARAM_INPUT_PREFIX`

    aleaCheckFileExists "$PARAM_INPUT_PREFIX".bam
    
    # Create a wig profile from the bam file
    $AL_BIN_BAM2WIG \
        -samtools $AL_BIN_SAMTOOLS \
        -bamFile "$PARAM_INPUT_PREFIX".bam \
        -out $PARAM_OUTPUT_DIR/ \
        -q $VAR_q \
        -F $VAR_F \
        -cp \
        -x $VAR_x

    mv "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".q"$VAR_q".F"$VAR_F".PET.wig.gz "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".wig.gz
}


### projects a wig profile to reference genome
function projectToReferenceGenome {
    local PARAM_WIG_FILE=$1
    local PARAM_REFMAP_FILE=$2
    local PARAM_BEDGRAPH_FILE=$3
    
    aleaCheckFileExists "$PARAM_WIG_FILE"
    aleaCheckFileExists "$PARAM_REFMAP_FILE"

    printProgress "Started projectToReferenceGenome"

    $AL_BIN_ALEA project\
        --input-wig=$PARAM_WIG_FILE\
        --input-refmap=$PARAM_REFMAP_FILE\
        --output-bedgraph=$PARAM_BEDGRAPH_FILE
    
    printProgress "Finished projectToReferenceGenome"
}

### merge a cytosine report into a CpG site report
function mergeTwoStrandMethylation {
    local PARAM_CYTESINE_REPORT_FILE=$1
    local PARAM_SITE_REPORT_FILE=$2
    
    aleaCheckFileExists "$PARAM_CYTESINE_REPORT_FILE"

    printProgress "Started mergeTwoStrandMethylation"

    awk '
        BEGIN{
            FS = "\t"
            OFS = "\t"
            
            FIRST_POS = 0
            FIRST_METHYL = 0
            FIRST_UNMETHYL = 0
            METHYL = 0
            UNMETHYL = 0
            FIRST_TRI = ""
        }
        $3 == "+"{
            FIRST_POS = $2
            FIRST_METHYL = $4
            FIRST_UNMETHYL = $5
            FIRST_TRI = $7
        }
        $3 == "-"{
            if ($2 == FIRST_POS + 1) {
                METHYL = FIRST_METHYL + $4
                UNMETHYL = FIRST_UNMETHYL + $5
                
                if (METHYL + UNMETHYL > 0) {
                    printf $1 "\t" FIRST_POS "\t" $2 "\t"
                    printf "%6f\t", METHYL / (METHYL + UNMETHYL) * 100.0
                    print METHYL, UNMETHYL, $6, FIRST_TRI
                }
                else {
                    print $1, FIRST_POS, $2, "NA", METHYL, UNMETHYL, $6, FIRST_TRI
                }
            }
            FIRST_POS = 0
            FIRST_METHYL = 0
            FIRST_UNMETHYL = 0
            FIRST_TRI = ""
            METHYL = 0
            UNMETHYL = 0
        }
    ' "$PARAM_CYTESINE_REPORT_FILE" > "$PARAM_SITE_REPORT_FILE"
    
    printProgress "Finished mergeTwoStrandMethylation"
}

### detect the allelic cytosines in the concatenated genome method
function detectAllelicCytoConcatenated {
    
    printProgress "Started detectAllelicCytoConcatenated"

    local PARAM_INPUT_METHYL=$1
    local PARAM_STRAIN=$2
    local PARAM_OUT_PREFIX=$3
    
    aleaCheckFileExists "$PARAM_INPUT_METHYL"
    
    cat "$PARAM_INPUT_METHYL" \
        | awk -v ref="$PARAM_STRAIN" '($0 ~ ref) {print $0}' \
        | sed 's/'"$PARAM_STRAIN"'_chr//g' \
        >> "$PARAM_OUT_PREFIX".CpG_report.txt
    
    printProgress "Finished detectAllelicCytoConcatenated"
}

### convert a CpG site report into a Wig
function convertMethylationToWig {
    local PARAM_SITE_REPORT_FILE=$1
    local PARAM_WIG_FILE=$2

    local VAR_MIN_DEPTH=$AL_METH2WIG_PARAM_MIN_DEPTH
    
    aleaCheckFileExists "$PARAM_SITE_REPORT_FILE"

    printProgress "Started convertMethylationToWig"

    awk -v MIN_DEPTH=$VAR_MIN_DEPTH '
        BEGIN{
            FS = "\t"
            OFS = "\t"
            
            if (MIN_DEPTH < 1)
                MIN_DEPTH = 1
            
            print "track type=wiggle_0"
        }
        $5 + $6 >= MIN_DEPTH{
            if ($1 ~ /^chr/)
                print "fixedStep chrom=" $1 " start=" $2 " step=1 span=1"
            else
                print "fixedStep chrom=chr" $1 " start=" $2 " step=1 span=1"
	        for(i = 0; i < $3-$2+1; i++) {
	        	print $4
	        }
        }
    ' "$PARAM_SITE_REPORT_FILE" > "$PARAM_WIG_FILE"
    
    printProgress "Finished convertMethylationToWig"
}


VAR_OPTION=$1
shift

#function generateAllelicTracks {
    PARAM_BAM_PREFIX=$1
    PARAM_STRAIN1=$2
    PARAM_STRAIN2=$3
    PARAM_REFMAP_FILE1=$4
    PARAM_REFMAP_FILE2=$5
    PARAM_CHROM_SIZES=$6
    PARAM_OUTPUT_DIR=$7
    if [ $AL_USE_BISMARK = 0 ]
    then
        PARAM_OUTPUT_DIR=$7
    else
        PARAM_REFERECE_DIR=$7
        PARAM_OUTPUT_DIR=$8
    fi
    
    aleaCheckFileExists "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1".bam
    aleaCheckFileExists "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2".bam
    aleaCheckFileExists "$PARAM_REFMAP_FILE1"
    aleaCheckFileExists "$PARAM_REFMAP_FILE2"
    aleaCheckFileExists "$PARAM_CHROM_SIZES"
    aleaCreateDir "$PARAM_OUTPUT_DIR"
    
    
    
    if [ "$VAR_OPTION" = "-s" ]; then
        convertBam2WigSE "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1" "$PARAM_OUTPUT_DIR"
        convertBam2WigSE "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2" "$PARAM_OUTPUT_DIR"
    elif [ "$VAR_OPTION" = "-p" ]; then
        convertBam2WigPE "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1" "$PARAM_OUTPUT_DIR"
        convertBam2WigPE "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2" "$PARAM_OUTPUT_DIR"
    else
        echo "Invalid option $VAR_OPTION"
        exit 1
    fi
    
    VAR_OUTPUT_BASENAME=`basename $PARAM_BAM_PREFIX`
    VAR_OUTPUT_PREFIX1="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"
    VAR_OUTPUT_PREFIX2="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"

    projectToReferenceGenome "$VAR_OUTPUT_PREFIX1".wig.gz "$PARAM_REFMAP_FILE1" "$VAR_OUTPUT_PREFIX1".bedGraph
    $AL_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX1".bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX1".bw

    projectToReferenceGenome "$VAR_OUTPUT_PREFIX2".wig.gz "$PARAM_REFMAP_FILE2" "$VAR_OUTPUT_PREFIX2".bedGraph
    $AL_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX2".bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX2".bw

#}

if [ $AL_USE_BISMARK = 1 ]; then
#function generateAllelicMethylTracks {
    aleaCheckFileExists "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
    aleaCheckDirExists "$PARAM_REFERECE_DIR"
    
    mv "$VAR_OUTPUT_PREFIX1".wig.gz "$VAR_OUTPUT_PREFIX1"_coverage.wig.gz
    mv "$VAR_OUTPUT_PREFIX1".bedGraph "$VAR_OUTPUT_PREFIX1"_coverage.bedGraph
    mv "$VAR_OUTPUT_PREFIX1".bw "$VAR_OUTPUT_PREFIX1"_coverage.bw
    mv "$VAR_OUTPUT_PREFIX2".wig.gz "$VAR_OUTPUT_PREFIX2"_coverage.wig.gz
    mv "$VAR_OUTPUT_PREFIX2".bedGraph "$VAR_OUTPUT_PREFIX2"_coverage.bedGraph
    mv "$VAR_OUTPUT_PREFIX2".bw "$VAR_OUTPUT_PREFIX2"_coverage.bw
    
    if [ $AL_USE_CONCATENATED_GENOME = 1 ]; then
    
        samtools view -Sb -q 1 "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.bam
        
        if [ "$VAR_OPTION" = "-s" ]; then
            $AL_BIN_BISMARK_EXTRACT -s --comprehensive --cytosine_report --genome_folder $PARAM_REFERECE_DIR "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.bam
            
        elif [ "$VAR_OPTION" = "-p" ]; then
            $AL_BIN_BISMARK_EXTRACT -p --comprehensive --cytosine_report --genome_folder $PARAM_REFERECE_DIR "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.bam
            
        else
            echo "Invalid option $VAR_OPTION"
            exit 1
        fi
        
        sort -k1,1 -k2,2n "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.CpG_report.txt > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.CpG_report.txt.tmp
        if [ $AL_DEBUG = 0 ]
        then
            rm "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.CpG_report.txt
        else
            mv "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.CpG_report.txt "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.CpG_report.unsorted.txt
        fi
        mv "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.CpG_report.txt.tmp "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.CpG_report.txt
        
        detectAllelicCytoConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.CpG_report.txt "$PARAM_STRAIN1" "$VAR_OUTPUT_PREFIX1"
        detectAllelicCytoConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_filtered.CpG_report.txt "$PARAM_STRAIN2" "$VAR_OUTPUT_PREFIX2"
        
        mergeTwoStrandMethylation "$VAR_OUTPUT_PREFIX1".CpG_report.txt "$VAR_OUTPUT_PREFIX1".CpG_site_report.txt
        mergeTwoStrandMethylation "$VAR_OUTPUT_PREFIX2".CpG_report.txt "$VAR_OUTPUT_PREFIX2".CpG_site_report.txt
        
        convertMethylationToWig "$VAR_OUTPUT_PREFIX1".CpG_site_report.txt "$VAR_OUTPUT_PREFIX1"_methylation.wig
        convertMethylationToWig "$VAR_OUTPUT_PREFIX2".CpG_site_report.txt "$VAR_OUTPUT_PREFIX2"_methylation.wig
        
        projectToReferenceGenome "$VAR_OUTPUT_PREFIX1"_methylation.wig "$PARAM_REFMAP_FILE1" "$VAR_OUTPUT_PREFIX1"_methylation.bedGraph
        $AL_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX1"_methylation.bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX1"_methylation.bw
        
        projectToReferenceGenome "$VAR_OUTPUT_PREFIX2"_methylation.wig "$PARAM_REFMAP_FILE2" "$VAR_OUTPUT_PREFIX2"_methylation.bedGraph
        $AL_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX2"_methylation.bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX2"_methylation.bw
        
    
    else # [ $AL_USE_CONCATENATED_GENOME != 1 ] 
    
        echo -e "\nThe separate method using Bismark for alignment is not supported so far.\n"
        exit 1
    
    fi
#}
fi



### converts SORTED bam into wig format using bedtools (circumventing BAM2WIG)
# input sorted bam and filter variables
### output zipped wig for projection to reference

function convertBam2Wigbedtools {
    printProgress "Started convertBam2Wigbedtools"
    local PARAM_INPUT_PREFIX=$1

    local VAR_q=$AL_BAM2WIG_PARAM_MIN_QUALITY     # min read quality [0]
    local VAR_F=$AL_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag [0]
    aleaCheckFileExists "$PARAM_INPUT_PREFIX".bam
    
    $AL_BIN_SAMTOOLS view -bh -F $VAR_F -q $VAR_q "VAR_INPUT_BASENAME".bam > "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bam
    $AL_BIN_BEDTOOLS genomecov -ibam "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bam -bg -split > "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bedGraph
    
    awk '
        BEGIN {
                print "track type=wiggle_0"
        }
        NF == 4 {
                print "fixedStep chrom="$1" start="$2+1" step=1 span=1"
                for(i = 0; i < $3-$2; i++) {
                        print $4
                }
        }' "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bedGraph > "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".wig
    bgzip -c "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".wig > "$PARAM_INPUT_PREFIX".wig.gz
}


### input allelic projected bedGraphs
# counts allelic reads over exons and aggregates over genes
### output table with allelic read counts
function countAllelicReads {
    
    printProgress "[countAllelicReads] Started"
    local PARAM_PROJECTED_BEDGRAPH=$1
    local PARAM_EXON_COORDINATES=$2
    local PARAM_GENE_COORDINATES=$3
    local REF="C57BL6J"
    local ALT="CAST_EiJ"
    local AL_LOG="."

    echo "Input allelic bedGraphs: "$PARAM_PROJECTED_BEDGRAPH"_"$REF".bedGraph and "$PARAM_PROJECTED_BEDGRAPH"_"$ALT".bedGraph (must be projected onto reference)" >> "$AL_LOG"/run_log.tsv
    echo "Using exon coordinate file: $PARAM_EXON_COORDINATES (ensure no duplicate exons are present in file)" >> "$AL_LOG"/run_log.tsv
    echo "And gene coordinate file: $PARAM_GENE_COORDINATES" >> "$AL_LOG"/run_log.tsv

    $AL_BIN_BEDTOOLS intersect -a "$PARAM_PROJECTED_BEDGRAPH"_"$REF".bedGraph -b "$PARAM_EXON_COORDINATES" > READS_OVERLAPPING_EXONS_"$REF".bedGraph
    $AL_BIN_BEDTOOLS coverage -a "$PARAM_GENE_COORDINATES" -b READS_OVERLAPPING_EXONS_"$REF".bedGraph > "$PARAM_PROJECTED_BEDGRAPH"_"$REF"_count_tmp.tsv
    rm READS_OVERLAPPING_EXONS_"$REF".bedGraph

    $AL_BIN_BEDTOOLS intersect -a "$PARAM_PROJECTED_BEDGRAPH"_"$ALT".bedGraph -b "$PARAM_EXON_COORDINATES" > READS_OVERLAPPING_EXONS_"$ALT".bedGraph
    $AL_BIN_BEDTOOLS coverage -a "$PARAM_GENE_COORDINATES" -b READS_OVERLAPPING_EXONS_"$ALT".bedGraph > "$PARAM_PROJECTED_BEDGRAPH"_"$ALT"_count_tmp.tsv
    rm READS_OVERLAPPING_EXONS_"$ALT".bedGraph

    echo "#chr  start   end     strand  ID      name    exonLength      "$REF"_allelicReadCount" > tmp1.tsv
    awk 'OFS="\t" {print $1, $2, $3, $4, $5, $6, $7, $8}' "$PARAM_PROJECTED_BEDGRAPH"_"$REF"_count_tmp.tsv >> tmp1.tsv
    echo ""$ALT"_allelicReadCount" > tmp2.tsv
    cut -f8 "$PARAM_PROJECTED_BEDGRAPH"_"$ALT"_count_tmp.tsv >> tmp2.tsv
    cut -f1 tmp2.tsv | paste tmp1.tsv - > "$PARAM_PROJECTED_BEDGRAPH"_allelic_read_counts.tsv
    rm tmp1.tsv tmp2.tsv "$PARAM_PROJECTED_BEDGRAPH"_"$REF"_count_tmp.tsv "$PARAM_PROJECTED_BEDGRAPH"_"$ALT"_count_tmp.tsv

    echo "[countAllelicReads] Ended" >> "$AL_LOG"/run_log.tsv
    date >> "$AL_LOG"/run_log.tsv

    #make histogram of allelic reads (REF+ALT) over each gene
    #make 
	
}
countAllelicReads BC_ICM_H3K36me3 Refseq_exon_mm10.bed Refseq_VisRseq_Genes_mm10.bed







### converts SORTED bam into wig format using bedtools (circumventing BAM2WIG)
# input sorted bam and filter variables
### output zipped wig for visualization and RPKM values for analysis
function calculateTotalRPKM {
    printProgress "Started calculateTotalRPKM"
    local PARAM_INPUT_PREFIX=$1
    local PARAM_EXON_COORDINATES=$2
    local PARAM_GENE_COORDINATES=$3
    local VAR_q=$AL_BAM2WIG_PARAM_MIN_QUALITY     # min read quality
    local VAR_F=$AL_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag
 
    local AL_LOG="."
    aleaCheckFileExists "$PARAM_INPUT_PREFIX".bam
    
    #filter bam
    $AL_BIN_SAMTOOLS view -bh -F $VAR_F -q $VAR_q "$PARAM_INPUT_PREFIX".bam > "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bam
    #get scaling factor
    RPM_SCALING_FACTOR_tmp=`samtools view -c "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bam`
    RPM_SCALING_FACTOR=$(echo "scale=25;1000000/$RPM_SCALING_FACTOR_tmp" | bc)
    echo ""$PARAM_INPUT_PREFIX" scaling factor: "$RPM_SCALING_FACTOR"" >> "$AL_LOG"/run_log.tsv

    #convert bam to bedGraph
    $AL_BIN_BEDTOOLS genomecov -ibam "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bam -bg -split -scale "$RPM_SCALING_FACTOR" > "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bedGraph
    #convert bedGraph to bigWig for visualiztion
    $AL_BIN_BEDGRAPHTOBIGWIG "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bedGraph $AL_BIN_SIZES/sizes "$PARAM_INPUT_PREFIX"_total.q"$VAR_q".F"$VAR_F".bw #output this file for viz
    #get exon-overlapping intervals from bedGraph
    $AL_BIN_BEDTOOLS intersect -a "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bedGraph -b "$PARAM_EXON_COORDINATES" > READS_OVERLAPPING_EXONS_total.bedGraph
    #calculate coverage over genes
    $AL_BIN_BEDTOOLS coverage -a "$PARAM_GENE_COORDINATES" -b READS_OVERLAPPING_EXONS_total.bedGraph > "$PARAM_INPUT_PREFIX"_tmp.tsv

    echo "#chr  start   end     strand  ID      name    exonLength      RPKM" > tmp1.tsv
    awk 'OFS="\t" {print $1, $2, $3, $4, $5, $6, $7, $8/$7}' "$PARAM_INPUT_PREFIX"_tmp.tsv >> "$PARAM_INPUT_PREFIX"_RPKM.tsv

    rm tmp1.tsv "$PARAM_INPUT_PREFIX"_tmp.tsv READS_OVERLAPPING_EXONS_total.bedGraph "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bedGraph "$PARAM_INPUT_PREFIX".q"$VAR_q".F"$VAR_F".bam
    
    echo "[calculateTotalRPKM] Ended" >> "$AL_LOG"/run_log.tsv
    date >> "$AL_LOG"/run_log.tsv

}

calculateTotalRPKM BCliverRep1 Refseq_exon_mm10.bed Refseq_VisRseq_mm10.bed
