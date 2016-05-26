#------------------------------------------------------------------------------------
# This ALEA modification replaces the required BWA aligner with the new STAR aligner.
# It aligns F1 hybrid ChIP-seq datasets to a pseudogenome generated with ALEA and a
# list of known SNPs and INDELs unique to the parental strains.
# It then uses ALEA for creating parental BIGWIGs projected onto the reference for
# visualization.
#------------------------------------------------------------------------------------
module load bioinf/STAR/2.4.0i
module load bioinf/bedtools/2.22.1
mkdir $1 dupeMetrics WIGs
mkdir WIGs/mm9
cd $1
touch $1_stats.txt
echo "$1 alignment statistics file" >> $1_stats.txt #output metrics
#------------------------------------------------------------------------------------
# CONFIGURE
#------------------------------------------------------------------------------------
THREADS=16
RAM_MAX=72000000000 #72GB
PATH_TO_FASTQ="/brcwork/lorincz_lab/jbrind/ML102/forALLSTAR"
PATH_TO_GENOME="/brcwork/lorincz_lab/jrichardalbert/ALLSTAR/references/BL6_CAST"
PATH_TO_BEDGRAPH_TO_WIG="/brcwork/lorincz_lab/jrichardalbert/scripts"

#------------------------------------------------------------------------------------
# PAIRED END MODE FOR ChIP-SEQ SAMPLES
#------------------------------------------------------------------------------------
/brcwork/bioinf/tools/samtools-0.1.16/samtools sort -n $PATH_TO_FASTQ/$1.bam $1_sort.bam
mv $1_sort.bam.bam $1_sort.bam
bedtools bamtofastq -i $1_sort.bam -fq $1_1.fq -fq2 $1_2.fq
bgzip -c $1_1.fq > $1_1.fq.gz
bgzip -c $1_2.fq > $1_2.fq.gz
rm $1_sort.bam
rm $1_1.fq
rm $1_2.fq

STAR --runMode alignReads \
--genomeDir $PATH_TO_GENOME \
--runThreadN $THREADS \
--limitBAMsortRAM $RAM_MAX \
--readFilesCommand zcat \
--readFilesIn $1_1.fq.gz $1_2.fq.gz \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--outSAMtype SAM
echo "Pseudogenome_alignment_STAR_log" >> $1_stats.txt
cat Log.final.out >> $1_stats.txt #output metrics

#------------------------------------------------------------------------------------
# UNPACK ALLELIC READS
#------------------------------------------------------------------------------------
samtools view -SH Aligned.out.sam | awk '($0 ~ /C57BL6J/) {print $0}' | sed 's/C57BL6J_chr/chr/g' > $1_C57BL6J.sam
samtools view -SH Aligned.out.sam | awk '($0 ~ /CAST_EiJ/) {print $0}' | sed 's/CAST_EiJ_chr/chr/g' > $1_CAST_EiJ.sam

samtools view Aligned.out.sam | awk '(($3 ~ /C57BL6J/)&&($5==255)) {print $0}' | sed 's/C57BL6J_chr/chr/g' >> $1_C57BL6J.sam
samtools view Aligned.out.sam | awk '(($3 ~ /CAST_EiJ/)&&($5==255)) {print $0}' | sed 's/CAST_EiJ_chr/chr/g' >> $1_CAST_EiJ.sam
rm Aligned.out.sam

#------------------------------------------------------------------------------------
# CONVERT TO BAM, SORT, MARKDUPLICATES, INDEX
#------------------------------------------------------------------------------------
samtools view -@ 16 -bt $PATH_TO_GENOME/rousemat/mm10.fa $1_C57BL6J.sam > $1_C57BL6J_unsorted.bam
samtools view -@ 16 -bt $PATH_TO_GENOME/BL6_CAST/CAST_EiJ.fasta $1_CAST_EiJ.sam > $1_CAST_EiJ_unsorted.bam
rm $1_C57BL6J.sam
rm $1_CAST_EiJ.sam

samtools sort -@ 16 -m 3600000000 $1_C57BL6J_unsorted.bam $1_C57BL6J_sort
samtools sort -@ 16 -m 3600000000 $1_CAST_EiJ_unsorted.bam $1_CAST_EiJ_sort

java -jar -Xmx4G /brcwork/bioinf/tools/picard-tools-1.72/MarkDuplicates.jar \
INPUT=$1_C57BL6J_sort.bam OUTPUT=$1_C57BL6J.bam VERBOSITY=ERROR METRICS_FILE="$1_C57BL6J_dupeMetrics.txt" AS=true
echo "$1_C57BL6J_MarkDuplicates" >> $1_stats.txt
cat $1_C57BL6J_dupeMetrics.txt >> $1_stats.txt #output metrics
rm $1_C57BL6J_sort.bam
java -jar -Xmx4G /brcwork/bioinf/tools/picard-tools-1.72/MarkDuplicates.jar \
INPUT=$1_CAST_EiJ_sort.bam OUTPUT=$1_CAST_EiJ.bam VERBOSITY=ERROR METRICS_FILE="$1_CAST_EiJ_dupeMetrics.txt" AS=true
rm $1_CAST_EiJ_sort.bam
echo "$1_CAST_EiJ_MarkDuplicates" >> $1_stats.txt
cat $1_CAST_EiJ_dupeMetrics.txt >> $1_stats.txt #output metrics

echo "$1 C57BL6J flagstats" >> $1_stats.txt
samtools flagstat $1_C57BL6J.bam >> $1_stats.txt #output metrics
echo "$1 CAST_EiJ flagstats" >> $1_stats.txt
samtools flagstat $1_CAST_EiJ.bam >> $1_stats.txt #output metrics


#------------------------------------------------------------------------------------
# CONVERT BAM TO WIG
#------------------------------------------------------------------------------------
samtools view -bh -F 1540 $1_C57BL6J.bam > $1_C57BL6J_filtered.bam
bedtools genomecov -ibam $1_C57BL6J_filtered.bam -bg -split > $1_C57BL6J_preProject.bedGraph
echo "Converting C57BL6J bedGraph to WIG. Will tike time..."
perl $PATH_TO_BEDGRAPH_TO_WIG/bedGraphToFixedStepWig.pl --bedgraph $1_C57BL6J_preProject.bedGraph --step 1 --compact --wig $1_C57BL6J_preProject.wig
gzip -c $1_C57BL6J_preProject.wig > $1_C57BL6J_preProject.wig.gz

echo "Converting CAST_EiJ bedGraph to WIG. Will tike time..."
samtools view -bh -F 1540 $1_CAST_EiJ.bam > $1_CAST_EiJ_filtered.bam
bedtools genomecov -ibam $1_CAST_EiJ_filtered.bam -bg -split > $1_CAST_EiJ_preProject.bedGraph
perl $PATH_TO_BEDGRAPH_TO_WIG/bedGraphToFixedStepWig.pl --bedgraph $1_CAST_EiJ_preProject.bedGraph --step 1 --compact --wig $1_CAST_EiJ_preProject.wig
gzip -c $1_CAST_EiJ_preProject.wig > $1_CAST_EiJ_preProject.wig.gz

#------------------------------------------------------------------------------------
# PROJECT TO THE REFERENCE
#------------------------------------------------------------------------------------
java -Xms4G -Xmx8G -jar /brcwork/lorincz_lab/jrichardalbert/alea/bin/alea.jar project \
--input-wig=$1_C57BL6J_preProject.wig.gz \
--input-refmap=$PATH_TO_GENOME/C57BL6J.fasta.refmap \
--output-bedgraph=$1_C57BL6J.bedGraph #liftover to mm9, bigwig, output

java -Xms4G -Xmx8G -jar /brcwork/lorincz_lab/jrichardalbert/alea/bin/alea.jar project \
--input-wig=$1_CAST_EiJ_preProject.wig.gz \
--input-refmap=$PATH_TO_GENOME/CAST_EiJ.fasta.refmap \
--output-bedgraph=$1_CAST_EiJ.bedGraph #liftover to mm9, bigwig, output

/brcwork/bioinf/scripts/bedGraphToBigWig $1_C57BL6J.bedGraph $PATH_TO_GENOME/mm10_chr_sizes.sizes $1_C57BL6J.bw #output
/brcwork/bioinf/scripts/bedGraphToBigWig $1_CAST_EiJ.bedGraph $PATH_TO_GENOME/mm10_chr_sizes.sizes $1_CAST_EiJ.bw #output




#------------------------------------------------------------------------------------
# NON-ALLELIC ALIGNMENT
#------------------------------------------------------------------------------------
PATH_TO_MM10="/brcwork/lorincz_lab/jrichardalbert/STAR_bams/references/mm10"
#------------------------------------------------------------------------------------
# PAIRED END MODE
#------------------------------------------------------------------------------------
STAR --runMode alignReads \
--genomeDir $PATH_TO_MM10 \
--runThreadN $THREADS \
--limitBAMsortRAM $RAM_MAX \
--readFilesCommand zcat \
--readFilesIn $1_1.fq.gz $1_2.fq.gz \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--outSAMtype SAM
echo "reference_alignment_STAR_log" >> $1_stats.txt
cat Log.final.out >> $1_stats.txt #output metrics

#------------------------------------------------------------------------------------
# CONVERT TO BAM, SORT, MARKDUPLICATES, INDEX
#------------------------------------------------------------------------------------
samtools view -@ 16 -bt $PATH_TO_MM10/mm10.fa Aligned.out.sam > $1_total_unsorted.bam
samtools sort -@ 16 -m 3600000000 $1_total_unsorted.bam $1_total_sort
rm Aligned.out.sam
java -jar -Xmx4G /brcwork/bioinf/tools/picard-tools-1.72/MarkDuplicates.jar \
INPUT=$1_total_sort.bam OUTPUT=$1_total.bam VERBOSITY=ERROR METRICS_FILE="$1_dupeMetrics.txt" \
AS=true REMOVE_DUPLICATES=false
echo "$1_ignoreAllelic_MarkDuplicates" >> $1_stats.txt
cat $1_dupeMetrics.txt >> $1_stats.txt #output metrics

#------------------------------------------------------------------------------------
# FIX FALSE POSITIVE DUPLICATES (ALLELIC READS), FILTER OUT TRUE DUPLICATES, CONVERT TO BIGWIG
#------------------------------------------------------------------------------------
samtools view $1_C57BL6J.bam | awk '{print $1}' >> "readNames.txt"
samtools view $1_CAST_EiJ.bam | awk '{print $1}' >> "readNames.txt"
samtools view -H $1_total.bam > $1_total_fixedDupes.sam
#check if allelic read is marked as duplicate in total.
#these are false positives, so remove flag duplicate flag then filter on F1540.
samtools view $1_total.bam > $1_total.sam
awk 'BEGIN { OFS = "\t" } {
    if ( NR == FNR ) {
        a[$1]++;
    } else {
	if ( a[$1] > 0 ) {
            if( and($2,1024) == 1024 ) {
                $2-=1024;
            }
	}
	print $0;
    }
}' readNames.txt $1_total.sam >> $1_total_fixedDupes.sam
samtools view -@ 16 -bt /brcwork/lorincz_lab/jrichardalbert/STAR_bams/references/mm10/mm10.fa $1_total_fixedDupes.sam > $1_total_fixedDupes.bam
echo "$1 fixed duplicates flagstats" >> $1_stats.txt
samtools flagstat $1_total_fixedDupes.bam >> $1_stats.txt #output metrics
rm $1_total.sam
rm $1_fixedDupes.sam

samtools view -bh -F 1540 $1_total_fixedDupes.bam > $1_total_fixedDupes_F1540.bam
bedtools genomecov -ibam $1_total_fixedDupes_F1540.bam -bg -split > $1_total.bedGraph #lifover, convert to bw, output
/brcwork/bioinf/scripts/bedGraphToBigWig $1_total.bedGraph $PATH_TO_GENOME/mm10_chr_sizes.sizes $1_total.bw #output

#------------------------------------------------------------------------------------
# CONVERT BAM TO NORMALIZED BIGWIG
#------------------------------------------------------------------------------------
#RPM_SCALING_FACTOR=`samtools view -c $1_total_fixedDupes_F1540.bam` #get total read number
#RPM_SCALING_FACTOR_1=$(echo "scale=25;1000000/$RPM_SCALING_FACTOR" | bc) #Basic Calculator, 25 decimal points
#echo "Scaling Factor=" $RPM_SCALING_FACTOR_1 >> Log.final.out
#bedtools genomecov -ibam $1_total_fixedDupes_F1540.bam -bg -split -scale $RPM_SCALING_FACTOR_1 > $1_total_RPM.bedGraph
#/brcwork/bioinf/scripts/bedGraphToBigWig $1_total_RPM.bedGraph $PATH_TO_GENOME/mm10_chr_sizes.sizes $1_total_RPM.bw

#------------------------------------------------------------------------------------
# LIFTOVER FROM mm10 TO mm9
#------------------------------------------------------------------------------------

sed '1d' $1_total.bedGraph > $1_total_tmp.bedGraph #remove first line of bedGraph file
sort -k1,1 -k2,2n $1_total_tmp.bedGraph > $1_total_sort.bedGraph #sort by start site
/brcwork/bioinf/scripts/liftOver $1_total_sort.bedGraph \
/brcwork/lorincz_lab/jrichardalbert/ALLSTAR/references/BL6_CAST/mm10ToMm9.over.chain \
$1_total_mm9_tmp.bedGraph trash.txt #convert mm10 to mm9, discarding unmapped regions
grep -v 'random' $1_total_mm9_tmp.bedGraph > $1_total_mm9_tmp_noRandomChr.bedGraph #remove random chromosomes
sort -k1,1 -k2,2n $1_total_mm9_tmp_noRandomChr.bedGraph > $1_total_mm9_tmp_noRandomChr_sort.bedGraph #sort by start site
bedtools merge -i $1_total_mm9_tmp_noRandomChr_sort.bedGraph -d -1 -c 4 -o sum > $1_total_mm9_tmp_merge.bedGraph #add scores at overlapping intervals
/brcwork/bioinf/scripts/bedGraphToBigWig $1_total_mm9_tmp_merge.bedGraph /brcwork/bioinf/reference/genomes/Mus_musculus/mm9/mm9.chrom.sizes $1_total_mm9.bw #output
rm *tmp*
sed '1d' $1_C57BL6J.bedGraph > $1_C57BL6J_tmp.bedGraph #remove first line of bedGraph file
sort -k1,1 -k2,2n $1_C57BL6J_tmp.bedGraph > $1_C57BL6J_sort.bedGraph #sort by start site
/brcwork/bioinf/scripts/liftOver $1_C57BL6J_sort.bedGraph \
/brcwork/lorincz_lab/jrichardalbert/ALLSTAR/references/BL6_CAST/mm10ToMm9.over.chain \
$1_C57BL6J_mm9_tmp.bedGraph trash.txt #convert mm10 to mm9, discarding unmapped regions
grep -v 'random' $1_C57BL6J_mm9_tmp.bedGraph > $1_C57BL6J_mm9_tmp_noRandomChr.bedGraph #remove random chromosomes
sort -k1,1 -k2,2n $1_C57BL6J_mm9_tmp_noRandomChr.bedGraph > $1_C57BL6J_mm9_tmp_noRandomChr_sort.bedGraph #sort by start site
bedtools merge -i $1_C57BL6J_mm9_tmp_noRandomChr_sort.bedGraph -d -1 -c 4 -o sum > $1_C57BL6J_mm9_tmp_merge.bedGraph #add scores at overlapping intervals
/brcwork/bioinf/scripts/bedGraphToBigWig $1_C57BL6J_mm9_tmp_merge.bedGraph /brcwork/bioinf/reference/genomes/Mus_musculus/mm9/mm9.chrom.sizes $1_C57BL6J_mm9.bw #output
rm *tmp*
sed '1d' $1_CAST_EiJ.bedGraph > $1_CAST_EiJ_tmp.bedGraph #remove first line of bedGraph file
sort -k1,1 -k2,2n $1_CAST_EiJ_tmp.bedGraph > $1_CAST_EiJ_sort.bedGraph #sort by start site
/brcwork/bioinf/scripts/liftOver $1_CAST_EiJ_sort.bedGraph \
/brcwork/lorincz_lab/jrichardalbert/ALLSTAR/references/BL6_CAST/mm10ToMm9.over.chain \
$1_CAST_EiJ_mm9_tmp.bedGraph trash.txt #convert mm10 to mm9, discarding unmapped regions
grep -v 'random' $1_CAST_EiJ_mm9_tmp.bedGraph > $1_CAST_EiJ_mm9_tmp_noRandomChr.bedGraph #remove random chromosomes
sort -k1,1 -k2,2n $1_CAST_EiJ_mm9_tmp_noRandomChr.bedGraph > $1_CAST_EiJ_mm9_tmp_noRandomChr_sort.bedGraph #sort by start site
bedtools merge -i $1_CAST_EiJ_mm9_tmp_noRandomChr_sort.bedGraph -d -1 -c 4 -o sum > $1_CAST_EiJ_mm9_tmp_merge.bedGraph #add scores at overlapping intervals
/brcwork/bioinf/scripts/bedGraphToBigWig $1_CAST_EiJ_mm9_tmp_merge.bedGraph /brcwork/bioinf/reference/genomes/Mus_musculus/mm9/mm9.chrom.sizes $1_CAST_EiJ_mm9.bw #output
rm *tmp*
#------------------------------------------------------------------------------------
# MOVE AROUND
#------------------------------------------------------------------------------------
mv $1_stats.txt ../dupeMetrics/
mv $1_total.bw ../WIGs/
mv $1_C57BL6J.bw ../WIGs/
mv $1_CAST_EiJ.bw ../WIGs/
mv $1_total_mm9.bw ../WIGs/mm9/
mv $1_C57BL6J_mm9.bw ../WIGs/mm9/
mv $1_CAST_EiJ_mm9.bw ../WIGs/mm9/

mkdir ../Track_Hub
mkdir ../Track_Hub/bedGraphs
mv $1_C57BL6J.bedGraph ../Track_Hub/bedGraphs/
mv $1_CAST_EiJ.bedGraph ../Track_Hub/bedGraphs/
mv $1_total_fixedDupes_F1540.bam ../Track_Hub/
cd ..
rm -r $1
#------------------------------------------------------------------------------------
# GENERATE TRACK HUB
#------------------------------------------------------------------------------------

