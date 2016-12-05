filename=$1

# remove adaptors
#echo "...removing adaptors" > run_"$filename".log
#java -jar /home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 4 \
#-trimlog trimmolog1.txt \
#"$filename"_R1_raw.fastq \
#"$filename"_R2_raw.fastq \
#"$filename"_R1_raw_trimmo_paired_truseq3-PE-2_2_10_5_1.fastq \
#"$filename"_R1_raw_trimmo_unpaired_truseq3-PE-2_2_10_5_1.fastq \
#"$filename"_R2_raw_trimmo_paired_truseq3-PE-2_2_10_5_1.fastq \
#"$filename"_R2_raw_trimmo_unpaired_truseq3-PE-2_2_10_5_1.fastq \
#ILLUMINACLIP:/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:10:5:1 SLIDINGWINDOW:4:20 MINLEN:50 2> trimmolog2.txt

#  map with bowtie
echo "...mapping with bowtie" >> run_"$filename".log
bowtie2 -x /home/Reference_genomes/Arabidopsis_thaliana_Ensembl_TAIR10/Ensembl/TAIR10/Sequence/Bowtie2Index/genome \
-p 4 \
-U "$filename"_R1_raw.fastq \
-S "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed.sam \
--maxins 800 \
--no-mixed --no-discordant --no-unal 2>> run_"$filename".log

# bam and sort
echo "...converting sam to bam" >> run_"$filename".log
samtools view -bS "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed.sam \
> "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed.bam 2>> run_"$filename".log
echo "...sorting mapped reads" >> run_"$filename".log
samtools sort "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed.bam "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted  2>> run_"$filename".log

# mark and remove duplicates 
echo "...removing duplicates" >> run_"$filename".log
java -Xmx4g -jar /home/Program_NGS_sl-pw-srv01/picard-tools-1.103/MarkDuplicates.jar \
INPUT= "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted.bam \
OUTPUT="$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
METRICS_FILE=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> markdup_stderr.txt

# indexing
samtools index "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam

# word count raw data
raw_line=$(wc -l < "$filename"_R1_raw.fastq) 
raw_count=$(echo "$raw_line/4" | bc -l)
echo "number of raw reads: $raw_count" >> run_"$filename".log

# get flagstat
echo "number of clean reads mapped:" >> run_"$filename".log
samtools flagstat "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam >> run_"$filename".log

# estimate genome average and normalise
echo "...normalising reads" >> run_"$filename".log
genomeCoverageBed -split -bg -ibam "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
-g /home/Reference_genomes/Arabidopsis_thaliana_Ensembl_TAIR10/Ensembl/TAIR10/Annotation/Genes/ChromInfo.txt \
> "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bedgraph

sum=$(samtools depth "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam| awk '{sum+=$3;cnt++}END{printf "%.0f", sum}')
sum_norm=$(echo "$sum/119667750" | bc -l)
echo "genome normalised coverage: $sum_norm" >> run_"$filename".log

export MYVAR=$sum_norm
perl -e 'print $ENV{MYVAR}."\n"'

# normalise read counts by genome-wide coverage
perl -ne 'chomp($_); @a=split(/\t/,$_);print $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]/$ENV{MYVAR}."\n";' \
"$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bedgraph \
> "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bedgraph

# convert bedgraph to bigwig
bedGraphToBigWig "$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bedgraph \
/home/Reference_genomes/Arabidopsis_thaliana_Ensembl_TAIR10/Ensembl/TAIR10/Annotation/Genes/ChromInfo.txt \
"$filename"_raw_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bw
