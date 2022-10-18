cd /Share2/home/zhaolab/data/gastric_cancer/experiment/ATAC/MKN28_ATAC
trim_galore -q 20 --phred33 --stringency 3 --length 20 \
--paired mkn-28-ATAC_FKDL202617109-1a_1.fq.gz mkn-28-ATAC_FKDL202617109-1a_2.fq.gz  \
--gzip -o /Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/fastqc \

cd /Share2/home/zhaolab/data/healthy/chenwenyan
fastqc -o /Share2/home/zhaolab/jacklee/temp_script \
-t 8 *_1.fq.gz *_2.fq.gz

wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip##bowtie2索引
cd /Share2/home/lanxun3/jacklee/reference_genome/bowtie2_index
unzip hg19.zip


cd /Share2/home/zhaolab/data/gastric_cancer/experiment/ATAC/MKN28_ATAC
bowtie2 -x /Share2/home/lanxun3/jacklee/reference_genome/bowtie2_index/hg19 \
-1 *_1.fq.gz  -2 *_2.fq.gz | \
samtools sort -@ 10 -O bam -o \
/Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/MKN28_ATAC/MKN28_ATAC_sort.bam
cd /Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/MKN28_ATAC
java -jar /Share2/home/lanxun3/jacklee/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=MKN28_ATAC_sort.bam O=MKN28_ATAC.bam M=MKN28_ATAC.markdup_metrics.txt
samtools view -h -f 2 -q 30  MKN28_ATAC.bam | \
samtools sort  -O bam  -@ 10 -o MKN28_ATAC.last.bam
bedtools bamtobed -i MKN28_ATAC.last.bam > MKN28_ATAC.bed
#samtools view -h -F 1028
cd /Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/MKN28_ATAC
macs2 callpeak -f BEDPE -t MKN28_ATAC.bed -g hs -n MKN28_ATAC_peak -B -q 0.01 \
--outdir ./peaks
