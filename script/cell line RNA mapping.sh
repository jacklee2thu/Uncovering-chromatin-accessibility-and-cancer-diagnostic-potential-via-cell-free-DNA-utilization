cd /Share2/home/zhaolab/data/gastric_cancer/experiment/RNA/GSE-1_RNA
trim_galore -q 20 --phred33 --stringency 3 --length 20 \
--paired mkn-28-ATAC_FKDL202617109-1a_1.fq.gz mkn-28-ATAC_FKDL202617109-1a_2.fq.gz  \
--gzip -o /Share2/home/lanxun3/jacklee/cfDNA/experiment/ATAC/fastqc \

cd /Share2/home/lanxun3/jacklee/reference_genome
STAR --runMode genomeGenerate --runThreadN 8 \
--genomeDir /Share2/home/lanxun3/jacklee/reference_genome/star_index \
--genomeFastaFiles GRCh37.fa \
--sjdbGTFfile gencode.v36lift37.annotation.gtf \
--sjdbOverhang 99

cd /Share2/home/zhaolab/data/gastric_cancer/experiment/RNA/Sj_normal_RNA
STAR --quantMode GeneCounts --genomeDir /Share2/home/lanxun3/jacklee/reference_genome/star_index \
--runThreadN 2 \
--outFilterMismatchNmax 2 \
--readFilesIn *_1.fq.gz *_2.fq.gz \
--outFileNamePrefix /Share2/home/lanxun3/jacklee/cfDNA/experiment/RNA/Sj_normal_RNA/Sj_normal_RNA \
--outFilterMultimapNmax 1 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate
cd /Share2/home/lanxun3/jacklee/cfDNA/experiment/RNA/RNA_expression
paste GSE-1_RNAReadsPerGene.out.tab HGC_27_RNA_2ReadsPerGene.out.tab HGC_27_RNAReadsPerGene.out.tab GSE_1_RNA_2ReadsPerGene.out.tab \
MKN-28_RNAReadsPerGene.out.tab MKN_28_RNA_2ReadsPerGene.out.tab \
MKN45_RNAReadsPerGene.out.tab MKN_45_RNA_2ReadsPerGene.out.tab \
SGC7901_RNAReadsPerGene.out.tab SGC7901_RNA_2ReadsPerGene.out.tab| \
cut -f1,4,8,12,16,20,24,28,32,36,40 | \
tail -n +5 > tmpfile
cat tmpfile | sed "s/^gene://" >gene_count.txt
