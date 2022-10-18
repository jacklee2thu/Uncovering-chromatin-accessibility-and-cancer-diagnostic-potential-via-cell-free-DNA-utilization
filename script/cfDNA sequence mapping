cd /Share2/home/zhaolab/data/gastric_cancer/experiment/cfDNA/GSE_1_cfDNA
trim_galore -q 20 --phred33 --stringency 3 --length 20 \
--paired *_1.fq.gz *_2.fq.gz  \
--gzip -o /Share2/home/zhaolab/data/gastric_cancer/experiment/cfDNA/GSE_1_cfDNA \

rawdata='/fshare2/lijie/raw_data/PAN_CANCER/res/bladder'
datadir='/fshare2/lijie/cfDNA_transfactor/datasets/bladder_cancer'
for sample_name in TGY-2-TGYY000078_FKDL210224721-1a-D710-AK1545 TGYX000024_FKDL210180222-1a TGYY000033_FKDL210180212-1a \
TGYY000063_FKDL210263038-1a TGYY000063_FKDL210263038-1a_1 TGYY000193_FKDL220004223-1a_HYLWYDSX2_L1

do
mkdir -m 777 ${datadir}/${sample_name}
cd ${rawdata}/${sample_name}
trim_galore -q 20 --phred33 --stringency 3 --length 20 \
--paired *_1.fq.gz *_2.fq.gz  \
--gzip -o ${datadir}/${sample_name}

cd ${datadir}/${sample_name}
bwa mem -t 10 -R '@RG\tID:ST-E00318:816:H3JWJCCX2\tPL:illumina\tLB:library\tSM:${sample_name}' \
/Share2/home/lanxun3/jacklee/reference_genome/GRCh37.fa \
*_1.fq.gz \
*_2.fq.gz|samtools view -bS| \
samtools sort -@ 10 -O bam -o \
${datadir}/${sample_name}/${sample_name}_sort.bam
java -jar /Share2/home/lanxun3/jacklee/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${sample_name}_sort.bam O=${sample_name}.bam \
TMP_DIR=/fshare2/lijie/temp M=${sample_name}.markdup_metrics.txt
samtools view -q 30 ${sample_name}.bam|cut -f 3,4,9 > ${sample_name}.txt

done
rm *.fq.gz_trimming_report.txt
rm *.fq.gz
rm *_sort.bam
