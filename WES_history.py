#! /bin/bash
sample=$1

#mkdir
mkdir ${sample}
cd ${sample}

ln -s /BiO/data/raw_data/${sample}_1.fastq.gz.
ln -s /BiO/data/raw_data/${sample}_2.fastq.gz.

# trimming
/BiO/apps/sickle/sickle pe -t sanger -g -f ${sample}_1.fastq.gz -r ${sample}_2.fastq.gz -o ${sample}-trimmed_1.fastq.gz -p ${sample}-trimmed_2.fastq.gz -s ${sample}-trimmed_3.fastq.gz

#Alignment & BAM sorting
/BiO/apps/bwa-0.7.17/bwa mem -t 2 -M -R "@RG\tID:BioEdu\tSM:${sample}\tPL:illumina\tLB:WES" /BiO/data/reference/hg19.fa ${sample}-trimmed_1.fastq.gz ${sample}-trimmed_2.fastq.gz /BiO/apps/samtools/samtools view -bS -q 20 - | /BiO/apps/samtools/samtools sort -m 4000000000 -o ${sample}.sorted.bam 

#check the bam file 
/BiO/apps/samtools/samtools view -h ${sample}.sorted.bam | head -n 10

#Remove Duplicates
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar MarkDuplicates -| ${sample}.sorted.bam -O ${sample}.rmdup.bam -M ${sample}.rmdup.metrics --REMOVE_DUPLICATES=true

/BiO/apps/samtools/samtools index ${sample}.rmdup.bam

#Base quality score recalibration
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar BaseRecalibrator -R /BiO/data/reference/hg19.fasta -I ${sample}.rmdup.bam -L /BiO/data/target/target.bed --known-sites /BiO/data/DB/dbSnp151_chr.vcf.gz --known-sites /BiO/data/DB/mills_and_1000G_gold_standard.indels.hg19.vcf.gz -O Na12878_recal_data.table

#Apply BQSR
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyBQSR -bqsr ${sample}_recal_data.table -I ${sample}.remdup.bam -O ${sample}.recal.bam

/BiO/apps/samtools/samtools index ${sample}.recal.bam

#On-target coverage
java -Xmx4g -jar /BiO/apps/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T DepthOfCoverage -R /BiO/data/reference/hg19.fasta -I ${sample}.recal.bam -o ${sample}_target_cov -ct 1 -ct 5 -ct 10 -ct 20 -ct 30 -omitBaseOutput -L /BiO/data/target/target.bed

#remove off-target reads
bedtools intersect -wa -a ${sample}.recal.bam -b /BiO/data/target/target.bed > ${sample}.target.bam

# Variant calling GATK HaplotypeCaller
java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller -R /BiO/data/reference/hg19.fasta -I ${sample}.target.bam -O ${sample}.gatk.vcf


# Variant calling - samtools (bcftools)
/BiO/apps/bcftools/bcftools mpileup -f /BiO/data/reference/hg19.fasta ${sample}.target.bam | /BiO/apps/bcftools/bcftools call -mv -Ov -o ${sample}.samt.vcf

bgzip ${sample}.samt.vcf
tabix -p vcf ${sample}.samt.vcf.gz

# consensus VCF & filteration
vcf-isec -o -n +2 ${sample}.gatk.vcf.gz ${sample}.samt.vcf.gz > ${sample}.consensus.vcf

java -Xmx4g -jar /BiO/apps/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration -R /BiO/data/reference/hg19.fasta -O ${sample}.consensus.filt.vcf --variant ${sample}.consensus.vcf --filter-expression 'DP < 10 || FS > 60.0' --filter-name 'LOWQUAL'

cat ${sample}.consensus.filt.vcf | awk -F '\t'($7!="LOWQUAL"){print}' | bgzip > ${sample}.final.vcf.gz

tabix -p vcf ${sample}.final.vcf.gz


