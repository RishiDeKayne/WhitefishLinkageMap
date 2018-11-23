#!/bin/bash
#BSUB -J "GATK"
#BSUB -R "rusage[mem=10000]"
#BSUB -n 5
#BSUB -W 90:00

samples=""
for i in *.bam; do samples=$samples" -I "$i; done
java -Xmx16G \
-jar /cluster/project/gdc/shared/p129/bin/GenomeAnalysisTK.v3.7.jar \
-T HaplotypeCaller $samples \
-R /cluster/project/gdc/shared/p298/Rishi/parents/Cor_DenovoParentsQ30.fasta \
--genotyping_mode DISCOVERY \
--output_mode EMIT_ALL_CONFIDENT_SITES \
--min_base_quality_score 20 \
--standard_min_confidence_threshold_for_calling 20 \
-o ConcatParentsQ30Thresh20HaploCalled_Total.vcf
