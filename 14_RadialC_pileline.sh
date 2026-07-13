#!/bin/bash

# ============================================================================
# author: lxk
# date: 2026-02-09
# Radial-C 数据分析流程 v2.0
# 支持多样本并行处理，模块化设计
# Usage: bash radialc_pipeline.sh <sample_name> <linker_type> <read1_file> <read2_file> [options]
# ============================================================================

cutadapt=/home/lxk/private/software/miniconda3/envs/cutadapt_env/bin/cutadapt
trim_galore=/home/lxk/private/software/miniconda3/bin/trim_galore
flash=/home/share/glx/script/flash
seqkit=/mnt/disk1/6/lxk/miniconda3/bin/seqkit

# 输入参数
readonly SAMPLE_NAME=${1:-"DNaseC_K562_NaOH_R1r1"}
readonly LINKER=${2:-"Slinker"}
readonly READ1=${3:-"${SAMPLE_NAME}_1.fq.gz"}
readonly READ2=${4:-"${SAMPLE_NAME}_2.fq.gz"}
readonly OUTPUT_DIR=${5:-"$(pwd)/output_${SAMPLE_NAME}/${LINKER}"}

# ===========================================
if [ ${LINKER} == "Slinker" ]; then
    find_linker1="GCCCGGTNNACGCCCG"
    find_linker2="CGGGCGTNNACCGGGC"
    single_trim_linker1="GCCCGGTNNACGCCCGT"
    single_trim_linker2="CGGGCGTNNACCGGGCT"
    paired_trim_linker1="AGCCCGGTNNACGCCCGT"
    paired_trim_linker2="ACGGGCGTNNACCGGGCT"
elif [ ${LINKER} == "linker60bp" ]; then
    find_linker1="CGTCACCAACTTCCTTGGTTGAGTGTGAGCAGGTGGAGAAGGTTGTGGACCACTCATTG"
    find_linker2="CAATGAGTGGTCCACAACCTTCTCCACCTGCTCACACTCAACCAAGGAAGTTGGTGACG"
    single_trim_linker1="CGTCACCAACTTCCTTGGTTGAGTGTGAGCAGGTGGAGAAGGTTGTGGACCACTCATTGT"
    single_trim_linker2="CAATGAGTGGTCCACAACCTTCTCCACCTGCTCACACTCAACCAAGGAAGTTGGTGACGT"
    paired_trim_linker1="ACGTCACCAACTTCCTTGGTTGAGTGTGAGCAGGTGGAGAAGGTTGTGGACCACTCATTGT"
    paired_trim_linker2="ACAATGAGTGGTCCACAACCTTCTCCACCTGCTCACACTCAACCAAGGAAGTTGGTGACGT"
else
    find_linker1="GCCCGGTNNACGCCCG"
    find_linker2="CGGGCGTNNACCGGGC"
    single_trim_linker1="GCCCGGTNNACGCCCGT"
    single_trim_linker2="CGGGCGTNNACCGGGCT"
    paired_trim_linker1="AGCCCGGTNNACGCCCGT"
    paired_trim_linker2="ACGGGCGTNNACCGGGCT"
fi

# ============================================
mkdir -p "${OUTPUT_DIR}/01_fastqc"
fastqc ${READ1} -o "${OUTPUT_DIR}/01_fastqc/"
fastqc ${READ2} -o "${OUTPUT_DIR}/01_fastqc/"
rm "${OUTPUT_DIR}/01_fastqc/${SAMPLE_NAME}_1_fastqc.zip" "${OUTPUT_DIR}/01_fastqc/${SAMPLE_NAME}_2_fastqc.zip"

# ============================================
mkdir -p "${OUTPUT_DIR}/02_trimAdapt"
${trim_galore} -j 7 -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired ${READ1} ${READ2} --gzip -o "${OUTPUT_DIR}/02_trimAdapt/" --path_to_cutadapt ${cutadapt}
fastqc "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_1_val_1.fq.gz" -o "${OUTPUT_DIR}/01_fastqc/"
fastqc "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_2_val_2.fq.gz" -o "${OUTPUT_DIR}/01_fastqc/"
rm "${OUTPUT_DIR}/01_fastqc/${SAMPLE_NAME}_1_val_1_fastqc.zip" "${OUTPUT_DIR}/01_fastqc/${SAMPLE_NAME}_2_val_2_fastqc.zip"

# ============================================
mkdir -p "${OUTPUT_DIR}/03_merge"
# cd "${OUTPUT_DIR}/03_merge"
${flash} "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_1_val_1.fq.gz" "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_2_val_2.fq.gz" -z -M 150 -t 10 -o "output_${SAMPLE_NAME}/${LINKER}/03_merge/${SAMPLE_NAME}" 2>&1 | tee "${OUTPUT_DIR}/03_merge/${SAMPLE_NAME}.log"
# perl ../../${SAMPLE_NAME}_merge_S-lnknum.pl

# ${flash} output_RadialC_K562_2h_R1r1/02_trimAdapt/RadialC_K562_2h_R1r1_1_val_1.fq.gz output_RadialC_K562_2h_R1r1/02_trimAdapt/RadialC_K562_2h_R1r1_2_val_2.fq.gz -z -M 150 -t 10 -o output_RadialC_K562_2h_R1r1/03_merge/RadialC_K562_2h_R1r1 2>&1 | tee output_RadialC_K562_2h_R1r1/03_merge/RadialC_K562_2h_R1r1.log

# ============================================
mkdir -p "${OUTPUT_DIR}/04_trimLk"
${cutadapt} --pair-filter=both --discard-untrimmed --action=trim -O 16 -a ${find_linker1} -a ${find_linker2} -A ${find_linker1} -A ${find_linker2} -j 10 --minimum-length=0 -o "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_trimLk_1.fq.gz" -p "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_trimLk_2.fq.gz" "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_1_val_1.fq.gz" "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_2_val_2.fq.gz" > "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_report.txt"

${seqkit} fx2tab -n -l -j 80 "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_trimLk_1.fq.gz" > "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_1.txt"
${seqkit} fx2tab -n -l -j 80 "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_trimLk_2.fq.gz" > "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_2.txt"

# -------single-------
paste -d "\t" "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_1.txt" "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_2.txt" | awk -v OFS="\t" '{if($2<=10 && $4>10)print $1}' > "${OUTPUT_DIR}/04_trimLk/single_ID_1.txt"
paste -d "\t" "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_1.txt" "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_2.txt" | awk -v OFS="\t" '{if($2>10 && $4<=10)print $3}' > "${OUTPUT_DIR}/04_trimLk/single_ID_2.txt"

${seqkit} grep -j 80 -n -f "${OUTPUT_DIR}/04_trimLk/single_ID_1.txt" "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_1_val_1.fq.gz" -o "${OUTPUT_DIR}/04_trimLk/single_1.fq.gz"
${seqkit} grep -j 80 -n -f "${OUTPUT_DIR}/04_trimLk/single_ID_2.txt" "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_2_val_2.fq.gz" -o "${OUTPUT_DIR}/04_trimLk/single_2.fq.gz"

${cutadapt} --strip-suffix /1 -g "${single_trim_linker1};rightmost" -g "${single_trim_linker2};rightmost" --minimum-length=10 -j 20 "${OUTPUT_DIR}/04_trimLk/single_1.fq.gz" -o "${OUTPUT_DIR}/04_trimLk/single_trimLk_1.fq.gz" > "${OUTPUT_DIR}/04_trimLk/single_report1.txt"
${cutadapt} --strip-suffix /2 -g "${single_trim_linker1};rightmost" -g "${single_trim_linker2};rightmost" --minimum-length=10 -j 20 "${OUTPUT_DIR}/04_trimLk/single_2.fq.gz" -o "${OUTPUT_DIR}/04_trimLk/single_trimLk_2.fq.gz" > "${OUTPUT_DIR}/04_trimLk/single_report2.txt"
cat "${OUTPUT_DIR}/04_trimLk/single_trimLk_1.fq.gz" "${OUTPUT_DIR}/04_trimLk/single_trimLk_2.fq.gz" > "${OUTPUT_DIR}/04_trimLk/single_trimLk.fq.gz"

# -----paired-------
paste -d "\t" "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_1.txt" "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_2.txt" | awk -v OFS="\t" '{if($2>10 && $4>10)print $1}' > "${OUTPUT_DIR}/04_trimLk/paired_ID_1.txt"
paste -d "\t" "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_1.txt" "${OUTPUT_DIR}/04_trimLk/${SAMPLE_NAME}_length_2.txt" | awk -v OFS="\t" '{if($2>10 && $4>10)print $3}' > "${OUTPUT_DIR}/04_trimLk/paired_ID_2.txt"

${seqkit} grep -j 80 -n -f "${OUTPUT_DIR}/04_trimLk/paired_ID_1.txt" "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_1_val_1.fq.gz" -o "${OUTPUT_DIR}/04_trimLk/paired_1.fq.gz"
${seqkit} grep -j 80 -n -f "${OUTPUT_DIR}/04_trimLk/paired_ID_2.txt" "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_2_val_2.fq.gz" -o "${OUTPUT_DIR}/04_trimLk/paired_2.fq.gz"

${cutadapt} --strip-suffix /1 --strip-suffix /2 -a ${paired_trim_linker1} -a ${paired_trim_linker2} -A ${paired_trim_linker1} -A ${paired_trim_linker2} -j 10 --minimum-length=10 -o "${OUTPUT_DIR}/04_trimLk/paired_trimLk_1.fq.gz" -p "${OUTPUT_DIR}/04_trimLk/paired_trimLk_2.fq.gz" "${OUTPUT_DIR}/04_trimLk/paired_1.fq.gz" "${OUTPUT_DIR}/04_trimLk/paired_2.fq.gz" > "${OUTPUT_DIR}/04_trimLk/paired_report.txt"

# ============================================
mkdir -p "${OUTPUT_DIR}/05_mapping"
mkdir -p "${OUTPUT_DIR}/05_mapping/single"
mkdir -p "${OUTPUT_DIR}/05_mapping/paired"

# ---------------------
bowtie2 -t -p 80 --no-unal --very-sensitive --score-min L,0,-0.4 -x /ssd/index/bowtie2/hg38XX -U "output_${SAMPLE_NAME}/${LINKER}/04_trimLk/single_trimLk.fq.gz" 2> "${OUTPUT_DIR}/05_mapping/single/single.stat" | samtools view -@ 80 -bh -q 10 -o "${OUTPUT_DIR}/05_mapping/single/single_trimLk_q10.bam" -
bedtools bamtobed -i "${OUTPUT_DIR}/05_mapping/single/single_trimLk_q10.bam" | sort -k1,1 -k2,2n -k3,3n -k5,5n -k6,6 -S 20% --parallel 80 | awk '{if($1!~/chr[CLMT]/) print $0}' > "${OUTPUT_DIR}/05_mapping/single/single_trimLk_rmdup.bed"

# ---------------------
mkdir -p "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/data/paired"
# mkdir HiCPro_Analysis/data
# mkdir HiCPro_Analysis/data/trimLk

ln -sf "${OUTPUT_DIR}/04_trimLk/paired_trimLk_1.fq.gz" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/data/paired/paired_trimLk_1.fq.gz"
ln -sf "${OUTPUT_DIR}/04_trimLk/paired_trimLk_2.fq.gz" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/data/paired/paired_trimLk_2.fq.gz"

cp /home/lxk/private/optionData/DNaseC_config_hicpro_XX.txt "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/config_hicpro.txt"

# cd "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis"
sed -i 's/trimLk3/trimLk/g' "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/config_hicpro.txt"

hic-pro -c "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/config_hicpro.txt" -i "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/data" -o "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results" -s mapping -s quality_checks
hic-pro -c "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/config_hicpro.txt" -i "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/bowtie_results/bwt2" -o "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results" -s proc_hic -s merge_persample -s quality_checks


# ============================================================
total=$(zcat ${READ1} | echo $((`wc -l`/4)))
trimAdapt=$(zcat "${OUTPUT_DIR}/02_trimAdapt/${SAMPLE_NAME}_1_val_1.fq.gz" | echo $((`wc -l`/4)))
single=$(cat "${OUTPUT_DIR}/04_trimLk/single_1.fq.gz" "${OUTPUT_DIR}/04_trimLk/single_2.fq.gz" | zcat | echo $((`wc -l`/4)))
paired=$(zcat "${OUTPUT_DIR}/04_trimLk/paired_1.fq.gz" | echo $((`wc -l`/4)))



single_trimLk=$(zcat "${OUTPUT_DIR}/04_trimLk/single_trimLk.fq.gz" | echo $((`wc -l`/4)))
single_unmapped=$(grep "aligned 0 times" "${OUTPUT_DIR}/05_mapping/single/single.stat")
single_lowq10=$(bedtools bamtobed -i "${OUTPUT_DIR}/05_mapping/single/single_trimLk_q10.bam" | wc -l)
single_rmdup=$(wc -l "${OUTPUT_DIR}/05_mapping/single/single_trimLk_rmdup.bed")



paired_trimLk=$(zcat "${OUTPUT_DIR}/04_trimLk/paired_trimLk_1.fq.gz" | echo $((`wc -l`/4)))
paired_unmapped=$(grep "Unmapped_pairs" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired.mpairstat")
paired_lowq10=$(grep "Low_qual_pairs" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired.mpairstat")
paired_singleton=$(grep "Pairs_with_singleton" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired.mpairstat")
paired_Reported=$(grep "Reported_pairs" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired.mpairstat")
paired_rmdup=$(grep "valid_interaction_rmdup" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired_allValidPairs.mergestat")
paired_cis=$(grep "cis_interaction" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired_allValidPairs.mergestat")
paired_trans=$(grep "trans_interaction" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired_allValidPairs.mergestat")
paired_cis_short=$(grep "cis_shortRange" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired_allValidPairs.mergestat")
paired_cis_long=$(grep "cis_longRange" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired_allValidPairs.mergestat")
paired_FR=$(grep "Valid_interaction_pairs_FR" "${OUTPUT_DIR}/05_mapping/paired/HiCPro_Analysis/hicpro_results/hic_results/stats/paired/paired.mRSstat")

echo -e "total\t${total}\ntrimAdapt\t${trimAdapt}\nsingle\t${single}\npaired\t${paired}" > "${OUTPUT_DIR}/technical_stat.txt"
echo "===============================" >> "${OUTPUT_DIR}/technical_stat.txt"

echo -e "single_trimLk\t${single_trimLk}\nsingle_unmapped\t${single_unmapped}\nsingle_lowq10\t${single_lowq10}\nsingle_rmdup\t${single_rmdup}" >> "${OUTPUT_DIR}/technical_stat.txt"

echo "===============================" >> "${OUTPUT_DIR}/technical_stat.txt"
echo -e "paired_trimLk\t${paired_trimLk}\npaired_unmapped\t${paired_unmapped}\npaired_lowq10\t${paired_lowq10}\npaired_singleton\t${paired_singleton}\npaired_Reported\t${paired_Reported}\npaired_rmdup\t${paired_rmdup}\npaired_cis\t${paired_cis}\npaired_trans\t${paired_trans}\npaired_cis_short\t${paired_cis_short}\npaired_cis_long\t${paired_cis_long}\npaired_FR\t${paired_FR}" >> "${OUTPUT_DIR}/technical_stat.txt"