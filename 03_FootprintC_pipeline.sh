#!/bin/bash
set -euo pipefail

# ============================================================================
# author: lxk
# date: 2026-02-09
# DNase-C 数据分析流程 v2.0
# 支持多样本并行处理，模块化设计
# Usage: ./dnasec_pipeline.sh <sample_name> <read1_file> <read2_file> [options]
# ============================================================================

# 配置部分（根据环境修改）
readonly THREADS=${THREADS:-40}
readonly MEMORY=${MEMORY:-100G}
readonly BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
readonly SOFTWARE_DIR="/home/lxk/private/software/miniconda3/bin"
readonly OPTION_DIR="/home/lxk/private/optionData"
readonly GENOME_SIZE="/ssd/genome/hg38_chromsize.txt"
readonly CTCF_MOTIF="${OPTION_DIR}/function_motif/HEK293T_CTCF_RAD21_DNase.motif"

# 工具路径
readonly TRIM_GALORE="${SOFTWARE_DIR}/trim_galore"
readonly CUTADAPT="${SOFTWARE_DIR}/cutadapt"
readonly FLASH="/home/lxk/share/software/flash"
readonly BAMTOBED="bamToBed"
readonly PICARD="/home/lxk/share/software/picard.jar"
readonly HIC_PRO="hic-pro"

# 输入参数
readonly SAMPLE_NAME=${1:-"DNaseC_K562_NaOH_R1r1"}
readonly READ1=${2:-"${SAMPLE_NAME}_1.fq.gz"}
readonly READ2=${3:-"${SAMPLE_NAME}_2.fq.gz"}
# readonly OUTPUT_DIR=${4:-"./output"}
readonly OUTPUT_DIR=${4:-"$(pwd)/output_${SAMPLE_NAME}"}

# 创建输出目录结构
init_directories() {
    local dirs=(
        "01_fastqc" "02_trim_adapt" "03_merge" "04_trim_linker"
        "05_hicpro" "06_optical_dup" "07_fragment_len" "08_dimerization"
        "09_aschip" "logs"
    )
    
    for dir in "${dirs[@]}"; do
        mkdir -p "${OUTPUT_DIR}/${dir}"
    done
    
    exec > >(tee -a "${OUTPUT_DIR}/logs/pipeline_${SAMPLE_NAME}.log")
    exec 2>&1
    
    echo "========================================"
    echo "DNase-C Pipeline Started: $(date)"
    echo "Sample: ${SAMPLE_NAME}"
    echo "R1: ${READ1}"
    echo "R2: ${READ2}"
    echo "Output: ${OUTPUT_DIR}"
    echo "Threads: ${THREADS}"
    echo "========================================"
}

# 步骤 1: FastQC 质控
run_fastqc() {
    echo "[$(date)] Step 1: FastQC Quality Control"
    local input_files=("$@")
    
    for file in "${input_files[@]}"; do
        if [[ ! -f "$file" ]]; then
            echo "ERROR: Input file not found: $file"
            exit 1
        fi
        
        fastqc "$file" -o "${OUTPUT_DIR}/01_fastqc" \
            --threads $((THREADS / 2)) \
            2>&1 | tee "${OUTPUT_DIR}/logs/fastqc_$(basename "$file").log"
    done
    
    # 清理 zip 文件节省空间
    rm -f "${OUTPUT_DIR}/01_fastqc"/*.zip
    echo "[$(date)] FastQC completed"
}

# 步骤 2: Trim Galore 去接头
run_trim_adapt() {
    echo "[$(date)] Step 2: Adapter Trimming with Trim Galore"
    
    ${TRIM_GALORE} \
        -j $((THREADS / 2)) \
        -q 20 \
        --phred33 \
        --stringency 3 \
        --length 20 \
        -e 0.1 \
        --paired "$READ1" "$READ2" \
        --gzip \
        -o "${OUTPUT_DIR}/02_trim_adapt" \
        --path_to_cutadapt "${CUTADAPT}" \
        2>&1 | tee "${OUTPUT_DIR}/logs/trim_adapt.log"
    
    # 定义输出文件名（Trim Galore 命名规则）
    TRIM_R1="${OUTPUT_DIR}/02_trim_adapt/${SAMPLE_NAME}_1_val_1.fq.gz"
    TRIM_R2="${OUTPUT_DIR}/02_trim_adapt/${SAMPLE_NAME}_2_val_2.fq.gz"
    
    echo "[$(date)] Adapter trimming completed"
}

# 步骤 3: FLASH 合并 reads
run_merge() {
    echo "[$(date)] Step 3: Merging Paired Reads with FLASH"
    
    cd "${OUTPUT_DIR}/03_merge"
    
    ${FLASH} \
        "${TRIM_R1}" "${TRIM_R2}" \
        -z \
        -M 150 \
        -t $((THREADS / 4)) \
        -o "${SAMPLE_NAME}" \
        2>&1 | tee "${SAMPLE_NAME}.log"
    
    # S-linker 统计（如果有对应的 perl 脚本）
    if [[ -f "${BASE_DIR}/${SAMPLE_NAME}_merge_S-lnknum.pl" ]]; then
        perl "${BASE_DIR}/${SAMPLE_NAME}_merge_S-lnknum.pl" \
            2>&1 | tee -a "${SAMPLE_NAME}.log"
    fi
    
    cd "${BASE_DIR}"
    echo "[$(date)] Merging completed"
}

# 步骤 4: 去除 S-linker（三轮 cutadapt）
run_trim_linker() {
    echo "[$(date)] Step 4: S-linker Trimming (3 rounds)"
    
    local input_r1="${TRIM_R1}"
    local input_r2="${TRIM_R2}"
    local output_dir="${OUTPUT_DIR}/04_trim_linker"
    
    # Round 1: 去除带数字后缀的 linker
    echo "  Round 1/3: Strip suffix and trim linker"
    ${CUTADAPT} \
        --strip-suffix /1 --strip-suffix /2 \
        -a "file:${OPTION_DIR}/MicroC_S-linker.fa" \
        -A "file:${OPTION_DIR}/MicroC_S-linker.fa" \
        --rename='{id};{r1.adapter_name}' \
        -j $((THREADS / 4)) \
        --minimum-length=10 \
        -o "${output_dir}/${SAMPLE_NAME}_trimLk1_1.fq.gz" \
        -p "${output_dir}/${SAMPLE_NAME}_trimLk1_2.fq.gz" \
        "$input_r1" "$input_r2" \
        > "${output_dir}/report1.txt" 2>&1
    
    # Round 2: 再次去除 linker
    echo "  Round 2/3: Second round trimming"
    ${CUTADAPT} \
        -a "file:${OPTION_DIR}/MicroC_S-linker.fa" \
        -A "file:${OPTION_DIR}/MicroC_S-linker.fa" \
        --rename='{id};{r1.adapter_name}' \
        -j $((THREADS / 4)) \
        --minimum-length=10 \
        -o "${output_dir}/${SAMPLE_NAME}_trimLk2_1.fq.gz" \
        -p "${output_dir}/${SAMPLE_NAME}_trimLk2_2.fq.gz" \
        "${output_dir}/${SAMPLE_NAME}_trimLk1_1.fq.gz" \
        "${output_dir}/${SAMPLE_NAME}_trimLk1_2.fq.gz" \
        > "${output_dir}/report2.txt" 2>&1
    
    # Round 3: 去除 noAT linker
    echo "  Round 3/3: Remove noAT linker"
    ${CUTADAPT} \
        -a "file:${OPTION_DIR}/MicroC_S-linker_noAT.fa" \
        -A "file:${OPTION_DIR}/MicroC_S-linker_noAT.fa" \
        --rename='{id};{r1.adapter_name}' \
        -j $((THREADS / 4)) \
        --minimum-length=10 \
        -o "${output_dir}/${SAMPLE_NAME}_trimLk3_1.fq.gz" \
        -p "${output_dir}/${SAMPLE_NAME}_trimLk3_2.fq.gz" \
        "${output_dir}/${SAMPLE_NAME}_trimLk2_1.fq.gz" \
        "${output_dir}/${SAMPLE_NAME}_trimLk2_2.fq.gz" \
        > "${output_dir}/report3.txt" 2>&1
    
    # 最终输出文件
    FINAL_R1="${output_dir}/${SAMPLE_NAME}_trimLk3_1.fq.gz"
    FINAL_R2="${output_dir}/${SAMPLE_NAME}_trimLk3_2.fq.gz"
    
    # FastQC 检查最终质量
    run_fastqc "$FINAL_R1" "$FINAL_R2"
    
    echo "[$(date)] S-linker trimming completed"
}

# 步骤 5: HiC-Pro 分析
run_hicpro() {
    echo "[$(date)] Step 5: HiC-Pro Analysis"
    
    local hicpro_dir="${OUTPUT_DIR}/05_hicpro"
    local data_dir="${hicpro_dir}/data/${SAMPLE_NAME}_trimLk3"
    
    mkdir -p "$data_dir"
    
    # 复制输入文件
    cp "$FINAL_R1" "$FINAL_R2" "$data_dir/"
    
    # 复制配置文件（或动态生成）
    local config_file="${hicpro_dir}/config_hicpro.txt"
    if [[ -f "${OPTION_DIR}/DNaseC_config_hicpro_XY.txt" ]]; then
        cp "${OPTION_DIR}/DNaseC_config_hicpro_XY.txt" "$config_file"
    else
        echo "ERROR: HiC-Pro config file not found"
        exit 1
    fi
    
    cd "$hicpro_dir"
    
    # Step 1: Mapping
    ${HIC_PRO} \
        -c "$config_file" \
        -i data \
        -o hicpro_results \
        -s mapping \
        -s quality_checks \
        2>&1 | tee "${OUTPUT_DIR}/logs/hicpro_mapping.log"
    
    # Step 2: Processing
    ${HIC_PRO} \
        -c "$config_file" \
        -i hicpro_results/bowtie_results/bwt2 \
        -o hicpro_results \
        -s proc_hic \
        -s merge_persample \
        -s quality_checks \
        2>&1 | tee "${OUTPUT_DIR}/logs/hicpro_processing.log"
    
    cd "${BASE_DIR}"
    
    # 定义 BAM 文件路径
    HICPRO_BAM="${hicpro_dir}/hicpro_results/bowtie_results/bwt2/${SAMPLE_NAME}_trimLk3/${SAMPLE_NAME}__hg38XY.bwt2pairs.bam"
    
    echo "[$(date)] HiC-Pro completed"
}

# 步骤 6: 光学重复去除
run_optical_dup() {
    echo "[$(date)] Step 6: Optical Duplicate Removal"
    
    local dup_dir="${OUTPUT_DIR}/06_optical_dup"
    local input_bam="${HICPRO_BAM}"
    
    if [[ ! -f "$input_bam" ]]; then
        echo "ERROR: HiC-Pro BAM file not found: $input_bam"
        exit 1
    fi
    
    # 排序
    echo "  Sorting BAM..."
    samtools sort -@ "$THREADS" -m 2G "$input_bam" \
        > "${dup_dir}/${SAMPLE_NAME}_sort.bam"
    
    # MarkDuplicates（光学距离 10000）
    echo "  Marking optical duplicates..."
    java -jar "${PICARD}" MarkDuplicates \
        -REMOVE_DUPLICATES false \
        --VALIDATION_STRINGENCY SILENT \
        -I "${dup_dir}/${SAMPLE_NAME}_sort.bam" \
        -O "${dup_dir}/${SAMPLE_NAME}_markdup_DIS10k.bam" \
        -M "${dup_dir}/${SAMPLE_NAME}_DIS10k.metrics" \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 10000 \
        --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 \
        --TAGGING_POLICY All \
        2>&1 | tee "${OUTPUT_DIR}/logs/optical_dup.log"
    
    # 清理中间文件
    rm -f "${dup_dir}/${SAMPLE_NAME}_sort.bam"
    
    echo "[$(date)] Optical duplicate removal completed"
}

# 步骤 7: 片段长度分析
run_fragment_length() {
    echo "[$(date)] Step 7: Fragment Length Analysis"
    
    local frag_dir="${OUTPUT_DIR}/07_fragment_len"
    local input_bam="${HICPRO_BAM}"
    
    # 转换为 bedpe 并处理
    echo "  Converting to BEDPE..."
    ${BAMTOBED} -bedpe -i "$input_bam" | \
        sort -k1,1 -k2,3n -k4,4 -k5,6n -S "${MEMORY}" | \
        awk '{if($1==$4 && $2>$5){print $4 "\t" $5 "\t" $6 "\t" $1 "\t" $2 "\t" $3 "\t" $7 "\t" $8 "\t" $10 "\t" $9}else{print $0}}' | \
        awk -F';' '{print $1 "\t" $2$3$4$5$6$7$8$9}' | \
        sort -k1,6 -k8,8 -k10,11 -u -S "${MEMORY}" | \
        awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/ && $8!="no_adapterno_adapterno_adapter") print $0}' \
        > "${frag_dir}/${SAMPLE_NAME}_rmdup.pair"
    
    # 提取长度
    awk '{print $3-$2"\t"$6-$5}' "${frag_dir}/${SAMPLE_NAME}_rmdup.pair" \
        > "${frag_dir}/${SAMPLE_NAME}_rmdup.length"
    
    # V-plot 分析（CTCF motif）
    echo "  Generating V-plot data..."
    awk '{print $1"\t"int(($2+$3)/2)"\t"int(($2+$3)/2)+1"\t"NR"\t"$3-$2"\t"$10"\n"$4"\t"int(($5+$6)/2)"\t"int(($5+$6)/2)+1"\t"NR"_2""\t"$6-$5"\t"$11}' \
        "${frag_dir}/${SAMPLE_NAME}_rmdup.pair" | \
        sort -k1,1 -k2,2n -S "${MEMORY}" \
        > "${frag_dir}/${SAMPLE_NAME}_sort.singlepair"
    
    closestBed -a "${frag_dir}/${SAMPLE_NAME}_sort.singlepair" \
        -b "$CTCF_MOTIF" -d -t first \
        > "${frag_dir}/${SAMPLE_NAME}_singlepair_closest_CTCFmotif.txt"
    
    awk -f "${OPTION_DIR}/FragmentLengthVsDistance_stat.awk" \
        "${frag_dir}/${SAMPLE_NAME}_singlepair_closest_CTCFmotif.txt" \
        > "${frag_dir}/FragmentLengthVsDistance_${SAMPLE_NAME}.csv"
    
    echo "[$(date)] Fragment length analysis completed"
}

# 步骤 8: Dimerization 分析
run_dimerization() {
    echo "[$(date)] Step 8: Dimerization Analysis"
    
    local dimer_dir="${OUTPUT_DIR}/08_dimerization/${SAMPLE_NAME}"
    mkdir -p "$dimer_dir"
    
    cd "$dimer_dir"
    
    # 分离 up/down stream
    awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$10}' \
        "../../07_fragment_len/${SAMPLE_NAME}_rmdup.pair" | \
        sort -k1,1 -k2,2n -S "${MEMORY}" \
        > "${SAMPLE_NAME}_rmdup_sort.uppair"
    
    awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$11}' \
        "../../07_fragment_len/${SAMPLE_NAME}_rmdup.pair" | \
        sort -k1,1 -k2,2n -S "${MEMORY}" \
        > "${SAMPLE_NAME}_rmdup_sort.downpair"
    
    # 找最近的 CTCF motif
    closestBed -a "${SAMPLE_NAME}_rmdup_sort.uppair" \
        -b "$CTCF_MOTIF" -d -t first \
        > "${SAMPLE_NAME}_rmdup_sort_uppair_closest_CTCFmotif.txt"
    
    closestBed -a "${SAMPLE_NAME}_rmdup_sort.downpair" \
        -b "$CTCF_MOTIF" -d -t first \
        > "${SAMPLE_NAME}_rmdup_sort_downpair_closest_CTCFmotif.txt"
    
    # 清理大文件
    rm -f "${SAMPLE_NAME}_rmdup_sort.uppair" "${SAMPLE_NAME}_rmdup_sort.downpair"
    
    # 合并配对
    sort -k4,4 -S "${MEMORY}" \
        "${SAMPLE_NAME}_rmdup_sort_uppair_closest_CTCFmotif.txt" \
        > "${SAMPLE_NAME}_rmdup_sort_uppair_closest_CTCFmotif_sort.txt"
    
    rm -f "${SAMPLE_NAME}_rmdup_sort_uppair_closest_CTCFmotif.txt"

    sort -k4,4 -S "${MEMORY}" \
        "${SAMPLE_NAME}_rmdup_sort_downpair_closest_CTCFmotif.txt" \
        > "${SAMPLE_NAME}_rmdup_sort_downpair_closest_CTCFmotif_sort.txt"
    
    rm -f "${SAMPLE_NAME}_rmdup_sort_downpair_closest_CTCFmotif.txt"

    paste "${SAMPLE_NAME}_rmdup_sort_uppair_closest_CTCFmotif_sort.txt" \
        "${SAMPLE_NAME}_rmdup_sort_downpair_closest_CTCFmotif_sort.txt" \
        > "${SAMPLE_NAME}_rmdup_pair_closet_CTCFmotif.txt"
    
    rm -f *_sort.txt
    
    # 统计
    awk -f "${OPTION_DIR}/script/DNase-C/stat_noextend.awk" \
        "${SAMPLE_NAME}_rmdup_pair_closet_CTCFmotif.txt" \
        > "${SAMPLE_NAME}_rmdup_pair_closet_CTCFmotif.stat"
    
    cd "${BASE_DIR}"
    echo "[$(date)] Dimerization analysis completed"
}

# 步骤 9: asChIP-seq 分析
run_aschip() {
    echo "[$(date)] Step 9: asChIP-seq Analysis"
    
    local aschip_dir="${OUTPUT_DIR}/09_aschip"
    local plot_dir="${aschip_dir}/plot_Average_Profile"
    mkdir -p "$plot_dir"
    
    # 生成 BigWig
    awk '{print($1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6)}' \
        "${OUTPUT_DIR}/07_fragment_len/${SAMPLE_NAME}_rmdup.pair" \
        > "${aschip_dir}/${SAMPLE_NAME}_rmdup.bed"
    
    local total_reads=$(wc -l < "${aschip_dir}/${SAMPLE_NAME}_rmdup.bed")
    
    sort -k1,1 -k2,2n -k3,3n -S "${MEMORY}" \
        "${aschip_dir}/${SAMPLE_NAME}_rmdup.bed" | \
        genomeCoverageBed -bg -i - -g "$GENOME_SIZE" | \
        awk -v total="$total_reads" '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4/total*1000000}' | \
        awk '{$4/=1;print}' OFS='\t' \
        > "${aschip_dir}/${SAMPLE_NAME}_rmdup.bdg"
    
    bedGraphToBigWig "${aschip_dir}/${SAMPLE_NAME}_rmdup.bdg" \
        "$GENOME_SIZE" \
        "${aschip_dir}/${SAMPLE_NAME}_rmdup.bw"
    
    rm -f "${aschip_dir}/${SAMPLE_NAME}_rmdup.bed" \
        "${aschip_dir}/${SAMPLE_NAME}_rmdup.bdg"
    
    # 绘制平均谱图
    cd "$plot_dir"
    
    computeMatrix reference-point \
        --referencePoint center \
        -S "${aschip_dir}/${SAMPLE_NAME}_rmdup.bw" \
        -R "$CTCF_MOTIF" \
        -a 500 -b 501 -bs 1 \
        --missingDataAsZero \
        -p max \
        -o "${SAMPLE_NAME}_CTCFmotif.gz"
    
    plotProfile \
        -m "${SAMPLE_NAME}_CTCFmotif.gz" \
        -o "${SAMPLE_NAME}_CTCFmotif_profile.png" \
        --outFileNameData "${SAMPLE_NAME}_CTCFmotif_profile.tab" \
        --plotHeight 6 --plotWidth 9 \
        --perGroup \
        --refPointLabel 'CTCF motif' \
        --yMin 0 \
        -z '' \
        --colors darkblue \
        --samplesLabel "${SAMPLE_NAME#DNaseC_}"
    
    cd "${BASE_DIR}"
    echo "[$(date)] asChIP-seq analysis completed"
}

run_rm() {
    echo "[$(date)] Cleaning up intermediate files"
    
    rm -f "${OUTPUT_DIR}/02_trim_adapt"/*.fq.gz
    rm -f "${OUTPUT_DIR}/03_merge"/*.fastq.gz
    rm -f "${OUTPUT_DIR}/04_trim_linker"/*.fq.gz
    
    rm -rf "${OUTPUT_DIR}/05_hicpro/data"
    rm -rf "${OUTPUT_DIR}/05_hicpro/hicpro_results/bowtie_results/bwt2_global"
    
    rm -f "${OUTPUT_DIR}/05_hicpro/hicpro_results/hic_results/data/${SAMPLE_NAME}_trimLk3/${SAMPLE_NAME}__hg38XY.bwt2pairs.validPairs"
    
    rm -f "${OUTPUT_DIR}/06_optical_dup"/*.bam
    
    rm -f "${OUTPUT_DIR}/07_fragment_len"/*.txt
    rm -f "${OUTPUT_DIR}/07_fragment_len"/*.singlepair
    
    rm -f "${OUTPUT_DIR}/08_dimerization/${SAMPLE_NAME}"/*.txt
    
    rm -f "${OUTPUT_DIR}/09_aschip/plot_Average_Profile"/*.gz
    
    echo "[$(date)] Cleanup completed"
}

# 主流程
main() {
    init_directories
    
    # 顺序执行各步骤（可根据需要注释掉某些步骤）
    run_fastqc "$READ1" "$READ2"
    run_trim_adapt
    run_merge  # 如果需要 FLASH 合并，取消注释
    run_trim_linker
    run_hicpro
    run_optical_dup
    run_fragment_length
    run_dimerization
    run_aschip
    run_rm  # 可选：清理中间文件
    
    echo "========================================"
    echo "Pipeline completed successfully: $(date)"
    echo "Output directory: ${OUTPUT_DIR}"
    echo "========================================"
}

# 运行
main "$@"