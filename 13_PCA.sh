#!/bin/bash

# 26/06/05
for i in DMSO4h dTAG4h; do
    for j in R1 R2; do
        zcat /mnt/disk5/1/DNaseC/241112_DNaseC_RC1A3_R1,R2_dATG4h_5per3T3_FA_0.5ulMNase_20m_S-link/fragment_Length_R1R2r0r3r4/DNaseC_RC1A3_${i}_UMI_rmdup_${j}r0r3r4.pair.gz | wc -l
    done
done
# 182118113
# 178177294
# 196530715
# 185237429

# RC1A3
for i in DMSO4h dTAG4h; do
    for j in R1 R2; do
        zcat /mnt/disk5/1/DNaseC/241112_DNaseC_RC1A3_R1,R2_dATG4h_5per3T3_FA_0.5ulMNase_20m_S-link/fragment_Length_R1R2r0r3r4/DNaseC_RC1A3_${i}_UMI_rmdup_${j}r0r3r4.pair.gz | awk '{if($10=="+"){if($11=="+"){print $1,$2,$4,$5,$10,$11}else{print $1,$2,$4,$6,$10,$11}}else{if($11=="+"){print $1,$3,$4,$5,$10,$11}else{print $1,$3,$4,$6,$10,$11}}}' OFS='\t' | awk '{if($1==$3 && $4-$2>0 || $1!=$4){print $0}}' > DNaseC_RC1A3_${i}_${j}.pair
        shuf -n 178170000 DNaseC_RC1A3_${i}_${j}.pair > DNaseC_RC1A3_${i}_${j}_shuf.pair
        cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 --assembly hg38 /ssd/genome/hg38_chromsize_min.txt:1000 DNaseC_RC1A3_${i}_${j}_shuf.pair DNaseC_RC1A3_${i}_${j}_1kb.cool
        cooler zoomify -p 80 DNaseC_RC1A3_${i}_${j}_1kb.cool -r 1000,5000,10000 -o DNaseC_RC1A3_${i}_${j}.mcool --balance
    done
done


awk '{print $1,$2,$3,$9,$10,$11,NR}' OFS='\t' /mnt/disk5/1/DNaseC/total/293T/regular_file/CTCF_loop/DNaseC_293T_loop_groupby.txt > CTCF_loop.bedpe

for i in DMSO4h dTAG4h; do
    for j in R1 R2; do
        chromosight quantify --pattern loops \
                         --subsample 178170000 \
                         --win-fmt npy \
                         CTCF_loop.bedpe \
                         DNaseC_RC1A3_${i}_${j}.mcool::/resolutions/5000 \
                         ./${i}_${j} \
                         --perc-zero 100 \
                         --perc-undetected 100 \
                         -t 80
                     
        chromosight quantify --pattern loops \
                         --subsample 178170000 \
                         --win-fmt npy \
                         CTCF_loop.bedpe \
                         DNaseC_RC1A3_${i}_${j}.mcool::/resolutions/5000 \
                         ./${i}_${j} \
                         --perc-zero 100 \
                         --perc-undetected 100 \
                         -t 80
    done
done


for i in DMSO4h dTAG4h; do
    for j in R1 R2; do
        awk 'NR>1{print $7"_"$8"\t"$9}' ${i}_${j}.tsv | sort -k1,1 > t_${i}_${j}.tab
    done
done

join -t $'\t' t_DMSO4h_R1.tab t_DMSO4h_R2.tab | join -t $'\t' - t_dTAG4h_R1.tab | join -t $'\t' - t_dTAG4h_R2.tab | cut -f2-5 | awk '{if(NF==4){print $0}}' > DMSO4h_dTAG4h_R1R2.tab

rm t_DMSO4h_R1.tab t_DMSO4h_R2.tab t_dTAG4h_R1.tab t_dTAG4h_R2.tab