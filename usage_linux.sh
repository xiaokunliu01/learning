#!/bin/bash

# 列转行
awk 'BEGIN{FS=OFS="\t"}{for(i=1;i<NF+1;i++){a[NR,i]=$i}}END{for(i=1;i<NF+1;i++){for(j=1;j<NR;j++){printf a[j,i]"\t"}print a[j,i]}}' file

# 去除bed里面有重叠的行
# 最后一行需自己判断
awk '{
if(NR==1){
all=$0; 
tmp=$3; 
chr=$1; 
s=$5
}
else if(NR>1 && $1==chr){
if($2-tmp<=0){if($5>s){all=$0; tmp=$3; chr=$1; s=$5}} else{print all; all=$0; tmp=$3; chr=$1; s=$5}
}
else if(NR>1 && $1!=chr){
print all; 
all=$0; 
tmp=$3; 
chr=$1; 
s=$5}
}' file