#!/bin/bash

quast_route="/public/home/hpc214712170/shixf/tools/quast-5.0.2"
cur_dir=$(realpath $(dirname $0))
asm=        #$1  use ref1 as asm
Reference=  #$2
out_dir=    #$3
threads=40    #$4
large=      #$5
# params=     #
## 大于100m选large
if [ "$large" == "no" ];then
  time python $quast_route/quast.py $asm -o $out_dir -t $threads -r $Reference $params
elif [ "$large" == "yes" ];then
  time python $quast_route/quast.py --large $asm -o ${out_dir}_large -t $threads -r $Reference $params
else
  echo Error
fi