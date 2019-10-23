#!/bin/bash
cd "$(dirname "$0")"
cd ..

mkdir -p log

#python -u core/ochem.py
thread_num=3
for property_id in 46  114  221  377  375  218  697  156  114  206 ;
do

  bin/main.sh $property_id

done






