#!/bin/bash
cd "$(dirname "$0")"
cd ..

mkdir -p log

#python -u core/ochem.py
thread_num=3
for xx in {0..3}:
do
  for property_id in {1..5};
  do
    bin/main.sh $property_id
  done
done






