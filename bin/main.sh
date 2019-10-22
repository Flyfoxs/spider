#!/bin/bash
cd "$(dirname "$0")"
cd ..

mkdir -p log

#python -u core/ochem.py
thread_num=3
python -u core/ochem.py  process_one_item $1 $thread_num > log/p_$1.log 2>&1
python -u core/ochem.py  fill_smiles      $1 $thread_num > log/s_$1.log 2>&1






