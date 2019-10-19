#!/bin/bash
cd "$(dirname "$0")"
cd ../src

python -u main.py $1

CUDA_VISIBLE_DEVICES=3  &&   python train.py --config resources/train_config_ce.yaml



