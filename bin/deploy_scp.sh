#!/usr/bin/env bash
cd "$(dirname "$0")"
cd ..

work_dir="$(basename "$(pwd)")"

ssh -t -oPort=39121 felix@c23a473940.51mypc.cn "mkdir -p $work_dir"

for fold in bin core notebook;
do
  scp -r  -oPort=39121   ../$work_dir/$fold       felix@c23a473940.51mypc.cn:~/$work_dir/$fold
done

scp -r  -oPort=39121   ../$work_dir/*.*        felix@c23a473940.51mypc.cn:~/$work_dir/*.*
