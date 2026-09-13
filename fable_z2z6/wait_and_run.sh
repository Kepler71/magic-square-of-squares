#!/bin/bash
cd /home/kep/magicKube/fable_z2z6
while [ $(cat grid/shard_200_*.jsonl | wc -l) -lt 12231 ]; do sleep 15; done
python3 coverage.py > coverage.log 2>&1
python3 explicit_formula.py 200 > explicit_200.log 2>&1
nohup sage sample_pari.sage > sample_pari.log 2>&1 &
echo done
