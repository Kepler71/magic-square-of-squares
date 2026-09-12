#!/bin/bash
cd /home/kep/magicKube/family
for sec in "23 7 17" "31 17 25" "41 1 29" "47 23 37" "49 31 41" "73 17 53" "71 49 61" "89 23 65" "79 47 65"; do
  set -- $sec
  echo "######## section ($1,$2,$3)"
  timeout 1200 sage pipeline.sage $1 $2 $3 150 2>&1 | grep -vE 'Deprecation|issues/' | grep -E '^section|\(1\)|\(2\)|\(3\)|class \(|control|\(4\)|UNKNOWN|Error|assert'
done
echo BATCH_DONE
