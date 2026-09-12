#!/bin/bash
cd /home/kep/magicKube/family
for spec in "71 49 61 671 36/25" "89 23 65 910 121/9" "79 47 65 455 64"; do
  set -- $spec
  dk=$(sage -c "A=($2^2+$3^2)/2; C=($1^2+$3^2)/2; n=ZZ((C-A)*(C+A)); print(sign(n)*prod(p^(e%2) for p,e in n.abs().factor()))" 2>/dev/null)
  for p in $(sage -c "print(' '.join(str(p) for p in primes(7,300) if kronecker($dk,p)==1))" 2>/dev/null); do
    out=$(timeout 900 sage chabauty.sage $1 $2 $3 $4 $5 $p 2>&1 | grep -vE 'Deprecation|issues/')
    N=$(echo "$out" | grep -oP 'N = \K\d+' | head -1); [ -z "$N" ] && continue
    if [ $N -le 12 ]; then echo "######## ($1,$2,$3) D=$4 t0=$5 p=$p N=$N"; echo "$out" | grep -E 'local saturation|classes:|rigorous|TOTAL|known P|exact zero|Strassmann bound'; break; fi
  done
done
echo BATCH2_DONE
