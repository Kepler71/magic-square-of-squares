# Разведка (не независимая проверка): критерий Fable, подставленный в (m^2, n^2), против ранга фактора B у G1.
import sys, json, subprocess
sys.path.insert(0, "linear_sections_codex/fable")
from criterion import criterion
from math import gcd
pairs=[(m,n) for m in range(2,21) for n in range(1,m) if gcd(m,n)==1]
out={}
for m,n in pairs:
    a,b=min(m*m,n*n),max(m*m,n*n)
    out[f"{m},{n}"]=criterion(a,b)['pure_rank0']
json.dump(out,open("unify_claude/crit_g1.json","w"))
print(sum(out.values()), "пар из", len(pairs), "критерий объявляет чистым рангом 0")
