import json, time, sys
sys.path.insert(0,'../linear_sections_codex/fable')
from selmer_exact import selmer, DepthExhausted
from criterion import S_phi as F_Sphi, S_phihat as F_Sphihat      # формула Fable — только для сравнения
from math import log2
rows=json.load(open('../linear_sections_codex/fable/grid_rows.json'))
rows=[x for x in rows if x['n']<=int(sys.argv[1]) ] if len(sys.argv)>1 else rows
t0=time.time(); pos_fail=[]; pari_bad=[]; parity_bad=[]; fable_diff=[]; errs=[]; r0=0
for x in rows:
    m,n=x['m'],x['n']
    try: S,T=selmer(m,n)
    except DepthExhausted as e: errs.append((m,n,str(e))); continue
    must={1,-1,m*n,-m*n}
    from math import isqrt
    must={d if isqrt(abs(d))**2!=abs(d) else (1 if d>0 else -1) for d in must}   # классы по модулю квадратов
    def sqfree(d):
        s=1 if d>0 else -1; d=abs(d); out=1; p=2
        while p*p<=d:
            while d%(p*p)==0: d//=p*p
            if d%p==0: out*=p; d//=p
            p+=1
        return s*out*d
    must={sqfree(d) for d in {1,-1,m*n,-m*n}}
    if not must<=set(S) or 1 not in T: pos_fail.append((m,n))
    rb=int(round(log2(len(S)*len(T))))-2
    if rb < x['r1']: pari_bad.append((m,n,rb,x['r1']))
    if (rb - (x['r2']+x['s'])) % 2: parity_bad.append((m,n))
    if rb==0: r0+=1
    if set(S)!=set(F_Sphi(m,n)) or set(T)!=set(F_Sphihat(m,n)): fable_diff.append((m,n))
print(f"пар: {len(rows)}, время {time.time()-t0:.0f} с, исчерпаний глубины {len(errs)}")
print(f"положительный контроль (±1, ±mn ∈ Sel^phi, 1 ∈ Sel^phihat): провалов {len(pos_fail)} {pos_fail[:3]}")
print(f"граница ранга ниже доказанной нижней границы PARI: {len(pari_bad)} {pari_bad[:3]}")
print(f"нечётная разность с Селмером PARI (r2+s): {len(parity_bad)} {parity_bad[:3]}")
print(f"доказан ранг 0 (точный Селмер): {r0}; у PARI ранг 0: {sum(1 for x in rows if x['r2']==0)}")
print(f"расхождений с формулой Fable (множества S^phi, S^phihat): {len(fable_diff)} {fable_diff[:8]}")
json.dump({'fable_diff':fable_diff,'errs':errs},open('selmer_controls.json','w'))
