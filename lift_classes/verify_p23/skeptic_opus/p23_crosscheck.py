# Перекрёстная проверка p23_full.py (скептик Claude/Opus, 2026-09-26).
# (1) Вторая реализация классификатора вычетов на чистом Python; сравнение счётчиков с numpy-JSON.
# (2) Корректность статусов: для КАЖДОГО класса вычетов - случайные точные рациональные подъёмы
#     (Fraction), прямая проверка квадратности в Q_p. EXCL обязан провалить тест 8 произведений;
#     PASS обязан его пройти, дать 9 квадратных клеток и v(b),v(c) >= 3 (p=2) / >= 1 (p=3).
# (3) Независимый случайный тест в форме b=rz, c=sz (r,s взаимно простые, r!=+-s).
import sys, json, random, time
from fractions import Fraction as F
from math import gcd
random.seed(20260926)
CELLS=[(i,j) for i in (-1,0,1) for j in (-1,0,1)]
LINES=[[(i,-1),(i,0),(i,1)] for i in (-1,0,1)]+[[(-1,j),(0,j),(1,j)] for j in (-1,0,1)]
LINES+=[[(-1,-1),(0,0),(1,1)],[(-1,1),(0,0),(1,-1)]]

def vint(n,p):
    v=0
    while n%p==0: n//=p; v+=1
    return v,n
def vfrac(x,p):
    a,va=vint(x.numerator,p)[0],None
    va,na=vint(x.numerator,p); vb,nb=vint(x.denominator,p)
    return va-vb,na,nb
def qp_square(x,p):
    assert x!=0
    v,a,b=vfrac(x,p)
    if v%2: return False
    if p==2: return (a*b)%8==1          # b^-1 = b mod 8 для нечётного b
    return (a*b)%3==1                    # b^-1 = b mod 3
def eight_ok(b,c,p):
    val={(i,j):1+i*b+j*c for (i,j) in CELLS}
    if any(x==0 for x in val.values()): return None
    for ln in LINES:
        pr=F(1)
        for cc in ln: pr*=val[cc]
        if not qp_square(pr,p): return False
    return True
def nine_ok(b,c,p):
    return all(qp_square(1+i*b+j*c,p) for (i,j) in CELLS)
def vval(x,p): return 10**9 if x==0 else vfrac(x,p)[0]

# --- (1) классификатор вычетов, чистый Python ---
def classify(p,k,m,B,C):
    mod=p**k; need=3 if p==2 else 1; um_mod=8 if p==2 else 3
    info={}
    for (i,j) in CELLS:
        if (i,j)==(0,0): info[(i,j)]=(False,m,1,10**6); continue
        N=(pow(p,m,mod)+i*B+j*C)%mod
        if N==0: info[(i,j)]=(True,None,None,None); continue
        v,u=vint(N,p); info[(i,j)]=(False,v,u%um_mod,k-v)
    excl=False; allsq=True
    for ln in LINES:
        if any(info[c][0] for c in ln): allsq=False; continue
        V=m+sum(info[c][1] for c in ln)
        U=1
        for c in ln: U=(U*info[c][2])%um_mod
        P=min(info[c][3] for c in ln)
        known=P>=need
        if V%2==1 or (known and U!=1): excl=True
        if not(V%2==0 and known and U==1): allsq=False
    return 'EXCL' if excl else ('PASS' if allsq else 'UND')

def rand_lift(R,p,k):
    # рациональный x in Z_p, x = R mod p^k: x = R + p^k * t/q, q взаимно просто с p
    t=random.randint(-10**6,10**6); q=random.randint(1,10**4)
    while q%p==0: q=random.randint(1,10**4)
    if random.random()<0.3: q=1
    return F(R)+F(p**k*t,q)

t0=time.time()
report={'levels':{},'soundness_violations':0,'lifts_tested':0,'lifts_zero_cell_skipped':0}
LEVELS={2:[(k,m) for k in range(1,7) for m in range(0,8)], 3:[(k,m) for k in range(1,5) for m in range(0,6)]}
num={p:{(r['k'],r['m']):r for r in json.load(open(f'p23_full_p{p}.json'))} for p in (2,3)}
NLIFT=int(sys.argv[1]) if len(sys.argv)>1 else 6
for p in (2,3):
    tgt=3 if p==2 else 1
    for (k,m) in LEVELS[p]:
        mod=p**k; cnt={'EXCL':0,'PASS':0,'UND':0}
        for B in range(mod):
            for C in range(mod):
                if m>0 and B%p==0 and C%p==0: continue
                st=classify(p,k,m,B,C); cnt[st]+=1
                if st=='UND': continue
                for _ in range(NLIFT):
                    b=rand_lift(B,p,k)/p**m; c=rand_lift(C,p,k)/p**m
                    res=eight_ok(b,c,p)
                    if res is None: report['lifts_zero_cell_skipped']+=1; continue
                    report['lifts_tested']+=1
                    if st=='EXCL' and res: report['soundness_violations']+=1; print('VIOL EXCL',p,k,m,B,C,b,c)
                    if st=='PASS':
                        if not (res and nine_ok(b,c,p) and vval(b,p)>=tgt and vval(c,p)>=tgt):
                            report['soundness_violations']+=1; print('VIOL PASS',p,k,m,B,C,b,c)
        n=num[p][(k,m)]
        same=(cnt['EXCL']==n['EXCL'] and cnt['PASS']==n['PASS'] and cnt['UND']==n['UND'])
        report['levels'][f'p{p}_k{k}_m{m}']=dict(cnt,numpy_match=same)
        if not same: print('MISMATCH with numpy',p,k,m,cnt,n)
    print(f'p={p}: уровни пройдены [{time.time()-t0:.1f}s], подъёмов проверено {report["lifts_tested"]}, '
          f'нарушений {report["soundness_violations"]}',flush=True)
print('все счётчики совпали с numpy:',all(v['numpy_match'] for v in report['levels'].values()))

# --- (3) случайный тест в исходной форме b=rz, c=sz ---
for p in (2,3):
    tgt=3 if p==2 else 1
    tested=passes=bad=0
    for it in range(300000):
        r=random.randint(-60,60); s=random.randint(-60,60)
        if r==0 or s==0 or abs(r)==abs(s) or gcd(r,s)!=1: continue
        e=random.randint(-4,6)
        u=F(random.randint(-500,500),random.randint(1,500))
        if u==0 or vval(u,p)!=0: continue
        z=u*F(p)**e
        res=eight_ok(r*z,s*z,p)
        if res is None: continue
        tested+=1
        if res:
            passes+=1
            if not(vval(z,p)>=tgt and nine_ok(r*z,s*z,p)): bad+=1; print('COUNTEREXAMPLE',p,r,s,z)
    report[f'rsz_p{p}']=dict(tested=tested,passes=passes,bad=bad)
    print(f'p={p}: форма rz,sz: проверено {tested}, прошли 8 условий {passes}, контрпримеров {bad}',flush=True)
json.dump(report,open('p23_crosscheck.json','w'),indent=1)
print('ИТОГ: нарушений корректности статусов =',report['soundness_violations'],f'[{time.time()-t0:.1f}s]')
