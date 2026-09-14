"""Fable, 14.09.2026. Девятикратное расширение теоремы о наклонах 1/ell (RESULT_SLOPES_1_OVER_L):
пара (1, ell) входит в список 12 пар (RESULT_ODD_CRITERION, теорема §3) у ДЕВЯТИ наклонов, а не двух.
Свой код: критерий через нечётные простые реализован заново (Sphi, Sphihat по формулировке теоремы §1),
сверка с переписью (coverage_exact.json, census_300/*.jsonl, lpass_*.json, magma_calc/last*.jsonl) — только как контроль.
Отрицательный контроль — квадрат Бремнера–Саллоуса (наклон 247/825, семь квадратов): пары, у которых все три клетки
тройки {a,b,a+b} — квадраты, ОБЯЗАНЫ не проходить критерий (у E_{a,b} есть невырожденная точка).
"""
import json, sys, os
from math import gcd
from fractions import Fraction
ROOT='/home/kep/magicKube'

def primes_of(x):
    x=abs(x); out=[]; d=2
    while d*d<=x:
        if x%d==0:
            out.append(d)
            while x%d==0: x//=d
        d+=1 if d==2 else 2
    if x>1: out.append(x)
    return out
def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
def squarefree_divisors(P):
    out=[1]
    for p in P: out=out+[d*p for d in out]
    return out
def odd_criterion(m,n):
    """Sphi, Sphihat из теоремы §1 RESULT_ODD_CRITERION (моя реализация). Возвращает (|Sphi|,|Sphihat|, граница)."""
    assert gcd(m,n)==1 and m!=n
    Pmn=primes_of(m*n); Pd=primes_of(n*n-m*m)
    q1=[q for q in Pd if q%4==1]
    Sphi=[]
    for d0 in squarefree_divisors(Pmn):
        for d in (d0,-d0):
            if all(legendre(d,q)!=-1 for q in q1): Sphi.append(d)
    Sphihat=[]
    for d in squarefree_divisors(Pd):
        if any(q%2==1 and q%4==3 for q in primes_of(d)): continue
        if any(p%2==1 and legendre(d,p)==-1 for p in Pmn): continue
        Sphihat.append(d)
    from math import log2
    return len(Sphi),len(Sphihat),int(round(log2(len(Sphi)*len(Sphihat))))-2

def twelve_pairs(r,s):
    """12 пар (m,n) теоремы §3 RESULT_ODD_CRITERION для наклона r/s, после сокращения на НОД."""
    raw=[(r,s),(r,s-r),(s,s-r),(s,r+s),(r,r+s),(s-r,r+s),(s,abs(2*r-s)),(s,2*r+s),(r,2*s-r),(r,2*s+r),(s,2*r),(r,2*s)]
    out=[]
    for a,b in raw:
        g=gcd(a,b); a,b=a//g,b//g
        out.append(tuple(sorted((a,b))))
    return out

def ell_condition(l):
    if l%8 not in (3,5): return False
    return all(q%4==3 for q in primes_of(l*l-1) if q%2==1)
def is_prime(n):
    if n<2: return False
    return primes_of(n)==[n]

def nine_slopes(l):
    cands=[Fraction(1,l),Fraction(l-1,l),Fraction(1,l+1),Fraction(1,l-1),Fraction(l-1,l+1),
           Fraction(l-1,2*l),Fraction(l+1,2*l),Fraction(2,l+1),Fraction(2,l-1)]
    out=[]
    for k in cands:
        if 0<k<1 and 2*k!=1 and k not in out: out.append(k)
    return out

# --- перепись: закрыт ли наклон где-либо ---
def load_census():
    closed={}
    cov=json.load(open(f'{ROOT}/criterion_proof/coverage_exact.json'))
    openset=set()
    for rng in cov.values(): openset|=set(rng['open'])
    fs=json.load(open(f'{ROOT}/census_genus1/final_status_200.json'))
    for sl,src in fs['sources'].items(): closed[sl]=src
    for f in ['c201_300','c301_400','c401_500']:
        for line in open(f'{ROOT}/census_300/{f}.jsonl'):
            rec=json.loads(line)
            if rec.get('closed'): closed.setdefault(rec['slope'],'census_300:'+str(rec['by'][0]))
    for f in ['lpass_201_300','lpass_301_400','lpass_401_500']:
        p=f'{ROOT}/census_300/{f}.json'
        if os.path.exists(p):
            for rec in json.load(open(p)):
                if rec.get('closed'): closed.setdefault(rec['slope'],'lpass:L_ratio')
    for f in ['last300','last400','last2','rest3','extra_48','stage2']:
        p=f'{ROOT}/magma_calc/{f}.jsonl'
        if os.path.exists(p):
            for line in open(p):
                try: rec=json.loads(line)
                except Exception: continue
                if 'slope' in rec and ('PROVEN: true' in str(rec.get('out','')) or rec.get('closed')):
                    closed.setdefault(rec['slope'],'magma:'+f)
    return closed, openset

def status(k,closed,openset):
    sl=f'{k.numerator}/{k.denominator}'; s=k.denominator
    if s<=200:
        if sl in closed: return closed[sl]
        if sl not in openset: return 'точный Селмер 12 пар (coverage_exact)'
        return 'НЕ НАЙДЕН'
    if s<=500:
        return closed.get(sl,'НЕ НАЙДЕН (201–500)')
    return 'вне переписи (s>500)'

if __name__=='__main__':
    closed,openset=load_census()
    L=[l for l in range(3,3001) if is_prime(l) and ell_condition(l)]
    print('простые ell<=3000 с условием теоремы:',L)
    total=set(); report=[]
    for l in L:
        sp,sh,b=odd_criterion(1,l); assert sp*sh<=4, (l,sp,sh)
        for k in nine_slopes(l):
            r,s=k.numerator,k.denominator
            assert (1,l) in twelve_pairs(r,s), (l,k,twelve_pairs(r,s))
            total.add(k)
            st=status(k,closed,openset)
            report.append((l,str(k),st))
    print('всего различных наклонов от девятикратного расширения:',len(total))
    from collections import Counter
    c=Counter(st for _,_,st in report); print('статусы в переписи:');
    for k_,v in c.items(): print('  ',v,k_)
    print('первые 30 строк:')
    for row in report[:30]: print('  ',row)
    bad=[r for r in report if 'НЕ НАЙДЕН' in r[2]]
    print('противоречий с переписью (наклон s<=500, не закрытый там):',len(bad), bad[:10])
    # --- отрицательный контроль: Бремнер–Саллоус ---
    print('\n--- отрицательный контроль: квадрат Бремнера–Саллоуса ---')
    sq=[[373**2,289**2,565**2],[360721,425**2,23**2],[205**2,527**2,222121]]
    S=sum(sq[0]); assert all(sum(row)==S for row in sq) and all(sum(sq[i][j] for i in range(3))==S for j in range(3))
    assert sq[0][0]+sq[1][1]+sq[2][2]==S and sq[0][2]+sq[1][1]+sq[2][0]==S
    c0=sq[1][1]; x=sq[0][2]-c0; y=sq[2][2]-c0
    assert sq[1][0]-c0==x+y and sq[2][1]-c0==x-y, (x,y)
    k=Fraction(y,x); r,s=k.numerator,k.denominator
    print('наклон Бремнера–Саллоуса k = y/x =',k,' s-r =',s-r,' s+r =',s+r)
    z=Fraction(x,s*c0)
    cells={c:1+c*z for c in (s,-s,r,-r,s-r,-(s-r),s+r,-(s+r))}
    def is_sq(q):
        from math import isqrt
        return q>0 and isqrt(q.numerator)**2==q.numerator and isqrt(q.denominator)**2==q.denominator
    print('клетки-квадраты по коэффициентам c (1+cz):',{c:is_sq(v) for c,v in cells.items()})
    # тройки {a,b,a+b} c коэффициентами, где все три клетки — квадраты, дают невырожденную точку на E_{|a|,|b|}
    coefs=[s,-s,r,-r,s-r,-(s-r),s+r,-(s+r)]
    triples=set()
    for a in coefs:
        for b in coefs:
            if a+b in coefs and a<b and all(is_sq(cells[c]) for c in (a,b,a+b)):
                m,n=sorted((abs(a),abs(b))); g=gcd(m,n); triples.add((m//g,n//g))
    print('пары (m,n) с невырожденной точкой:',sorted(triples))
    for m,n in sorted(triples):
        sp,sh,b=odd_criterion(m,n)
        print(f'  ({m},{n}): |Sphi|={sp} |Sphihat|={sh} граница ранга {b}  -> критерий {"ПРОШЁЛ (ОШИБКА!)" if sp*sh<=4 else "не прошёл, как и должен"}')
        assert sp*sh>4
    json.dump(dict(primes=L,slopes=sorted(str(k) for k in total),report=report,bs_slope=str(k),bs_pairs=sorted(triples)),open(f'{ROOT}/fable_infinity/ninefold.json','w'),ensure_ascii=False,indent=1)
