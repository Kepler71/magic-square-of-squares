# Fable, 15.09.2026. Единый разбор переписи наклонов: что общего во всех закрытиях.
# Для каждого наклона r/s (s ≤ 500): (i) какие из 12 классов E_{m,n} проходят точный 2-изогенный Селмер (свой пересчёт,
# решатель criterion_proof/selmer_exact.py Claude); (ii) фактический закрывающий метод по файлам переписи;
# (iii) для закрытий кривыми рода 1 из 3–4 клеток — принадлежит ли подмножество к 22 «генерическим» (j над Q(k) совпадает
# с j одного из 12 классов) или это «случайное» закрытие; (iv) арифметика: число простых ≡1 (mod 4) в N=rs(s−r)(s+r)(2r±s)(2s±r).
import sys, json, os
from math import gcd, isqrt
from collections import Counter, defaultdict
sys.path.insert(0,'/home/kep/magicKube/criterion_proof')
from selmer_exact import selmer
ROOT='/home/kep/magicKube/'
def red(a,b):
    g=gcd(a,b); a//=g; b//=g; return (min(a,b),max(a,b))
PAIRS=lambda r,s: [(r,s),(r,s-r),(s,s-r),(s,r+s),(r,r+s),(s-r,r+s),(s,abs(2*r-s)),(s,2*r+s),(r,2*s-r),(r,2*s+r),(s,2*r),(r,2*s)]
PAIRNAMES=['(r,s)','(r,s-r)','(s,s-r)','(s,r+s)','(r,r+s)','(s-r,r+s)','(s,|2r-s|)','(s,2r+s)','(r,2s-r)','(r,2s+r)','(s,2r)','(r,2s)']
cache={}
def rank0(m,n):
    if (m,n) not in cache:
        S,T=selmer(m,n); cache[(m,n)]=(len(S),len(T))
    a,b=cache[(m,n)]; return a*b<=4
def sq(x): return x>=0 and isqrt(x)**2==x
def primes_of(x):
    x=abs(x); out=[]; d=2
    while d*d<=x:
        if x%d==0:
            out.append(d)
            while x%d==0: x//=d
        d+=1 if d==2 else 2
    if x>1: out.append(x)
    return out
# --- фактические закрытия по файлам
method={}; detail={}
fs=json.load(open(ROOT+'census_genus1/final_status_200.json'))['sources']
for sl,src in fs.items(): method[sl]=src
for l in open(ROOT+'census_six_cells/census.jsonl'):
    d=json.loads(l)
    if d['closed']: detail[d['slope']]=('six',d['by'])
for l in open(ROOT+'census_six_cells/kolyvagin.jsonl'):
    d=json.loads(l)
    if d.get('closed'): detail[d['slope']]=('sixL',d.get('by'))
for l in open(ROOT+'census_genus1/g1census.jsonl'):
    d=json.loads(l)
    if d.get('closed'): detail[d['slope']]=('g1',d['by'])
for f in ('magma_calc/stage2.jsonl','magma_calc/last2.jsonl','magma_calc/extra_48.jsonl','magma_calc/chab_all.jsonl'):
    if not os.path.exists(ROOT+f): continue
    for l in open(ROOT+f):
        d=json.loads(l)
        if d.get('closed') and d['slope'] not in detail: detail[d['slope']]=('g2',{k:v for k,v in d.items() if k in ('c','pair','rb','curve','singles','pairs')})
for f in ('census_300/c201_300.jsonl','census_300/c301_400.jsonl','census_300/c401_500.jsonl'):
    for l in open(ROOT+f):
        d=json.loads(l)
        if d.get('closed'): method[d['slope']]='300:'+d['by'][0]; detail[d['slope']]=('300',d['by'])
        else: method.setdefault(d['slope'],'open')
for f in ('census_300/lpass_201_300.json','census_300/lpass_301_400.json','census_300/lpass_401_500.json'):
    if os.path.exists(ROOT+f):
        for d in json.load(open(ROOT+f)):
            if d.get('closed'): method[d['slope']]='300:L'; detail[d['slope']]=('300L',d['by'])
# --- классификация подмножеств по j над Q(k): считаем j для всех 3- и 4-подмножеств символьно
from sage.all import QQ, FunctionField, PolynomialRing, EllipticCurve
K=FunctionField(QQ,'k'); k=K.gen()
LAM={'1':K(1),'-1':K(-1),'k':k,'-k':-k,'1+k':1+k,'-1-k':-1-k,'1-k':1-k,'k-1':k-1}
lamnames=list(LAM)
def curve_from(S):
    # кривая y^2 = prod (1+λ z) (3 или 4 клетки), через точку z=0 -> Вейерштрасс; для 3: кубика; для 4: квартика с рац. точкой
    R=PolynomialRing(K,'z'); z=R.gen()
    f=R.prod(1+LAM[s]*z for s in S)
    if len(S)==3:
        # y^2=c(z-e1)(z-e2)(z-e3): x=c z -> y^2/c^2... используем EllipticCurve из кубики: Sage умеет Jacobian кубики
        co=f.list(); c3,c2,c1,c0=co[3],co[2],co[1],co[0]
        return EllipticCurve(K,[0,c2,0,c1*c3,c0*c3**2])
    else:
        pass
        # квартика y^2=f(z), f(0)=1: Jacobian через формулу инвариантов (I,J): E: y^2=x^3-27I x-27J
        co=f.list()+[0]*(5-len(f.list()))
        a,b,c,d,e=co[4],co[3],co[2],co[1],co[0]
        I=12*a*e-3*b*d+c**2
        J=72*a*c*e+9*b*c*d-27*a*d**2-27*e*b**2-2*c**3
        return EllipticCurve(K,[0,0,0,-27*I,-27*J])
import itertools
jclass={}
Jmn={}
for i,(m,n) in enumerate([(k,1),(k,1-k),(1,1-k),(1,1+k),(k,1+k),(1-k,1+k),(1,2*k-1),(1,2*k+1),(k,2-k),(k,2+k),(1,2*k),(k,2)]):
    E=EllipticCurve(K,[0,m**2+n**2,0,m**2*n**2,0]); Jmn[i]=E.j_invariant()
subsets=[]
for r_ in (3,4):
    for S in itertools.combinations(lamnames,r_):
        try:
            E=curve_from(S); j=E.j_invariant()
        except Exception as ex:
            jclass[S]=('err',str(ex)[:40]); continue
        if j in QQ: jclass[S]=('const',None); continue
        hit=[i for i in range(12) if Jmn[i]==j]
        jclass[S]=('E_mn',hit[0]) if hit else ('other',None)
gen22=[S for S,v in jclass.items() if v[0]=='E_mn']
print("подмножеств 3–4 клеток с j одного из 12 классов:",len(gen22),"; прочих:",Counter(v[0] for v in jclass.values()))
def signed_subset_to_names(S,r,s):
    vals={r:'k',-r:'-k',s:'1',-s:'-1',s+r:'1+k',-(s+r):'-1-k',s-r:'1-k',-(s-r):'k-1'}
    return tuple(vals[c] for c in S)
# --- главный цикл
def passes_of(rs):
    r,s=rs; passes=[]
    for i,(a,b) in enumerate(PAIRS(r,s)):
        if a<=0 or a==b: continue
        m,n=red(a,b)
        if rank0(m,n) and not (sq(m) and sq(n) and sq(isqrt(m)**2+isqrt(n)**2)): passes.append(i)
    return passes
import multiprocessing as mp
allrs=[(r,s) for s in range(2,501) for r in range(1,s) if gcd(r,s)==1 and 2*r!=s]
need=[(r,s) for r,s in allrs if s<=200]
if os.path.exists(ROOT+'fable_bl/census_rows.json'):
    PASS={}
    for d in json.load(open(ROOT+'fable_bl/census_rows.json')):
        if d['s']<=200: PASS[(d['r'],d['s'])]=d['passes']
else:
    with mp.get_context('fork').Pool(6) as pool:
        PASS=dict(zip(need,pool.map(passes_of,need,chunksize=200)))
rows=[]
for s in range(2,501):
    for r in range(1,s):
        if gcd(r,s)!=1 or 2*r==s: continue
        sl=f"{r}/{s}"
        passes=PASS.get((r,s),[-1])
        N=r*s*(s-r)*(s+r)*(2*r-s)*(2*r+s)*(2*s-r)*(2*s+r)
        P=primes_of(N)
        w1=sum(1 for p in P if p%4==1); w3=sum(1 for p in P if p%4==3)
        P4=primes_of(r*s*(s-r)*(s+r)); w1_4=sum(1 for p in P4 if p%4==1)
        cls=None
        det=detail.get(sl)
        if det and det[0] in ('g1','300','300L'):
            S=det[1][1] if det[0]=='300' else (det[1][0] if det[0]=='g1' else det[1][1])
            if det[0]=='300' and det[1][0]!='g1': cls='six:'+str(det[1][1:])
            else:
                names=signed_subset_to_names(S,r,s)
                # нормализуем по знаку z: подмножество и его отрицание — одна кривая
                key=tuple(sorted(names,key=lamnames.index))
                v=jclass.get(key)
                if v is None:
                    neg={'1':'-1','-1':'1','k':'-k','-k':'k','1+k':'-1-k','-1-k':'1+k','1-k':'k-1','k-1':'1-k'}
                    key2=tuple(sorted((neg[x] for x in names),key=lamnames.index)); v=jclass.get(key2)
                cls=('g1:'+v[0]+(':'+PAIRNAMES[v[1]] if v[0]=='E_mn' else '')) if v else 'g1:?'
        if det and det[0] in ('six','sixL'):
            by=det[1] if isinstance(det[1],list) else [det[1]]
            def flat(x):
                if isinstance(x,list): return [y for z in x for y in flat(z)]
                return [x]
            fl=flat(by); ints=[x for x in fl if isinstance(x,int)]; kind=[x for x in fl if x in ('+','-')]
            kind=kind[0] if kind else '?'
            if len(ints)>=3:
                a,b,c=sorted(ints[:3]); cls='six'+kind+(':c=a+b' if c==a+b else ':c≠a+b')
            else: cls='six?:'+str(by)[:60]
        if det and det[0]=='300' and det[1][0]!='g1':
            a,b,c=det[1][1]; kind=det[1][2] if len(det[1])>2 else '?'
            cls='six'+str(kind)+(':c=a+b' if c==a+b else ':c≠a+b')
        if det and det[0]=='g2':
            g=det[1]; cls='g2:rb='+str(g.get('rb','?'))+(':c='+str(g.get('c'))+',pair='+str(g.get('pair')) if 'pair' in g else ':'+str(g.get('curve')))
        rows.append(dict(slope=sl,r=r,s=s,passes=passes,npass=len(passes),w1=w1,w3=w3,w1_4=w1_4,method=method.get(sl,'?'),cls=cls,detail=det))
json.dump(rows,open(ROOT+'fable_bl/census_rows.json','w'),default=str)
# --- сводки
def summary(rows,label):
    print("\n===",label,"наклонов:",len(rows))
    c=Counter(r['npass'] for r in rows); print("число проходящих классов Селмера:",sorted(c.items()))
    zero=[r for r in rows if r['npass']==0]
    print("наклонов без единого класса ранга 0 по Селмеру:",len(zero))
    print("методы у них:",Counter(r['method'] for r in zero).most_common())
    print("классы закрытий (g1 по j над Q(k); six по c=a+b; g2):",Counter(r['cls'] for r in zero if r['cls']).most_common(40))
    # арифметика
    by=defaultdict(lambda:[0,0])
    for r in rows:
        by[r['w1']][0]+=1; by[r['w1']][1]+=(r['npass']>0)
    print("ω1(N) -> (всего, закрыто Селмером):",{w:tuple(v) for w,v in sorted(by.items())})
    by4=defaultdict(lambda:[0,0])
    for r in rows:
        by4[r['w1_4']][0]+=1; by4[r['w1_4']][1]+=(r['npass']>0)
    print("ω1(rs(s²−r²)) -> (всего, закрыто Селмером):",{w:tuple(v) for w,v in sorted(by4.items())})
    # какие классы срабатывают чаще всего
    cc=Counter(i for r in rows for i in r['passes']); print("частота классов:",{PAIRNAMES[i]:cc[i] for i in range(12)})
    # минимальный ω1 среди наклонов, потребовавших рода 2 / L-значения
    hard=[r for r in rows if 'род 2' in str(r['method']) or 'L' in str(r['method'])]
    print("трудные (род 2 / L(E,1)):",len(hard),"ω1 у них:",Counter(r['w1'] for r in hard))
summary([r for r in rows if r['s']<=200],'s ≤ 200')
rows2=[r for r in rows if 200<r['s']<=500]
print("\n=== 200 < s ≤ 500: наклонов",len(rows2),"методы:",Counter(r['method'] for r in rows2).most_common())
print("классы закрытий (201..500):",Counter(r['cls'] for r in rows2 if r['cls']).most_common(30))
hard=[r for r in rows2 if r['method'] in ('open','300:L') or (r['cls'] and r['cls'].startswith('six'))]
print("трудные при 201..500:",len(hard),"ω1:",Counter(r['w1'] for r in hard))
# наклоны, где ни один из 12 классов не прошёл И закрытие не через класс E_mn (случайные закрытия)
acc=[r for r in rows if r['s']<=200 and r['npass']==0 and r['cls'] and 'other' in r['cls']]
acc2=[r for r in rows if r['s']>200 and r['cls'] and 'other' in r['cls']]
print('при 201..500 закрытий g1 вне 22 классов:',len(acc2))
print("\nзакрытия кривой рода 1 НЕ из 22 генерических классов (случайный ранг 0):",len(acc),[r['slope'] for r in acc][:40])
open_=[r['slope'] for r in rows if r['method']=='open']
print("открытых при s ≤ 500:",len(open_),open_)
