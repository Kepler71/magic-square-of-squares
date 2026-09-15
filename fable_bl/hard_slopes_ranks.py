# Fable, 15.09.2026. Для незакрытых наклонов (s ≤ 500) и трёх QC-наклонов: ранги (PARI ellrank), знаки функциональных
# уравнений и кручение всех 22 генерических эллиптических множителей и 8 множителей E± шести клеток.
# Проверка гипотезы «всегда найдётся эллиптический множитель ранга 0» (заведомо ложна на открытых) и
# «всегда найдётся биэллиптическая шестиклеточная кривая с rank E+ + rank E- ≤ 1».
from sage.all import *
import json, itertools, sys
sys.path.insert(0,'/home/kep/magicKube/census_genus1')
from g1census import weier
ROOT='/home/kep/magicKube/'
rows=json.load(open(ROOT+'fable_bl/census_rows.json'))
slopes=[r['slope'] for r in rows if r['method']=='open']+['265/298','73/362','126/451']
# 22 подмножества: берём из classification в census_pattern (пересчёт j) — проще: все 3-4 подмножества, но помечаем класс по j
K=FunctionField(QQ,'k'); kk=K.gen()
LAM={'1':K(1),'-1':K(-1),'k':kk,'-k':-kk,'1+k':1+kk,'-1-k':-1-kk,'1-k':1-kk,'k-1':kk-1}
names=list(LAM)
def Emn_j():
    J=[]
    for m,n in [(kk,1),(kk,1-kk),(1,1-kk),(1,1+kk),(kk,1+kk),(1-kk,1+kk),(1,2*kk-1),(1,2*kk+1),(kk,2-kk),(kk,2+kk),(1,2*kk),(kk,2)]:
        J.append(EllipticCurve(K,[0,m**2+n**2,0,m**2*n**2,0]).j_invariant())
    return J
JM=Emn_j()
def generic_curve(S):
    R=PolynomialRing(K,'z'); z=R.gen(); f=R.prod(1+LAM[s]*z for s in S); co=f.list()
    if len(S)==3:
        return EllipticCurve(K,[0,co[2],0,co[1]*co[3],co[0]*co[3]**2])
    a,b,c,d,e=co[4],co[3],co[2],co[1],co[0]
    I=12*a*e-3*b*d+c**2; Jj=72*a*c*e+9*b*c*d-27*a*d**2-27*e*b**2-2*c**3
    return EllipticCurve(K,[0,0,0,-27*I,-27*Jj])
gen22=[]
for r_ in (3,4):
    for S in itertools.combinations(names,r_):
        j=generic_curve(S).j_invariant()
        if j not in QQ and any(j==x for x in JM): gen22.append(S)
print("генерических подмножеств:",len(gen22),flush=True)
def e_pm(a,b,c,kind):
    X=PolynomialRing(QQ,'X').gen()
    pol=(1-QQ(a)**2*X)*(1-QQ(b)**2*X)*(1-QQ(c)**2*X) if kind=='+' else (X-QQ(a)**2)*(X-QQ(b)**2)*(X-QQ(c)**2)
    lc=pol.leading_coefficient(); g=(pol(X/lc)*lc**2).monic() if lc!=1 else pol
    co=g.list(); return EllipticCurve([0,co[2],0,co[1],co[0]])
out={}
for sl in slopes:
    r,s=map(int,sl.split('/')); val={'1':s,'-1':-s,'k':r,'-k':-r,'1+k':s+r,'-1-k':-(s+r),'1-k':s-r,'k-1':r-s}
    rec={'ell':[],'six':[]}
    for S in gen22:
        E,L,c0=weier([QQ(val[x]) for x in S])
        lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]
        w=int(E.root_number())
        rec['ell'].append([list(S),lo,hi,w])
    vals=sorted(set([r,s-r,s,s+r]))
    for a,b,c in itertools.combinations(vals,3):
        for kind in '+-':
            E=e_pm(a,b,c,kind); lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]; w=int(E.root_number())
            rec['six'].append([[a,b,c],kind,lo,hi,w,c==a+b])
    minell=min(x[1] for x in rec['ell']); minell_hi=min(x[2] for x in rec['ell'])
    # биэллиптические шестиклеточные: rank J = rank E+ + rank E- (нижние границы)
    biel=[]
    for a,b,c in itertools.combinations(vals,3):
        lo=sum(x[2] for x in rec['six'] if x[0]==[a,b,c]); lo_lower=sum(x[2] for x in rec['six'] if x[0]==[a,b,c])
        loL=sum(x[2] for x in rec['six'] if x[0]==[a,b,c])
        lower=sum(x[2] for x in rec['six'] if x[0]==[a,b,c])
        biel.append([[a,b,c],sum(x[2] for x in rec['six'] if x[0]==[a,b,c]),sum(x[3] for x in rec['six'] if x[0]==[a,b,c])])
    rec['min_ell_lower']=minell; rec['min_ell_upper']=minell_hi
    rec['bielliptic_rank_lower_upper']=biel
    rec['min_biel_lower']=min(x[1] for x in biel)
    rec['w_plus_count']=sum(1 for x in rec['ell'] if x[3]==1)
    out[sl]=rec
    print(sl,"мин. нижн. ранг эллипт.:",minell,"| мин. верхн.:",minell_hi,"| w=+1 у",rec['w_plus_count'],"из 22 | биэллиптические (нижн,верхн) ранги J:",[(x[0],x[1],x[2]) for x in biel],flush=True)
json.dump(out,open(ROOT+'fable_bl/hard_slopes_ranks.json','w'),indent=1)
print("готово")
