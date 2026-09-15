# Fable 15.09. Генерическая структура над ℚ(k): 126 кривых = 4-подмножества сетки Λ₀={0,±1,±k,±(1−k),±(1+k)}
# (3-клеточные T ↔ 4-подмножества с 0: y²=∏(1+λz), λ∈T ⇔ в u=1/z точки ветвления −λ и 0).
# Группируем по j ∈ ℚ(k); внутри группы — класс твиста (c4,c6); между группами — Φ_2(j,j')=0.
from sage.all import *
import itertools, json, sys
sys.path.insert(0,'/home/kep/magicKube/census_genus1'); from g1census import weier
K = FractionField(PolynomialRing(QQ,'k')); k = K.gen()
Lam = [K(1),K(-1),k,-k,1-k,k-1,1+k,-1-k]
def curve(T):   # модель weier над ℚ(k) (та же, что в переписи)
    c0=T[0]; others=T[1:]; R=PolynomialRing(K,'w'); w=R.gen()
    g=R.prod((c0-c)*w + c*c0 for c in others)
    if len(T)==3: g=c0*w*g
    L=g.leading_coefficient(); h=(g(w/L)*L**2).monic(); co=h.list()
    return EllipticCurve(K,[0,co[2],0,co[1],co[0]])
subs=[T for m in (3,4) for T in itertools.combinations(Lam,m)]
assert len(subs)==126
data=[]
for T in subs:
    E=curve(list(T)); data.append((T,E,E.j_invariant()))
groups={}
for T,E,j in data: groups.setdefault(j,[]).append((T,E))
print('различных j над ℚ(k):',len(groups))
# твист: E ≅ E' над ℚ(k) ⇔ (c4'/c4)^3 = (c6'/c6)^2 и d = (c6'/c6)/(c4'/c4) квадрат ⇒ иначе твист на d
def twist_factor(E,F):
    c4,c6=E.c4(),E.c6(); d4,d6=F.c4(),F.c6()
    if c4==0 or c6==0: return None
    d=(d6/c6)/(d4/c4); assert (d4/c4)**3==(d6/c6)**2
    return d
def sqfree(d):  # квадратная часть в ℚ(k): выделяем sqfree по многочленам
    n,dn=d.numerator(),d.denominator(); f=n*dn
    out=K(1); c=f.leading_coefficient(); out*=squarefree_part(QQ(c))
    for p,e in f.factor():
        if e%2: out*=p
    return out
report=[]
for j,lst in groups.items():
    T0,E0=lst[0]
    tw=[str(sqfree(twist_factor(E0,E))) for T,E in lst]
    report.append(dict(j=str(j), n=len(lst), Ts=[[str(t) for t in T] for T,E in lst], twists=tw))
report.sort(key=lambda d:-d['n'])
for d in report: print(d['n'], d['Ts'], 'твисты:', d['twists'])
json.dump(report,open('/home/kep/magicKube/isogeny_mult/fable/generic_j.json','w'),indent=0)
