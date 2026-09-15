# Fable, 15.09.2026. (a) закрытие 96/407 множителем ранга 0; (b) E8+: v^2=prod_{λ∈{r,s-r,s,s+r}}(1-λ^2 u) — фактор
# восьмиклеточной кривой рода 3 по z↦-z: тип семейства над Q(k), слои, ранги на 49 трудных наклонах; (c) классификация
# 32 закрытий L(E,1) шести клеток по c=a+b.
from sage.all import *
import json, itertools, sys
sys.path.insert(0,'/home/kep/magicKube/census_genus1')
from g1census import weier, zvals, nine_ok
ROOT='/home/kep/magicKube/'
out={}
# (a)
H=json.load(open(ROOT+'fable_bl/hard_slopes_ranks.json'))
for sl in ('96/407','95/349'):
    r,s=map(int,sl.split('/')); val={'1':s,'-1':-s,'k':r,'-k':-r,'1+k':s+r,'-1-k':-(s+r),'1-k':s-r,'k-1':r-s}
    res=[]
    for S,lo,hi,w in H[sl]['ell']:
        if hi==0 or lo==0:
            E,L,c0=weier([QQ(val[x]) for x in S])
            how=None
            if hi==0: how='ellrank'
            else:
                try:
                    if E.minimal_model().lseries().L_ratio()!=0: how='L_ratio'
                except Exception as e: how='L_err:'+str(e)[:40]
            zs=zvals(E,L,c0); bad=[str(z) for z in zs if z!=0 and nine_ok(z,r,s)]
            res.append(dict(S=S,lo=lo,hi=hi,w=w,how=how,z=[str(z) for z in zs],bad=bad,cond=int(E.conductor()) if hi==0 else None))
    out[sl]=res; print(sl,res,flush=True)
# (b) E8+
K=FunctionField(QQ,'k'); k=K.gen()
R=PolynomialRing(K,'u'); u=R.gen()
f=(1-k**2*u)*(1-(1-k)**2*u)*(1-u)*(1-(1+k)**2*u)
co=f.list(); a,b,c,d,e=co[4],co[3],co[2],co[1],co[0]
I=12*a*e-3*b*d+c**2; J=72*a*c*e+9*b*c*d-27*a*d**2-27*e*b**2-2*c**3
E8=EllipticCurve(K,[0,0,0,-27*I,-27*J])
disc=E8.discriminant(); fac=disc.numerator().factor(); degd=disc.numerator().degree()-disc.denominator().degree()
c4=E8.c4()
# минимальная модель по местам: ищем места, где v(Δ)>=12 и v(c4)>=4 — тогда не минимальна; для грубой оценки печатаем
info=dict(disc=[(str(g),int(m)) for g,m in fac],deg_disc=int(degd),c4=str(c4.factor()),j=str(E8.j_invariant().factor()) if E8.j_invariant()!=0 else '0')
out['E8plus_generic']=info; print("E8+ над Q(k):",info,flush=True)
# ранги E8+ на трудных наклонах
def E8_at(r,s):
    X=PolynomialRing(QQ,'X').gen()
    pol=(1-QQ(r)**2*X)*(1-QQ(s-r)**2*X)*(1-QQ(s)**2*X)*(1-QQ(s+r)**2*X)
    co=pol.list(); a,b,c,d,e=co[4],co[3],co[2],co[1],co[0]
    I=12*a*e-3*b*d+c**2; J=72*a*c*e+9*b*c*d-27*a*d**2-27*e*b**2-2*c**3
    return EllipticCurve(QQ,[0,0,0,-27*I,-27*J])
res8={}
for sl in H:
    r,s=map(int,sl.split('/')); E=E8_at(r,s)
    lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]; w=int(E.root_number()); tor=E.torsion_order()
    res8[sl]=(lo,hi,w,int(tor))
print("E8+ на трудных наклонах (lo,hi,w,#tors):",res8,flush=True)
out['E8plus_hard']=res8
print("E8+ с верхней границей 0:",[sl for sl,v in res8.items() if v[1]==0])
# контроль: E8+ на нескольких закрытых Селмером наклонах — ранг обычно?
ctrl={}
for sl in ('1/3','2/5','3/7','1/4','3/8','5/9','4/11','7/12'):
    r,s=map(int,sl.split('/')); E=E8_at(r,s); lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]; ctrl[sl]=(lo,hi,int(E.torsion_order()))
print("контроль E8+ на малых наклонах:",ctrl); out['E8plus_ctrl']=ctrl
# (c)
cnt={'c=a+b':0,'c≠a+b':0}; kinds={}
for l in open(ROOT+'census_six_cells/kolyvagin.jsonl'):
    dd=json.loads(l)
    if dd.get('closed'):
        by=dd['by']; a_,b_,c_=sorted(by['abc']); key=('c=a+b' if c_==a_+b_ else 'c≠a+b'); cnt[key]+=1
        kinds[by['kind']+key]=kinds.get(by['kind']+key,0)+1
print("L(E,1)-закрытия шести клеток:",cnt,kinds); out['sixL']=[cnt,kinds]
json.dump(out,open(ROOT+'fable_bl/extra_checks.json','w'),indent=1,default=str)
print("готово")
