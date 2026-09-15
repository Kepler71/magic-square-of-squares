# Fable, 15.09.2026. Геометрические проверки для FABLE_BOMBIERI_LANG_V:
#  (1) D = {x1^2+x2^2=2t^2, y1^2+y2^2=2t^2} ⊂ P^4 — дель Пеццо степени 4 с 4 узлами; V → D — (Z/2)^4-накрытие,
#      ветвление по четырём кривым B_ij: x_i^2+y_j^2=t^2 (класс -2K_D), их роды, попарные касания, счёт 256 ОДТ;
#  (2) слой рода 49 расслоения V → C_1 и его неизотривиальность (j трёхклеточного фактора непостоянен по t);
#  (3) E_{m,n}: y^2=x(x+m^2)(x+n^2) как рациональная эллиптическая поверхность над P^1_{m/n}: слои Кодаиры, ранг над Qbar(t) по Шиоде–Тейту;
#  (4) семейство E^+ шести клеток с c=a+b: тип поверхности (K3), слои, граница ранга над Qbar(a/b).
from sage.all import *
import json, sys
out={}
R=PolynomialRing(QQ,'x1,x2,y1,y2,t'); x1,x2,y1,y2,t=R.gens()
Q1=x1**2+x2**2-2*t**2; Q2=y1**2+y2**2-2*t**2
D=R.ideal([Q1,Q2])
# (1a) особые точки D: якобиан ранга <2
J=jacobian([Q1,Q2],[x1,x2,y1,y2,t])
sing=D+R.ideal(J.minors(2))
sing_rad=sing.radical()
pd=sing_rad.primary_decomposition()
out['D_singular_components']=[str(p.gens()) for p in pd]
out['D_singular_degree']=int(sing_rad.hilbert_polynomial())  # число точек над Qbar (сумма степеней нульмерных компонент)
print("особые точки D:",out['D_singular_degree'],[str(p.gens()) for p in pd],flush=True)
# (1b) кривые ветвления B_ij, роды и пересечения
B={}
for i,xi in ((1,x1),(2,x2)):
    for j,yj in ((1,y1),(2,y2)):
        B[(i,j)]=xi**2+yj**2-t**2
gen={}
Ku=FunctionField(QQ,'u'); uu=Ku.gen()
# C_1 ∋ (x1,x2): x1^2+x2^2=2; B_11: x1^2+y1^2=1 параметризуем x1=(1-u^2)/(1+u^2), y1=2u/(1+u^2); x2^2=2-x1^2, y2^2=2-y1^2
X1=(1-uu**2)/(1+uu**2); Y1=2*uu/(1+uu**2)
Pw=PolynomialRing(Ku,'w'); w=Pw.gen()
A=2-X1**2; Bb=2-Y1**2
# примитивный элемент θ=x2+y2: θ^4-2(A+B)θ^2+(A-B)^2=0
L2=Ku.extension(w**4-2*(A+Bb)*w**2+(A-Bb)**2,'th')
gen['B11']=int(L2.genus())
# B_12: x1^2+y2^2=1: y2=2u/(1+u^2), y1^2=2-y2^2 — та же пара классов (A, B): тождественно та же кривая с переименованием y1<->y2
gen['B12']='=B11 по симметрии y1<->y2'
# контроль формулы: y^2=2-X1^2 над Q(u) должен иметь род 0 (коника)
gen['ctrl_conic']=int(Ku.extension(w**2-A,'c').genus())
print("род B_ij:",gen,flush=True)
inter={}
keys=list(B)
for a in range(4):
    for b in range(a+1,4):
        I=D+R.ideal([B[keys[a]],B[keys[b]]])
        Ir=I.radical()
        npts=int(Ir.hilbert_polynomial()); deg=int(I.hilbert_polynomial())
        inter[str((keys[a],keys[b]))]=(npts,deg)
        print("B",keys[a],"∩ B",keys[b],": точек",npts,", с кратностью",deg,flush=True)
out['B_intersections']=inter
# тройные пересечения
trip=0
from itertools import combinations
for c3 in combinations(keys,3):
    I=D+R.ideal([B[k] for k in c3])
    trip+= int(I.hilbert_polynomial()) if I.dimension()>=1 else 0
out['B_triple_points']=trip
# ветвится ли накрытие в узлах D: значения B_ij в узлах
nodes_on_B=0
for p in pd:
    for k,f in B.items():
        if f in p: nodes_on_B+=1
out['nodes_on_branch']=nodes_on_B
# счёт ОДТ V: над каждой точкой касания двух B (кратность 2, третья не проходит) — (Z/2)^4/(Z/2)^2=4 ОДТ; над узлом D — 16
tang=sum(v[0] for v in inter.values())
out['ODP_count']=4*tang+16*out['D_singular_degree']
print("точек касания:",tang,"; ОДТ =",out['ODP_count'],flush=True)

# (2) слой рода 49: t фиксировано, u^2=1+t, v^2=1-t; клетки 1±q, u^2±q, v^2±q. Род через Риман–Гурвиц: (Z/2)^6, 7 точек ветвления -> 49.
# Проверка прямым счётом рода гиперэллиптических факторов (Кани–Розен): сумма floor((|T|-1)/2) по T ⊂ 6 = 20*1+15*1+6*2+1*2=49
kr=sum(binomial(6,k)*((k-1)//2) for k in range(1,7))
out['genus49_KaniRosen']=kr
# неизотривиальность: j трёхклеточного фактора y^2=(1-q^2)(u^2+q) как функция u (v не входит)
K=FunctionField(QQ,'u'); u=K.gen()
S=PolynomialRing(K,'q'); q=S.gen()
E=EllipticCurve(K,[0,-1,0,-u**2,u**2]) if False else None
# y^2=(1-q)(1+q)(u^2+q) = -q^3 - u^2 q^2 + q + u^2 -> умножим на -1: (-y)^2... возьмём Y^2 = -(...): чтобы monic, положим q=-X: Y^2 = X^3 - u^2 X^2 - X + u^2
E=EllipticCurve(K,[0,-u**2,0,-1,u**2])
jE=E.j_invariant()
out['fiber_factor_j_constant']=bool(jE in QQ)
print("j трёхклеточного фактора слоя постоянен?",out['fiber_factor_j_constant'],flush=True)

# (3) E_{m,n} над Q(t), t=m/n: y^2=x(x+t^2)(x+1)
E3=EllipticCurve(K,[0,u**2+1,0,u**2,0])
kod={}
for pl in E3.discriminant().numerator().factor():
    pass
# Sage: local data по местам функционального поля
try:
    Emin=E3
    disc=E3.discriminant()
    fac=disc.numerator().factor()
    kd=[]
    tot=0
    for f,e in fac:
        kd.append((str(f),int(e)))
    # степень дискриминанта и слой на бесконечности
    degdisc=disc.numerator().degree()-disc.denominator().degree()
    kd.append(('inf',12-degdisc if degdisc<12 else 0))
    kod['Emn_disc_valuations']=kd
    print("E_{m,n}: нули Δ и кратности",kd,flush=True)
except Exception as e:
    kod['err']=str(e)
# все слои мультипликативные? проверим c4 в нулях Δ: если c4 не делится, слой I_n
c4=E3.c4()
kod['Emn_c4']=str(c4.factor())
out['Emn']=kod
# Шиода–Тейт для рациональной поверхности: rank = 8 - sum(m_v-1); при I_4,I_4,I_2,I_2: 8-(3+3+1+1)=0
# (4) E^+ с c=a+b, a/b=u: v^2=(1-a^2 X)(1-b^2 X)(1-c^2 X); в весах: замена X -> дает кубику; возьмём Y^2=(X-a^2)(X-b^2)(X-c^2) (E^-) и E^+
a=u; b=K(1); c=a+b
Eminus=EllipticCurve(K,[0,-(a**2+b**2+c**2),0,a**2*b**2+a**2*c**2+b**2*c**2,-(a*b*c)**2])
Eplus=EllipticCurve(K,[0,-(a**2*b**2+a**2*c**2+b**2*c**2),0,(a**2+b**2+c**2)*(a*b*c)**2,-(a*b*c)**4]) # v^2 = prod(1-λ^2 X): X=1/W -> W^3 v^2 = prod(W-λ^2); замена
for name,EE in (('Eplus',Eplus),('Eminus',Eminus)):
    d=EE.discriminant(); fac=d.numerator().factor(); degd=d.numerator().degree()-d.denominator().degree()
    info=[(str(f),int(e)) for f,e in fac]
    c4v=EE.c4()
    info_c4=str(c4v.numerator().factor())
    # минимизируем по каждому месту: Sage умеет local_data? используем E.minimal_model? Для функциональных полей нет; считаем v(Δ) и v(c4)
    out[name]={'disc_factors':info,'deg_disc':degd,'c4_factors':info_c4,'torsion_order_bound':None}
    print(name,"Δ:",info,"deg",degd,"c4:",info_c4,flush=True)
json.dump(out,open('/home/kep/magicKube/fable_bl/geom_check.json','w'),indent=1,default=str)
print("готово")
