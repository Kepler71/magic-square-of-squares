"""Таблица 2-адических условий: для C_d и C'_d проверить, что разрешимость в Q_2 зависит только от
класса d в Q_2^*/Q_2^{*2} и (m,n) mod 2^j; найти минимальное j."""
import json
from collections import defaultdict
from isogeny_descent import qp_soluble,qp_is_square,squarefree_signed_divisors
from math import gcd
def q2class(d):
    e=0
    while d%2==0: d//=2; e+=1
    return (e%2, d%8)
def sol_C(d,m,n):
    A2,B2=m*m,n*n; c=[d*A2*B2,0,d*d*(A2+B2),0,d**3]
    return qp_soluble(c,2,qp_is_square(-d*A2,2) or qp_is_square(-d*B2,2))
def sol_Cp(d,m,n):
    A2,B2=(m+n)**2,(m-n)**2; c=[d*A2*B2,0,-d*d*(A2+B2),0,d**3]
    return qp_soluble(c,2,qp_is_square(d,2))
recC=defaultdict(set); recCp=defaultdict(set)
MOD=64
for n in range(2,130):
    for m in range(1,n):
        if gcd(m,n)!=1: continue
        for d in squarefree_signed_divisors(m*n):
            recC[(q2class(d),m%MOD,n%MOD)].add(sol_C(d,m,n))
        for d in squarefree_signed_divisors((m+n)*(m-n)):
            if d<0: continue
            recCp[(q2class(d),m%MOD,n%MOD)].add(sol_Cp(d,m,n))
def minimal_modulus(rec,name):
    for j in (1,2,4,8,16,32,64):
        agg=defaultdict(set)
        for (cl,a,b),v in rec.items(): agg[(cl,a%j,b%j)]|=v
        amb=sum(1 for v in agg.values() if len(v)>1)
        print(f'{name}: modulus {j}: классов {len(agg)}, неоднозначных {amb}')
        if amb==0: return j,agg
    return None,None
jC,aggC=minimal_modulus(recC,'C_d ')
jCp,aggCp=minimal_modulus(recCp,'C\'_d')
json.dump({'C':{str(k):list(v)[0] for k,v in aggC.items()},'Cp':{str(k):list(v)[0] for k,v in aggCp.items()},'modC':jC,'modCp':jCp},open('local2_table.json','w'))
# компактная печать: для каждого класса d — при каких (m mod j, n mod j) разрешимо
CL={(0,1):'1',(0,3):'3',(0,5):'5',(0,7):'7',(1,1):'2',(1,3):'6',(1,5):'10',(1,7):'14'}
for name,agg,j in (('C_d (d | mn)',aggC,jC),("C'_d (d | m^2-n^2, d>0)",aggCp,jCp)):
    print('\n==',name,' модуль',j)
    for cl in sorted(CL):
        yes=sorted((a,b) for (c,a,b),v in agg.items() if c==cl and True in v)
        no=sorted((a,b) for (c,a,b),v in agg.items() if c==cl and False in v)
        print(f' d ~ {CL[cl]:>2} (mod Q2^2): разрешимо при (m,n) mod {j} in {yes}; НЕразрешимо при {no}')
