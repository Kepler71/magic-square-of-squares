#!/usr/bin/env python3
# Fable, 14.09.2026. Локальная разрешимость с девятью различными ненулевыми квадратами: R, Q_p (все p), сравнения по модулю.
from sage.all import *
import json, itertools
out={}
# 1. Необходимые сравнения для целого примитивного магического квадрата из квадратов: перебор по модулю 8, 3, 5, 7
def magic_mod(N):
    sq=sorted(set((x*x)%N for x in range(N)))
    res=set()
    for c in sq:
        for b in range(N):
            for cc in range(N):
                cells=[c-b,c+b,c-cc,c+cc,c-b-cc,c+b+cc,c-b+cc,c+b-cc]
                if all(v%N in sq for v in cells):
                    res.add((c,b%N,cc%N))
    return res
for N in [8,3,5,7,16,9]:
    R=magic_mod(N)
    prim=[(c,b,cc) for (c,b,cc) in R if gcd(gcd(c,b),gcd(cc,N))==1]
    ub=set();
    for (c,b,cc) in prim:
        ub.add((b%N,cc%N))
    print(f"mod {N}: решений {len(R)}, примитивных (gcd(c,b,c')=1 mod N) {len(prim)}, возможные (b,c) mod N: {sorted(ub)[:12]}{'...' if len(ub)>12 else ''}")
# ожидание: mod 8: b,c ≡ 0 (mod 8) при нечётном c; mod 3: b,c ≡ 0; mod 5: либо b≡c≡0, либо 5|c
# 2. Явные Q_p-точки для всех p<=200: u=(b,c,b+c,b-c) с b=p^k, c=3p^k (p нечётное), b=8,c=24 (p=2)
def is_sq_Qp(v,p,prec=60):
    return Qp(p,prec)(v).is_square()
cert={}
for p in prime_range(2,201):
    if p==2: b,c=8,24
    else: b,c=p,3*p
    us=[b,c,b+c,b-c]
    cells=[1]+[1-v for v in us]+[1+v for v in us]
    ok=all(is_sq_Qp(v,p) for v in cells) and len(set(cells))==9 and 0 not in cells
    cert[p]=(b,c,ok)
    assert ok, (p,cells)
print("Q_p-точки с девятью различными ненулевыми квадратами: все p<=200 — OK; сертификат (b,c):", {p:v[:2] for p,v in list(cert.items())[:6]}, "...")
out['Qp']={int(p):[int(v[0]),int(v[1])] for p,v in cert.items()}
# доказательство для всех p: 1+pZ_p ⊂ Q_p^2 при p нечётном, 1+8Z_2 ⊂ Q_2^2; b=p, c=3p дают |u| = p,3p,4p,2p попарно различные и ненулевые.
# 3. Решения по модулю p^k с девятью различными квадратами целых чисел (явные корни): p^k ~ 10^6
modcert={}
for p in prime_range(2,60):
    k=1
    while p**k<10**6: k+=1
    N=p**k
    b,c=cert[p][0],cert[p][1]
    us=[b,c,b+c,b-c]
    cells=[1]+[1-v for v in us]+[1+v for v in us]
    roots=[]
    for v in cells:
        r=Zp(p,k+2)(v).sqrt()
        r=ZZ(r)%N if p!=2 else ZZ(r)%N
        assert (r*r-v)%N==0
        roots.append(int(r))
    assert len(set(v%N for v in cells))==9
    modcert[int(p)]=dict(N=int(N),cells=[int(v%N) for v in cells],roots=roots)
print("по модулю p^k>10^6, p<60: девять различных ненулевых квадратов целых чисел, магические суммы — OK, напр. p=7:", modcert[7])
out['mod_pk']=modcert
# 4. Совместно по модулю M=2^7*3^4*5^3*7^3*11^2*13^2 (КТО): один квадрат из целых чисел, магический по модулю M с девятью различными квадратами
ps=[(2,7),(3,4),(5,3),(7,3),(11,2),(13,2)]
M=prod(p**k for p,k in ps)
Bs=crt([cert[p][0] for p,k in ps],[p**k for p,k in ps]); Cs=crt([cert[p][1] for p,k in ps],[p**k for p,k in ps])
us=[Bs,Cs,Bs+Cs,Bs-Cs]; cells=[1]+[1-v for v in us]+[1+v for v in us]
roots=[]
for v in cells:
    rs=[ZZ(Zp(p,k+2)(v).sqrt())%p**k for p,k in ps]
    r=crt(rs,[p**k for p,k in ps]); assert (r*r-v)%M==0; roots.append(int(r%M))
print(f"по модулю M={M}: клетки {[int(v%M) for v in cells]}, корни {roots}, различны: {len(set(v%M for v in cells))==9}")
out['crt']=dict(M=int(M),b=int(Bs),c=int(Cs),roots=roots)
# 5. Вещественная точка: b=3/10, c=1/10 (см. s3) — очевидно.
# 6. Целочисленное следствие: у примитивного магического квадрата из квадратов все u ≡ 0 (mod 24) — проверка перебором выше (mod 8, mod 3)
json.dump(out, open('/home/kep/magicKube/fable_symmetry/s4_local.json','w'), indent=1)
print("OK s4")
