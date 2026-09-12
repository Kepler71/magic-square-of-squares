# -*- coding: utf-8 -*-
# РЕШАЮЩИЙ РАСЧЁТ для (m,n)=(11,4).
# 2-спуск (eclib) даёт только 1 <= rank(E) <= 3, rk Sel^2 = 5, #Sha[2] <= 4.
# Поэтому «8 из 8 классов» НЕ доказывает полноту образа, пока ранг не равен 1 ДОКАЗАННО.
#
# Здесь считается САМ полный 2-Селмер (локальные образы во всех плохих местах)
# и проверяется, лежит ли требуемый класс (1,274,274) в Sel^2(E/Q).
#   - если НЕТ  -> исключение (11,4) верно БЕЗУСЛОВНО (ранг не нужен вовсе);
#   - если ДА   -> вывод целиком висит на ранге; при rank=3 образ = Селмер
#                  и требуемый класс ДОСТИГАЕТСЯ, т.е. заявление Codex ложно.
import sys, itertools
from sage.all import *
def flush(): sys.stdout.flush()

m=11; n=4
s=QQ(m^2+n^2)/2; b=s*m^2*n^2
f=[ZZ(-4*b), ZZ(-4*s*m^4), ZZ(-4*s*n^4)]     # целая модель, x=4X (4 — квадрат)
print("f =", f)
xv=polygen(QQ,'x'); cub=(xv-f[0])*(xv-f[1])*(xv-f[2])
c2=ZZ(cub[2]); c1=ZZ(cub[1]); c0=ZZ(cub[0])
E=EllipticCurve(QQ,[0,c2,0,c1,c0])
BAD = sorted(set(E.conductor().prime_factors()) | {2})
print("плохие простые:", BAD, " N =", factor(E.conductor()))

def sqfree(q):
    q=QQ(q); v=ZZ(q.numerator()*q.denominator()); sg=1 if v>0 else -1
    r=ZZ(1)
    for p,e in factor(abs(v)):
        if e%2==1: r*=p
    return sg*r

TARGET = (ZZ(1), sqfree(s), sqfree(s))
print("ТРЕБУЕМЫЙ КЛАСС:", TARGET)

# --------- локальные квадратные классы ---------
def cls_p(a, p):
    """класс a in Q_p^*/(Q_p^*)^2: (v mod 2, единичная часть)"""
    a=QQ(a); assert a!=0
    v=a.valuation(p)
    u=a/QQ(p)**v
    num=ZZ(u.numerator()); den=ZZ(u.denominator())
    if p==2:
        uu=(num*inverse_mod(den%8,8))%8
        return (v%2, uu)
    uu=(num*inverse_mod(den%p,p))%p
    return (v%2, 1 if kronecker(uu,p)==1 else -1)
def cls_R(a):
    return (1 if QQ(a)>0 else -1,)
def mul_cls(x,y,p):
    if p=='R': return ((x[0]*y[0]),)
    if p==2:  return ((x[0]+y[0])%2, (x[1]*y[1])%8)
    return ((x[0]+y[0])%2, x[1]*y[1])
def trip(a3, p):
    if p=='R': return tuple(cls_R(a) for a in a3)
    return tuple(cls_p(a,p) for a in a3)
def trip_mul(u,v,p):
    return tuple(mul_cls(u[i],v[i],p) for i in range(3))
def trip_one(p):
    if p=='R': return ((1,),(1,),(1,))
    if p==2:  return ((0,1),(0,1),(0,1))
    return ((0,1),(0,1),(0,1))

def delta_of_x(x):
    """тройка (x-f1, x-f2, x-f3) как рациональные числа; x не корень"""
    return [QQ(x)-f[0], QQ(x)-f[1], QQ(x)-f[2]]
def delta_tors(i):
    j,k=[t for t in range(3) if t!=i]
    out=[QQ(0)]*3
    for t in range(3):
        out[t] = (f[i]-f[j])*(f[i]-f[k]) if t==i else QQ(f[i]-f[t])
    return out

def is_sq_p(a,p):
    a=QQ(a)
    if a==0: return True
    if p=='R': return a>0
    c=cls_p(a,p)
    if p==2: return c==(0,1)
    return c==(0,1)

# --------- локальные образы delta_p ---------
expected = {'R':2, 2:8}
for p in BAD:
    if p!=2: expected[p]=4

def gen_group(elts, p):
    G={trip_one(p)}
    changed=True
    while changed:
        changed=False
        for a in list(G):
            for g in elts:
                h=trip_mul(a,g,p)
                if h not in G: G.add(h); changed=True
    return G

def local_image(p):
    gens=[]
    # кручение
    for i in range(3):
        gens.append(trip(delta_tors(i), p))
    cand=[]
    if p=='R':
        cand=[QQ(x) for x in [f[2]+1, f[2]+10**6, 0, 10**8, f[1]+1, f[1]+10, (f[1]+f[0])/2]]
    else:
        cand=[]
        for a in range(-8, 9):
            for u in range(1, min(p**3, 2000)+1):
                if u % p == 0: continue
                for sg in [1,-1]:
                    cand.append(sg*QQ(p)**a*u)
        for i in range(3):
            for j in range(-6,7):
                for u in [1,2,3,5,7,-1,-2,-3]:
                    cand.append(QQ(f[i])+u*QQ(p)**j)
        cand += [QQ(x) for x in range(-5000,5001)]
    got=set(gen_group(gens,p))
    for x in cand:
        d=delta_of_x(x)
        if any(t==0 for t in d): continue
        prod=d[0]*d[1]*d[2]
        if not is_sq_p(prod,p): continue           # нужна Q_p-точка над x
        t=trip(d,p)
        if t not in got:
            got=gen_group(list(got)+[t],p)
            if len(got)>=expected[p]: break
    return got

print()
print("=== локальные образы delta_p: E(Q_p)/2E(Q_p) ===")
LOC={}
for p in ['R']+BAD:
    im=local_image(p)
    LOC[p]=im
    print("  место %-4s: |образ| = %2d  (теоретически %d)  %s"
          % (str(p), len(im), expected[p],
             "OK" if len(im)==expected[p] else "<<< НЕ НАБРАН, расчёт неполон"))
    flush()

print()
print("=== лежит ли требуемый класс в локальных образах? ===")
in_all=True
for p in ['R']+BAD:
    t=trip([QQ(TARGET[0]),QQ(TARGET[1]),QQ(TARGET[2])], p)
    ok = t in LOC[p]
    print("  место %-4s: (1,274,274) в образе ? %s" % (str(p), ok))
    if not ok: in_all=False
print()
print(">>> ТРЕБУЕМЫЙ КЛАСС В 2-СЕЛМЕРЕ ?", in_all)
if not in_all:
    print(">>> ЗНАЧИТ исключение (11,4) ВЕРНО БЕЗУСЛОВНО, ранг не нужен.")
else:
    print(">>> ЗНАЧИТ вывод зависит ТОЛЬКО от того, что rank(E)=1.")
    print(">>> При rank(E)=3 образ = Селмер и класс ДОСТИГАЕТСЯ => заявление Codex было бы ложным.")
flush()

# --------- полный Селмер: сверка с mwrank (rk S^2 = 5) ---------
print()
print("=== полный перебор 2-Селмера (сверка размера с mwrank) ===")
supp=[-1]+BAD
def from_bits(bits):
    r=ZZ(1)
    for i,v in enumerate(bits):
        if v: r*=supp[i]
    return r
allsq=[from_bits(bits) for bits in itertools.product([0,1],repeat=len(supp))]
sel=[]
for d1 in allsq:
    for d2 in allsq:
        d3=sqfree(d1*d2)
        ok=True
        for p in ['R']+BAD:
            if trip([QQ(d1),QQ(d2),QQ(d3)],p) not in LOC[p]:
                ok=False; break
        if ok: sel.append((d1,d2,d3))
print("  |Sel^2(E/Q)| (мой расчёт) =", len(sel), " -> rk =", log(len(sel),2))
print("  mwrank: rk(S^2(E)) = 5 -> |Sel| = 32.  СОВПАДАЕТ ?", len(sel)==32)
print("  требуемый класс в моём Селмере ?", (TARGET[0],TARGET[1],TARGET[2]) in
      set((ZZ(a),ZZ(b),ZZ(c)) for a,b,c in sel))
print()
print("  элементы Sel^2 (d1,d2,d3):")
for e in sorted(sel):
    mark = "  <<< ТРЕБУЕМЫЙ" if (ZZ(e[0]),ZZ(e[1]),ZZ(e[2]))==TARGET else ""
    print("   ", e, mark)
