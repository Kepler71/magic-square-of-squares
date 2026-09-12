# -*- coding: utf-8 -*-
# Часть 3: перечисление классов delta СВОЕЙ точкой, своей арифметикой. Без вызова E.rank().
from sage.all import *
import sys

def hdr(t):
    print("\n" + "="*78); print(t); print("="*78); sys.stdout.flush()

def sqclass(x):
    x = QQ(x)
    if x == 0: raise ValueError("zero")
    z = ZZ(x.numerator())*ZZ(x.denominator())
    out = 1
    for (p,e) in factor(abs(z)):
        if e % 2: out *= p
    return out if z > 0 else -out

def tmul(A,B): return tuple(sqclass(QQ(a)*QQ(c)) for a,c in zip(A,B))

m, n = 19, 5
s = QQ(m**2+n**2)/2; b = s*m**2*n**2
e1, e2, e3 = -b, -s*m**4, -s*n**4
A2 = -(e1+e2+e3); A4 = e1*e2+e1*e3+e2*e3; A6 = -e1*e2*e3
EE = [e1,e2,e3]
req = (sqclass((m*n)**2), sqclass(s*m**2), sqclass(s*n**2))
print("s =", s, " b =", b, " e =", (e1,e2,e3), " требуемый класс =", req)

# ---- своя арифметика
def on_curve(P):
    if P is None: return True
    x,y = P; return y*y == x**3 + A2*x*x + A4*x + A6
def neg(P): return None if P is None else (P[0], -P[1])
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1 = P; x2,y2 = Q
    if x1 == x2:
        if y1 + y2 == 0: return None
        lam = (3*x1*x1 + 2*A2*x1 + A4)/(2*y1)
    else:
        lam = (y2-y1)/(x2-x1)
    x3 = lam*lam - A2 - x1 - x2
    return (x3, lam*(x1-x3) - y1)
def mul(k,P):
    k = ZZ(k)
    if k < 0: return mul(-k, neg(P))
    R_, Q_ = None, P
    while k > 0:
        if k & 1: R_ = add(R_, Q_)
        Q_ = add(Q_, Q_); k >>= 1
    return R_
def delta(P):
    if P is None: return (1,1,1)
    x,y = P; out = []
    for i in range(3):
        if x == EE[i]:
            j,k = [q for q in range(3) if q != i]
            out.append(sqclass((EE[i]-EE[j])*(EE[i]-EE[k])))
        else:
            out.append(sqclass(x - EE[i]))
    return tuple(out)

T1 = (e1, QQ(0)); T2 = (e2, QQ(0)); T3 = (e3, QQ(0))
assert on_curve(T1) and on_curve(T2) and on_curve(T3)
assert add(T1,T2) == T3

# ---- МОЯ точка: беру ту, что выдал МОЙ прогон PARI при effort=2 (она ОТЛИЧАЕТСЯ от точки Codex)
MY = (QQ(3557429375)/QQ(14641), QQ(7584055446240000)/QQ(1771561))
CODEX = (QQ(-28561963657)/QQ(1369), QQ(-2089086742828800)/QQ(50653))
ALT = (QQ(11775817140519)/QQ(105625), QQ(45108971173988070272)/QQ(34328125))

hdr("[1] ПРОВЕРКА ТОЧЕК ПОДСТАНОВКОЙ")
for nm, P in [("MY (effort2)", MY), ("CODEX (effort1)", CODEX), ("ALT (effort4)", ALT)]:
    print("  %-16s на кривой: %s" % (nm, on_curve(P)))

hdr("[2] СООТНОШЕНИЯ МЕЖДУ ТРЕМЯ ТОЧКАМИ (своей арифметикой)")
Tors = [None, T1, T2, T3]
def find_rel(P, Q, K=12):
    """ищем k,T: Q = k*P + T"""
    for k in range(-K, K+1):
        if k == 0: continue
        kP = mul(k, P)
        for T in Tors:
            if add(kP, T) == Q:
                return (k, T)
    return None
print("  CODEX = k*MY + T ?  ->", find_rel(MY, CODEX))
print("  ALT   = k*MY + T ?  ->", find_rel(MY, ALT))
print("  MY    = k*CODEX+T ? ->", find_rel(CODEX, MY))
print("  ALT   = k*CODEX+T ? ->", find_rel(CODEX, ALT))

hdr("[3] МОЯ ТАБЛИЦА ВОСЬМИ КЛАССОВ (точка MY)")
rows = []
for a in [0,1]:
    for c in [0,1]:
        for d in [0,1]:
            P = add(mul(a,MY), add(mul(c,T1), mul(d,T2)))
            assert on_curve(P)
            rows.append(((a,c,d), P, delta(P)))
print("  | (a,c,d) |     d1 |     d2 |     d3 |")
for (k,P,D) in rows:
    print("  | %s | %6s | %6s | %6s |" % (k, D[0], D[1], D[2]))
cls = [D for (_,_,D) in rows]
print("\n  различных:", len(set(cls)), "из", len(cls))
print("  требуемый", req, "присутствует:", req in set(cls))

hdr("[4] ТАБЛИЦА ПО ТОЧКЕ CODEX (мой пересчёт)")
rows_c = []
for a in [0,1]:
    for c in [0,1]:
        for d in [0,1]:
            P = add(mul(a,CODEX), add(mul(c,T1), mul(d,T2)))
            rows_c.append(((a,c,d), delta(P)))
for (k,D) in rows_c: print("  ", k, D)
cls_c = [D for (_,D) in rows_c]
print("  различных:", len(set(cls_c)))
print("  требуемый присутствует:", req in set(cls_c))
codex_table = [(1,1,1),(-1,4053,-4053),(-4053,386,-42),(4053,42,386),
               (-3,386,-1158),(3,42,14),(1351,1,1351),(-1351,4053,-3)]
print("  множество == таблице Codex:", set(cls_c) == set(codex_table))
print("  построчно  == таблице Codex:", cls_c == codex_table)
print("  множество МОИХ == множеству Codex:", set(cls) == set(cls_c))

hdr("[5] ТАБЛИЦА ПО ТОЧКЕ ALT (третья, независимая выборка)")
rows_a = []
for a in [0,1]:
    for c in [0,1]:
        for d in [0,1]:
            P = add(mul(a,ALT), add(mul(c,T1), mul(d,T2)))
            rows_a.append(((a,c,d), delta(P)))
for (k,D) in rows_a: print("  ", k, D)
cls_a = [D for (_,D) in rows_a]
print("  различных:", len(set(cls_a)), " требуемый присутствует:", req in set(cls_a))
print("  множество == моему:", set(cls_a) == set(cls))

hdr("[6] ГОМОМОРФНОСТЬ delta НА 64 ПАРАХ (моя точка)")
P8 = {k:P for (k,P,_) in rows}; D8 = {k:D for (k,_,D) in rows}
bad = 0
for k1 in P8:
    for k2 in P8:
        if delta(add(P8[k1],P8[k2])) != tmul(D8[k1], D8[k2]):
            bad += 1; print("   НАРУШЕНИЕ", k1, k2)
print("  нарушений:", bad, "из 64")
print("  delta(2P) для восьми P:", [delta(mul(2,P)) for (_,P,_) in rows])

hdr("[7] ГРУППОВАЯ СТРУКТУРА ОБРАЗА")
G8 = set(cls)
closed = all(tmul(x,y) in G8 for x in G8 for y in G8)
print("  множество восьми классов замкнуто относительно умножения:", closed)
print("  элементарная 2-группа (x*x=1):", all(tmul(x,x)==(1,1,1) for x in G8))
print("  требуемый класс req умноженный на каждый элемент образа:")
for x in sorted(G8):
    print("     ", x, "*", req, "=", tmul(x, req))
print("  req лежит в образе:", req in G8)

hdr("[8] ЛЕЖИТ ЛИ req В 2-СЕЛМЕРЕ? (тогда ранговая граница НЕСУЩАЯ)")
# торсор: z1^2 - 193 z2^2 = e2-e1 ; z1^2 - 193 z3^2 = e3-e1
c12 = e2 - e1; c13 = e3 - e1
print("  e2-e1 =", c12, "=", factor(c12))
print("  e3-e1 =", c13, "=", factor(c13))
print("  явная связь с C:  z1=95*u4, z2=19*u0, z3=5*u8, x=b t^2")
# проверим на локальных точках C (ищем t в Z_p с квадратными F_i)
def loc_C(p, prec=40):
    K = Qp(p, prec=prec)
    for tt in range(0, p**3 if p < 20 else 4000):
        F0 = K(m**2 + n**2*tt**2); F4 = K(s*(1+tt**2)); F8 = K(n**2 + m**2*tt**2)
        if F0 == 0 or F4 == 0 or F8 == 0: continue
        if F0.is_square() and F4.is_square() and F8.is_square():
            return tt
    return None
bad_p = [2,3,5,7,19,193]
print("  локальные свидетели C(Q_p) (целое t):")
for p in bad_p + [11,13,17,23,29,31,37,41,43]:
    w = loc_C(p)
    print("    p=%-4s t=%s" % (p, w))
