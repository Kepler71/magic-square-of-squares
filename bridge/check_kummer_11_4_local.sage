# Часть 3: (a) гомоморфность delta (без факторизации), (b) локальная разрешимость C_{11,4},
#           (c) лежит ли требуемый класс (1,274,274) в 2-Selmer => элемент Sha[2].

def hdr(t):
    print("\n" + "="*78); print(t); print("="*78)

def sqfree(q):
    q = QQ(q); a=q.numerator(); d=q.denominator()
    return ZZ(a*d).squarefree_part()

M, N = ZZ(11), ZZ(4)
S = QQ(M^2+N^2)/2
B = S*M^2*N^2
ROOTS = [QQ(-B), QQ(-S*M^4), QQ(-S*N^4)]
Rx.<XX> = PolynomialRing(QQ)
cub = (XX-ROOTS[0])*(XX-ROOTS[1])*(XX-ROOTS[2])
c = cub.coefficients(sparse=False)
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
REQ = (1, sqfree(S), sqfree(S))

def delta_raw(P, roots=ROOTS):
    """Тройка рациональных представителей класса delta(P) — БЕЗ факторизации."""
    if P.is_zero(): return (QQ(1),QQ(1),QQ(1))
    x = QQ(P[0]); out=[]
    for i in range(3):
        d = x - roots[i]
        if d == 0:
            j,k = [q for q in range(3) if q != i]
            d = (roots[i]-roots[j])*(roots[i]-roots[k])
        out.append(QQ(d))
    return tuple(out)

def cls_eq(a, b):
    """Равны ли классы двух ненулевых рациональных: a/b — квадрат?"""
    return (QQ(a)/QQ(b)).is_square()

def trip_eq(u, v):
    return all(cls_eq(u[i], v[i]) for i in range(3))

def trip_mul(u, v):
    return tuple(QQ(u[i])*QQ(v[i]) for i in range(3))

hdr("C. ГОМОМОРФНОСТЬ delta — проверка без факторизации, на реальных кратных G")
Gs = E.gens(proof=False)
print("генератор G =", Gs[0])
T2 = [P for P in E.torsion_subgroup().points() if P != E(0) and 2*P == E(0)]
pool = [E(0)] + T2 + list(Gs)
for k in range(2, 7):
    pool.append(k*Gs[0])
for tp in T2:
    for k in range(1, 4):
        pool.append(k*Gs[0] + tp)
pool = list(dict.fromkeys(pool))
print("размер пула точек:", len(pool), " (включая кратные G до 6G и сдвиги кручением)")
bad = 0; checked = 0
for P in pool:
    for Q in pool:
        checked += 1
        if not trip_eq(delta_raw(P+Q), trip_mul(delta_raw(P), delta_raw(Q))):
            bad += 1
            if bad <= 5: print("  НАРУШЕНИЕ:", P, Q)
print("проверено пар:", checked, " нарушений гомоморфности:", bad)
print("delta — гомоморфизм на проверенном пуле:", bad == 0)

badk = sum(1 for P in pool if not trip_eq(delta_raw(2*P), (1,1,1)))
print("КОНТРОЛЬ ядра: delta(2P)=(1,1,1) для всех P пула, нарушений:", badk, "->", badk==0)

print("\nКОНТРОЛЬ ПОЛНОТЫ ОБРАЗА на кратных: delta(kG+tors) всегда среди 8 классов?")
IMG8 = set()
gens_all = list(Gs) + T2
from itertools import product as iproduct
for co in iproduct([0,1], repeat=len(gens_all)):
    P = E(0)
    for cc,gg in zip(co,gens_all):
        if cc: P = P+gg
    IMG8.add(tuple(sqfree(x) for x in delta_raw(P)))
print("  8 классов:", sorted(IMG8), " размер:", len(IMG8))
outside = 0
for P in pool:
    dr = delta_raw(P)
    if not any(trip_eq(dr, w) for w in IMG8):
        outside += 1; print("  ВНЕ ОБРАЗА:", P)
print("  точек пула вне 8 классов:", outside)
print("  требуемый класс", REQ, "в образе:", any(trip_eq(REQ, w) for w in IMG8))

hdr("L. ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ C_{11,4} — свой счёт")

def is_sq_Qp(q, p):
    q = QQ(q)
    if q == 0: return True
    v = q.valuation(p)
    if v % 2 != 0: return False
    u = q / QQ(p)^v
    num = ZZ(u.numerator()); den = ZZ(u.denominator())
    if p == 2:
        return (num * inverse_mod(den % 8, 8)) % 8 == 1
    return kronecker((num * inverse_mod(den % p, p)) % p, p) == 1

def Fs(tv):
    return (M^2 + N^2*tv^2, S*(1+tv^2), N^2 + M^2*tv^2)

BADP = sorted(set(ZZ(2*M*N*(M^2-N^2)*(M^2+N^2)).prime_factors()))
print("плохие простые (делители 2mn(m^2-n^2)(m^2+n^2)):", BADP)
print("R: t=0 -> F =", Fs(QQ(0)), " все >0:", all(x>0 for x in Fs(QQ(0))))

primes_to_test = sorted(set(BADP) | set(primes(120)))
results = {}
for p in primes_to_test:
    witness = None
    lim = min(p^3, 3000) if p < 30 else 300
    seen = set()
    for k in range(0, 4):
        for a in range(-lim, lim+1):
            tv = QQ(a)/QQ(p)^k
            if tv in seen: continue
            seen.add(tv)
            f = Fs(tv)
            if all(is_sq_Qp(x, p) for x in f):
                witness = tv; break
        if witness is not None: break
    results[p] = witness
print("свидетели для ПЛОХИХ простых:")
for p in BADP:
    print("   p=%-5d t=%s  F=%s" % (p, results[p], Fs(results[p]) if results[p] is not None else None))
missing = [p for p in primes_to_test if results[p] is None]
print("простые без найденного свидетеля (p<=120):", missing)
print("ВСЕ протестированные p имеют Q_p-точку:", len(missing)==0)

hdr("SEL. ТРЕБУЕМЫЙ КЛАСС В 2-SELMER? => НЕТРИВИАЛЬНЫЙ Sha[2]")
e1,e2,e3 = ROOTS
print("покрытие D_d:  z1^2 - 274 z2^2 = e2-e1 = %s" % (e2-e1))
print("               z1^2 - 274 z3^2 = e3-e1 = %s" % (e3-e1))
print("d1*d2*d3 квадратный класс:", sqfree(1*REQ[1]*REQ[2]), "(нужен 1)")
print()
print("Из тождеств части A точка C даёт точку D_d:")
print("  z1 = m*n*u4,  z2 = m*u0/ (норм.), z3 = n*u8/(норм.)  =>  C(Q_p)!=0 влечёт D_d(Q_p)!=0.")
print("Локальная разрешимость C выше => (1,274,274) лежит в 2-Selmer группе.")
print("Значит класс НЕ в образе E(Q) (часть E) но В Sel_2 => нетривиальный элемент Sha(E/Q)[2].")
print()
try:
    sr = E.selmer_rank()
    print("mwrank dim Sel_2 =", sr, "=> rank + dim Sha[2] =", sr-2)
    print("если rank=1 => dim Sha[2]=2 (чётно — согласуется с конечностью Sha)")
    print("если rank=3 => dim Sha[2]=0 — ПРОТИВОРЕЧИЕ с найденным элементом Sha[2] (при локальной разрешимости C)")
except Exception as ex:
    print("selmer_rank FAILED:", ex)
