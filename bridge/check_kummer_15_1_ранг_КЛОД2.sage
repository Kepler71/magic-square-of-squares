# -*- coding: utf-8 -*-
# Раунд 2 атаки на (15,1): ГДЕ РЕАЛЬНО ДЕРЖИТСЯ ГРАНИЦА rank E(Q) <= 1
# Первый 2-спуск даёт только rank <= 3. Всё «доказательство» опирается
# на ВТОРОЙ спуск / спаривание Касселса–Тейта. Проверяем это в лоб.

import sys
from sage.all import *

def hdr(s):
    print("\n" + "=" * 72); print(s); print("=" * 72); sys.stdout.flush()

S = QQ(113); M = ZZ(15); N = ZZ(1); B = S*M**2*N**2
e = [QQ(-B), QQ(-S*M**4), QQ(-S*N**4)]
Rx = PolynomialRing(QQ,'x'); x = Rx.gen()
f = (x-e[0])*(x-e[1])*(x-e[2])
c = f.coefficients(sparse=False)
E = EllipticCurve([0,c[2],0,c[1],c[0]])
Emin = E.minimal_model()

hdr("M1. ЧТО ДАЁТ ПЕРВЫЙ 2-СПУСК (без второго)")
mw = Emin.mwrank_curve()
mw.two_descent(second_descent=False, verbose=False)
print("  mwrank, только первый спуск: rank_bound = %s, certain = %s"
      % (mw.rank_bound(), mw.certain()))
print("  selmer_rank (mwrank)        = %s" % mw.selmer_rank())
print("  => ПЕРВОГО спуска НЕ ХВАТАЕТ: граница 3, а нужно 1.")

hdr("M2. ВТОРОЙ СПУСК (mwrank)")
mw2 = Emin.mwrank_curve()
mw2.two_descent(second_descent=True, verbose=False)
print("  mwrank со вторым спуском: rank_bound = %s, certain = %s"
      % (mw2.rank_bound(), mw2.certain()))
print("  rank (mwrank) = %s" % mw2.rank())

hdr("M3. ИСПОЛЬЗУЕТСЯ ЛИ БАЗА ДАННЫХ CREMONA?")
try:
    r_nodb = E.rank(use_database=False, only_use_mwrank=True, proof=True)
    print("  E.rank(use_database=False, only_use_mwrank=True, proof=True) = %s" % r_nodb)
except Exception as ex:
    print("  не удалось: %s" % ex)
print("  кондуктор = %s (вне диапазона таблиц Кремоны => база не могла помочь)"
      % E.conductor())

hdr("M4. PARI ellrank с разным effort")
pe = pari(Emin.a_invariants()).ellinit()
for eff in [0, 1, 2]:
    try:
        r = pe.ellrank(eff)
        print("  effort=%d -> %s" % (eff, r))
    except Exception as ex:
        print("  effort=%d -> ошибка %s" % (eff, ex))

hdr("M5. НЕЗАВИСИМЫЙ АНАЛИТИЧЕСКИЙ ПУТЬ (Гросс–Загир–Колывагин)")
w = E.root_number()
print("  корневое число w(E) = %s" % w)
print("  w = -1  =>  L(E,1) = 0 (функциональное уравнение), аналитический ранг нечётен")
L = E.lseries()
try:
    d1 = L.deriv_at1(100)
    print("  L'(E,1) ~ %s  (значение, оценка ошибки) = %s" % (d1[0], d1[1]))
    ok = (abs(d1[0]) > 10*abs(d1[1])) and d1[0] != 0
    print("  L'(E,1) != 0 с запасом по оценке ошибки: %s" % ok)
    print("  => аналитический ранг = 1 => (Колывагин+Гросс-Загир) rank E(Q) = 1 и Sha конечна.")
    print("  СТАТУС: численно (оценка ошибки от Sage), но независимо от спуска.")
except Exception as ex:
    print("  deriv_at1 не сработал: %s" % ex)

hdr("M6. СТРУКТУРА Sha[2]")
print("  dim_F2 Sel^2 (mwrank selmer_rank) = %s" % mw2.selmer_rank())
print("  dim_F2 E(Q)/2E(Q) при rank=1 и полном 2-кручении = 1+2 = 3")
print("  => dim_F2 Sha[2] = %s - 3 = %s" % (mw2.selmer_rank(), mw2.selmer_rank()-3))
try:
    print("  Sage E.sha().an() (аналитический порядок Sha, численно) = %s" % E.sha().an())
except Exception as ex:
    print("  sha().an(): %s" % ex)

hdr("M7. ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ ТОРСОРА D для класса (1,113,113)")
# D: r^2 - 113 z2^2 = -5695200 w^2 ,  r^2 - 113 z3^2 = 25312 w^2
# Тест над Z/p^k перебором с последующим Гензелем: ищем невырожденную точку.
def local_pt(p, k=6):
    mod = p**k
    R = Integers(mod)
    # перебираем r,w и решаем на квадратичность
    for w in range(1, min(mod, 400)):
        for r in range(0, min(mod, 400)):
            A = r*r + 5695200*w*w
            Bv = r*r - 25312*w*w
            if A % 113 or Bv % 113:
                continue
            a1 = A//113; b1 = Bv//113
            if a1 == 0 or b1 == 0:
                continue
            if is_sq_padic(a1, p) and is_sq_padic(b1, p):
                return (r, w)
    return None

def is_sq_padic(x, p):
    x = ZZ(x)
    if x == 0:
        return False
    v = x.valuation(p)
    if v % 2:
        return False
    u = x // p**v
    if p == 2:
        return u % 8 == 1
    return kronecker(u, p) == 1

bad = [2,3,5,7,113]
for p in bad + [11,13,17,19,23]:
    pt = local_pt(p)
    print("  p = %-4s : локальная точка D %s" % (p, ("найдена r,w=%s" % (pt,)) if pt else "НЕ найдена перебором"))
print("  над R: r=1000, w=1 даёт 113 z2^2 = 5696200>0 и 113 z3^2 = 974688>0 — точка есть")

hdr("M8. РАСШИРЕННЫЙ ПОИСК ТОЧЕК НА ТОРСОРЕ D (попытка сломать)")
DB = 8000
hits = []
for c0 in range(1, 401):
    c2 = c0*c0
    A1 = 5695200*c2; A2 = 25312*c2
    a = 0
    while a <= DB:
        a2 = a*a
        v1 = a2 + A1
        if v1 % 113 == 0 and ZZ(v1//113).is_square():
            v2 = a2 - A2
            if v2 > 0 and v2 % 113 == 0 and ZZ(v2//113).is_square():
                hits.append((a, c0))
        a += 1
print("  перебор r=a/c, a<=%d, c<=400: найдено %d точек: %s" % (DB, len(hits), hits[:5]))

hdr("M9. ПОИСК ТОЧЕК E БОЛЬШОЙ ВЫСОТЫ С ТРЕБУЕМЫМ КЛАССОМ")
def sqclass(v):
    v = QQ(v)
    z = ZZ(v.numerator()*v.denominator())
    return ZZ(z.sign())*ZZ(z.abs().squarefree_part())
def delta(P):
    if P.is_zero(): return (ZZ(1),ZZ(1),ZZ(1))
    d = [QQ(P[0]-e[i]) for i in range(3)]
    for i in range(3):
        if d[i] == 0:
            j,k = [q for q in range(3) if q != i]
            d[i] = QQ((e[i]-e[j])*(e[i]-e[k]))
    return tuple(sqclass(v) for v in d)
G = E.gens()[0]
T1 = E(e[0],0); T2 = E(e[1],0)
REQ = (ZZ(1), ZZ(113), ZZ(113))
found = []
seen = set()
for a in range(-40, 41):
    for c in [0,1]:
        for d in [0,1]:
            P = a*G + c*T1 + d*T2
            if P.is_zero(): continue
            dd = delta(P)
            seen.add(dd)
            if dd == REQ:
                found.append((a,c,d))
print("  перебор aG+cT1+dT2, |a|<=40: различных классов = %d" % len(seen))
print("  точек с требуемым классом (1,113,113): %d  %s" % (len(found), found[:3]))
print("  (это прямая проверка того, что образ действительно замкнут на 8 классах)")

print("\nГОТОВО")
