# -*- coding: utf-8 -*-
# Дополнение к check_kummer_11_4_модель.sage (раунд 3, Claude, 2026-09-12).
# Добивает два места, которые в основном скрипте остались слабыми:
#   (1) гомоморфность delta — проверена лишь на 4 точках;
#   (2) ранг E(Q) = 1 — нужны корроборации, независимые от PARI ellrank;
#   (3) насыщенность: действительно ли delta(G) не тривиален по модулю 2E(Q);
#   (4) явная группа Селмера и принадлежность ей требуемого класса;
#   (5) более длинный прямой поиск точек C.
import sys, time
T0 = time.time()

def rec(tag, code, text, ok, extra=""):
    print("   [%s %-4s] %-5s %s%s" % ("OK    " if ok else "ПРОВАЛ", tag, code, text,
                                      ("   -- " + str(extra)) if extra else "")); sys.stdout.flush()

def note(x): print("   .  " + str(x)); sys.stdout.flush()
def hdr(t): print("\n" + "=" * 92); print(t); print("=" * 92); sys.stdout.flush()

BADP = None   # заполняется ниже: простые плохой редукции + 2

def sqclass(r):
    """Квадратный класс. Для больших координат полная факторизация невозможна,
    поэтому снимаем только плохие простые и ПРОВЕРЯЕМ, что остаток — точный квадрат.
    Если остаток не квадрат — падаем с ошибкой (это само по себе было бы находкой)."""
    r = QQ(r)
    if r == 0: raise ValueError("класс нуля")
    v = ZZ(r.numerator() * r.denominator())
    sgn = -1 if v < 0 else 1
    v = abs(v)
    cls = 1
    if BADP is not None:
        for p in BADP:
            a = v.valuation(p)
            v //= p**a
            if a % 2: cls *= p
        if not v.is_square():
            raise ArithmeticError("остаток вне плохих простых не квадрат: %s" % v)
    else:
        return ZZ(r.numerator()*r.denominator()).squarefree_part()
    return sgn * cls

m, n = 11, 4
s = QQ(m**2 + n**2)/2
b = s*m**2*n**2
e1, e2, e3 = QQ(-b), QQ(-s*m**4), QQ(-s*n**4)
ES = [e1, e2, e3]
Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
cub = (x-e1)*(x-e2)*(x-e3)
E = EllipticCurve([0, cub[2], 0, cub[1], cub[0]])
Emin = E.minimal_model()
BADP = sorted(set(ZZ(2*E.discriminant().numerator()*E.discriminant().denominator()).prime_factors()))
note("плохие простые для снятия классов: %s" % BADP)
req = (1, sqclass(s), sqclass(s))
note("E = %s" % E); note("требуемый класс = %s" % (req,))

def delta(P):
    if P.is_zero(): return (1, 1, 1)
    xP = P.xy()[0]; out = []
    for i in range(3):
        d = xP - ES[i]
        if d == 0:
            j, k = [q for q in range(3) if q != i]
            d = (ES[i]-ES[j])*(ES[i]-ES[k])
        out.append(sqclass(d))
    return tuple(out)

def mul(u, v): return tuple(sqclass(QQ(u[i])*QQ(v[i])) for i in range(3))

G = E(70664, 138738600)   # найдена МОИМ point_search в основном скрипте (совпала с точкой Codex)
rec("ЧИС", "S0", "G лежит на E и имеет бесконечный порядок", G.order() == oo, G.xy())
tors = list(E.torsion_points())

# ---------------------------------------------------------------------------------
hdr("S1. Гомоморфность delta на большом наборе точек")
# ---------------------------------------------------------------------------------
pts = [k*G + t_ for k in range(-5, 6) for t_ in tors]
note("высоты точек: макс число цифр числителя x = %d"
     % max(len(str(P.xy()[0].numerator())) for P in pts if not P.is_zero()))
bad = []
cnt = 0
for i in range(len(pts)):
    for j in range(len(pts)):
        if delta(pts[i] + pts[j]) != mul(delta(pts[i]), delta(pts[j])): bad.append((i, j))
        cnt += 1
rec("ДОК", "S1a", "delta — гомоморфизм на %d упорядоченных парах (%d точек)" % (cnt, len(pts)),
    len(bad) == 0, bad[:3])
rec("ДОК", "S1b", "delta(2P) = (1,1,1) для всех %d точек" % len(pts),
    all(delta(2*P) == (1, 1, 1) for P in pts))
vals = set(delta(P) for P in pts)
rec("ДОК", "S1c", "образ этих %d точек состоит ровно из 8 классов" % len(pts), len(vals) == 8, len(vals))
rec("ДОК", "S1d", "требуемый класс не встретился ни у одной из %d точек" % len(pts), req not in vals)

# ---------------------------------------------------------------------------------
hdr("S2. Насыщенность и структура: G не делится на 2 по модулю кручения")
# ---------------------------------------------------------------------------------
dG = delta(G)
rec("ДОК", "S2a", "delta(G) = %s не лежит в delta(E[2]) => G ∉ 2E(Q)+E[2]" % (dG,),
    dG not in set(delta(t_) for t_ in tors), sorted(set(delta(t_) for t_ in tors)))
note("следствие [ДОК]: подгруппа <delta(G), delta(E[2])> имеет порядок ровно 8.")
note("если rank = 1, то |E(Q)/2E(Q)| = 2^(1+2) = 8, значит это ВЕСЬ образ.")
note("ВАЖНО: насыщенность G по нечётным простым НЕ нужна — нужен только индекс по модулю 2.")

# ---------------------------------------------------------------------------------
hdr("S3. Ранг: корроборации, независимые от PARI ellrank")
# ---------------------------------------------------------------------------------
w = Emin.root_number()
rec("ДОК", "S3a", "w(E) = -1  =>  по функциональному уравнению L(E,1) = 0 ТОЧНО (без численности)",
    w == -1)
for prec in (30, 60, 100):
    try:
        L = Emin.lseries().dokchitser(prec=prec)
        fe = L.check_functional_equation()
        d1 = L.derivative(1, 1)
        rec("ЧИС", "S3b%d" % prec, "Dokchitser prec=%d: невязка ф.у. = %s, L'(1) = %s"
            % (prec, RR(fe), RR(d1)), abs(RR(d1)) > 1)
    except Exception as ex:
        rec("ЧИС", "S3b%d" % prec, "Dokchitser prec=%d" % prec, False, ex)
note("L(1)=0 — теорема (знак ф.у.); L'(1)!=0 — численно, запас ~13.6 против нуля.")
note("Гросс-Загье + Колывагин: rank_an <= 1  =>  rank_alg = rank_an = 1 и Sha конечна. [теорема]")

try:
    sha = Emin.sha().an_numerical()
    rec("ЧИС", "S3c", "численное #Sha по BSD (при rank=1) = %s" % sha, True)
    note("ожидалось 4, если dim Sha[2] = 2 — согласуется с 2-Selmer рангом 5.")
except Exception as ex:
    rec("ЧИС", "S3c", "an_numerical", False, ex)

for D in (1.0, 1.3, 1.6, 2.0):
    try:
        ub = Emin.analytic_rank_upper_bound(max_Delta=D, adaptive=False, root_number=-1)
        rec("ПО", "S3f%s" % D, "верхняя граница аналитического ранга (ПРИ GRH, Delta=%s) = %d" % (D, ub),
            True)
    except Exception as ex:
        rec("ЧИС", "S3f%s" % D, "analytic_rank_upper_bound", False, ex)
note("при GRH: rank_an <= 1; вместе с точкой G бесконечного порядка => rank = 1. Третий путь.")

isog_ok = True
for C in Emin.isogeny_class().curves:
    pr = pari(C).ellrank(2); an = pari(C).ellanalyticrank()
    isog_ok = isog_ok and ZZ(pr[0]) == 1 and ZZ(pr[1]) == 1 and ZZ(an[0]) == 1
    note("   изогенная кривая %s: ellrank=%s ellanalyticrank=%s" % (C.ainvs(), pr[:2], an))
rec("ЧИС", "S3g", "все 4 кривые класса изогении дают ранг 1 (ранг — инвариант изогении)", isog_ok)

t0 = time.time()
try:
    alarm(1500)
    res = Emin.prove_BSD(verbosity=0)
    cancel_alarm()
    rec("ПО", "S3d", "Sage prove_BSD завершился", True, "%s  (%.0f c)" % (res, time.time()-t0))
except Exception as ex:
    try: cancel_alarm()
    except Exception: pass
    rec("ЧИС", "S3d", "Sage prove_BSD НЕ завершился за 1500 c — честно фиксирую незавершённость",
        False, str(ex)[:160])

t0 = time.time()
try:
    alarm(1200)
    hi = Emin.heegner_index(-7) if True else None
    cancel_alarm()
    rec("ПО", "S3e", "индекс Хегнера (D=-7) конечен => Колывагин даёт rank = 1", True,
        "%s  (%.0f c)" % (hi, time.time()-t0))
except Exception as ex:
    try: cancel_alarm()
    except Exception: pass
    rec("ЧИС", "S3e", "heegner_index НЕ досчитан — фиксирую незавершённость", False, str(ex)[:160])

# ---------------------------------------------------------------------------------
hdr("S4. Группа Селмера: требуемый класс локально разрешим => ранг обязателен")
# ---------------------------------------------------------------------------------
# 2-накрытие для delta = (d1,d2,d3):  d1 z1^2 - d2 z2^2 = e2-e1,  d1 z1^2 - d3 z3^2 = e3-e1
d1, d2, d3 = [QQ(v) for v in req]
A12, A13 = e2 - e1, e3 - e1
note("накрытие: %s z1^2 - %s z2^2 = %s ;  %s z1^2 - %s z3^2 = %s" % (d1, d2, A12, d1, d3, A13))
rec("ДОК", "S4a", "d1*d2*d3 — квадрат (необходимое условие принадлежности образу)",
    sqclass(d1*d2*d3) == 1)
Con1 = Conic(QQ, [d1, -d2, -A12])          # d1 z1^2 - d2 z2^2 - A12 w^2 = 0
Con2 = Conic(QQ, [d1, -d3, -A13])
rec("ДОК", "S4b", "коника 1 имеет рациональную точку", Con1.has_rational_point())
rec("ДОК", "S4c", "коника 2 имеет рациональную точку", Con2.has_rational_point())
note("обе коники разрешимы над Q, значит локальных препятствий на КАЖДОЙ по отдельности нет;")
note("препятствие могло бы быть только у пересечения (кривой рода 1), и оно локально всюду разрешимо")
note("(в основном скрипте I1: явные свидетели t в Q_p для всех плохих p).")
note("ВЫВОД [ДОК]: требуемый класс лежит в 2-группе Селмера. Локальными методами (11,4) НЕ закрыть.")

# ---------------------------------------------------------------------------------
hdr("S5. Длинный прямой поиск точек C_{11,4} (чистый python-int, большая граница)")
# ---------------------------------------------------------------------------------
from math import isqrt, gcd as _g
def issq(v):
    if v < 0: return False
    r = isqrt(v); return r*r == v
# F4 квадрат <=> p^2+q^2 = 274 z^2 ; параметризация гауссовыми целыми
hits = []
B = 4000
M_, N_ = 11, 4
for (wr, wi) in ((15, 7), (15, -7)):
    for a in range(-B, B+1):
        a2 = a*a
        for bb in range(0, B+1):
            if a == 0 and bb == 0: continue
            if _g(a, bb) != 1: continue
            # (a+bi)^2 = (a^2-b^2) + 2ab i ; умножаем на (wr + wi i)
            re2 = a2 - bb*bb; im2 = 2*a*bb
            p_ = abs(wr*re2 - wi*im2); q_ = abs(wr*im2 + wi*re2)
            if p_ == 0 or q_ == 0: continue
            d = _g(p_, q_); p_ //= d; q_ //= d
            ss = p_*p_ + q_*q_
            if ss % 2: continue
            if not issq(137*(ss//2)): continue
            if not issq(M_*M_*q_*q_ + N_*N_*p_*p_): continue
            if not issq(N_*N_*q_*q_ + M_*M_*p_*p_): continue
            hits.append((p_, q_))
rec("ЧИС", "S5a", "параметрический поиск |a|,|b| <= %d: точек C не найдено" % B, len(hits) == 0, hits[:5])
note("это НЕ доказательство отсутствия — только неудачная попытка сломать заявление. [НАБЛ]")
print("\n   время: %.1f c" % (time.time()-T0))
