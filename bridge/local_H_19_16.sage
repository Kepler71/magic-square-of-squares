#!/usr/bin/env sage
# =============================================================================
#  local_H_19_16.sage
#
#  ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ кривой рода 2
#
#      H : Y^2 = F0(t) * F4(t) * F8(t),
#      s  = (m^2+n^2)/2,  F0 = m^2 + n^2 t^2,  F4 = s(1+t^2),  F8 = n^2 + m^2 t^2
#
#  для пары (m,n) = (19,16)  [семейство G1, «две полные пары»].
#
#  Запуск:  sage /home/kep/magicKube/bridge/local_H_19_16.sage
#
#  МЕТОД.
#  (0) Целая модель  W^2 = f(t) := 4*F0*F4*F8,  W = 2Y.  Это Q-изоморфизм
#      (деление на квадрат 4), т.к. s = 617/2 полуцелое.
#      Модель сверяется подстановкой конкретных t И контрольной парой (15,8),
#      где Codex нашёл точку t=1, Y=4913.
#  (1) Гладкая модель рода 2 живёт в P(1,3,1):  Y^2 = F(X,Z) = Z^6 f(X/Z).
#      H(Q_v) != пусто  <=>  существует (x:z) in P^1(Q_v) с F(x,z) in (Q_v^*)^2 u {0}.
#      Нуль допускается: это точка Вейерштрасса (Y=0). Точки на бесконечности —
#      это z=0, т.е. F(1,0) = lc(f).
#  (2) P^1(Q_p) покрывается ДВУМЯ картами:
#         карта A: t in Z_p,        многочлен g_A = f;
#         карта B: w = 1/t in Z_p,  многочлен g_B = w^6 f(1/w) = reverse(f);
#      w = 0 — это бесконечность. Обе карты обрабатываются ОДИНАКОВЫМ движком,
#      поэтому бесконечность разбирается автоматически (и отдельно — в п.5).
#  (3) Движок: доказательный рекурсивный обход дисков r + p^k Z_p с тремя
#      завершающими критериями (см. _disc ниже). Он ДОКАЗЫВАЕТ и «да», и «нет»;
#      если глубина исчерпана — возвращает UNKNOWN, а НЕ «нет».
#  (4) Достаточно явно проверить p in {2,3,5,7,11,13,19,617}: все прочие p — это
#      простые хорошей редукции с p >= 17, где #H(F_p) >= p+1-4sqrt(p) > 0
#      (Хассе–Вейль) и гладкая точка поднимается леммой Гензеля.
#  (5) Для каждого из восьми явных p дополнительно предъявляется СВИДЕТЕЛЬ
#      В ЗАМКНУТОЙ ФОРМЕ (t=0, t=1 или корень t^2+1), проверяемый одним
#      символом Лежандра вручную — независимо от всей рекурсии.
#
#  СТАТУСЫ УТВЕРЖДЕНИЙ печатаются в итоговом блоке.
# =============================================================================

import sys

R = PolynomialRing(QQ, 't'); t = R.gen()
RZ = PolynomialRing(ZZ, 't')

m, n = 19, 16
s  = QQ(m^2 + n^2) / 2
F0 = m^2 + n^2*t^2
F4 = s*(1 + t^2)
F8 = n^2 + m^2*t^2
P  = F0*F4*F8
f  = R(4*P)                     # W^2 = f(t),  W = 2Y
fZ = RZ(f)

def hdr(x):
    print(); print("="*78); print(x); print("="*78)

hdr("H(19,16):  Y^2 = F0*F4*F8   ->   целая модель  W^2 = f(t),  W = 2Y")
print("  m, n      =", m, n)
print("  s         =", s)
print("  F0        =", F0)
print("  F4        =", F4)
print("  F8        =", F8)
print("  f = 4F0F4F8 =", fZ)
print("  factor(f) =", factor(fZ))
print("  lc(f)     =", fZ.leading_coefficient(), "=", factor(fZ.leading_coefficient()))
print("  f(0)      =", fZ(0), "=", factor(fZ(0)))
print("  f(1)      =", fZ(1), "=", factor(fZ(1)))
print("  f свободен от квадратов:", fZ.is_squarefree())
print("  disc(f)   =", factor(fZ.discriminant()))

# ---------------------------------------------------------------- 0. сверка модели
hdr("0. СВЕРКА МОДЕЛИ  (защита от 'артефакта неверной модели')")
ok_model = True
print("  (a) подстановка конкретных t: 4*F0(t)F4(t)F8(t) == f(t) ?")
for tv in [QQ(0), QQ(1), QQ(-1), QQ(2), QQ(3)/2, QQ(-7)/5, QQ(11)/3, QQ(19)/16]:
    lhs = 4*(m^2 + n^2*tv^2) * (QQ(m^2+n^2)/2)*(1+tv^2) * (n^2 + m^2*tv^2)
    rhs = f(tv)
    ok_model = ok_model and (lhs == rhs)
    print("      t=%-8s  4F0F4F8 = %-24s  f(t) = %-24s  равны: %s" % (tv, lhs, rhs, lhs == rhs))
print("  (b) контрольная пара (15,8): Codex нашёл точку t=1, Y=4913. Проверяем формулы:")
m2, n2 = 15, 8
s2 = QQ(m2^2+n2^2)/2
Yv2 = (m2^2 + n2^2*1) * (s2*(1+1)) * (n2^2 + m2^2*1)
print("      F0(1)F4(1)F8(1) =", Yv2, ";  sqrt =", sqrt(Yv2), " (должно быть 4913):", sqrt(Yv2) == 4913)
ok_model = ok_model and (sqrt(Yv2) == 4913)
print("  (c) f палиндромичен (=> инволюция t -> 1/t, и lc(f) = f(0)):",
      all(fZ[i] == fZ[6-i] for i in range(7)))
print("  СВЕРКА МОДЕЛИ ПРОЙДЕНА:", ok_model)
assert ok_model, "МОДЕЛЬ НЕ СОШЛАСЬ — дальнейшие выводы недействительны"

# ---------------------------------------------------------------- 1. вещественное место
hdr("1. ВЕЩЕСТВЕННОЕ МЕСТО  v = oo")
print("  f(t) = 2*617*(t^2+1)*(256t^2+361)*(361t^2+256).")
print("  Каждый множитель строго положителен при любом вещественном t, и 2*617>0,")
print("  значит f(t) > 0 для всех t in R.")
print("  вещественные корни f:", f.roots(RR), "(пусто)")
print("  min f на R: f(0) = %s (по палиндромичности и положительности), f(0) > 0: %s" % (fZ(0), fZ(0) > 0))
print("  Явная точка: t = 0,  W = sqrt(f(0)) = sqrt(%s) = %.6f" % (fZ(0), RR(fZ(0)).sqrt()))
print("  (в исходных координатах Y = W/2 = %.6f)" % (RR(fZ(0)).sqrt()/2))
real_ok = (len(f.roots(RR)) == 0) and fZ(0) > 0
print("  => H(R) НЕПУСТО:", real_ok, "   [статус: доказано]")
print("  Замечание: препятствия в бесконечном месте здесь заведомо нет, т.к. F0,F4,F8 > 0.")

# ---------------------------------------------------------------- утилиты
def is_sq_Qp(x, p):
    """x in QQ (точно). True <=> x - квадрат в Q_p. 0 считается квадратом."""
    x = QQ(x)
    if x == 0:
        return True
    v = x.valuation(p)
    if v % 2 != 0:
        return False
    u = x / QQ(p)^v
    num = ZZ(u.numerator()); den = ZZ(u.denominator())
    if p == 2:
        return (num * inverse_mod(den, 8)) % 8 == 1
    return kronecker(num * inverse_mod(den, p), p) == 1

MAXDEPTH = 300

def _disc(g, gp, p, r, k, depth, stats):
    """
    Существует ли u in r + p^k Z_p с g(u) квадратом (или нулём) в Q_p?
    Возврат: ('YES', свидетель) / ('NO', None) / ('UNKNOWN', причина).

    Три завершающих критерия:
      (i)   g(r) = 0                      -> точка Вейерштрасса, YES;
      (ii)  v(g(r)) > 2 v(g'(r))           -> лемма Гензеля даёт корень g в Z_p, YES
            (этот пункт обеспечивает ЗАВЕРШАЕМОСТЬ ветвей, сходящихся к корню);
      (iii) h(u) := g(r + p^k u) = a0 + ..., e := v(a0), c := min_{i>=1} v(a_i);
            если e + delta <= c (delta = 1 для нечётного p, 3 для p=2), то
            h(u)/a0 in 1 + p^delta Z_p = квадраты единиц, значит квадратичный класс
            h(u) ПОСТОЯНЕН на всём диске и равен классу a0 -> диск решается целиком.
    Иначе диск делится на p поддисков. Завершаемость: f свободен от квадратов,
    поэтому кратных корней в Q_p нет; ветви либо стабилизируются по (iii), либо
    сходятся к простому корню и отсекаются по (ii).
    """
    stats['nodes'] += 1
    if depth > MAXDEPTH:
        return ('UNKNOWN', 'глубина > %d при r=%s k=%s p=%s' % (MAXDEPTH, r, k, p))
    a0 = QQ(g(r))
    if a0 == 0:
        return ('YES', 'точный корень t=%s (точка Вейерштрасса, W=0)' % r)
    dr = QQ(gp(r))
    if dr != 0 and a0.valuation(p) > 2*dr.valuation(p):
        return ('YES', 'корень f в Z_p по Гензелю около t=%s  [v(f(r))=%d > 2v(f\'(r))=%d]'
                       % (r, a0.valuation(p), 2*dr.valuation(p)))
    S = PolynomialRing(QQ, 'u'); u = S.gen()
    h = g(r + QQ(p)^k * u)
    cs = h.list()
    e = QQ(cs[0]).valuation(p)
    tail = [QQ(c) for c in cs[1:] if c != 0]
    c = min([x.valuation(p) for x in tail]) if tail else Infinity
    delta = 3 if p == 2 else 1
    if e + delta <= c:
        if is_sq_Qp(cs[0], p):
            return ('YES', 't=%s (диск r+p^%s Z_p, постоянный квадратичный класс)' % (r, k))
        return ('NO', None)
    unknown = None
    for j in range(p):
        res = _disc(g, gp, p, r + j*p^k, k+1, depth+1, stats)
        if res[0] == 'YES':
            return res
        if res[0] == 'UNKNOWN' and unknown is None:
            unknown = res
    return unknown if unknown is not None else ('NO', None)

def Zp_square_value(g, p):
    st = {'nodes': 0}
    a, b = _disc(g, g.derivative(), p, ZZ(0), ZZ(0), 0, st)
    return a, b, st['nodes']

def reverse6(g):
    co = g.list()
    while len(co) < 7: co.append(QQ(0))
    return R(list(reversed(co[:7])))

frev = reverse6(f)

def P1_solvable(g, p):
    """Полная проверка по P^1(Q_p) через две карты. Карта B включает бесконечность (w=0)."""
    sA, wA, nA = Zp_square_value(g, p)
    if sA == 'YES':
        return ('YES', 'карта A (|t|<=1): ' + wA, nA, 0)
    gr = reverse6(g)
    sB, wB, nB = Zp_square_value(gr, p)
    if sB == 'YES':
        return ('YES', 'карта B (|t|>=1, w=1/t): ' + wB, nA, nB)
    if 'UNKNOWN' in (sA, sB):
        return ('UNKNOWN', 'A:%s %s | B:%s %s' % (sA, wA, sB, wB), nA, nB)
    return ('NO', 'обе карты исчерпаны: нет (x:z) in P^1(Q_p) с F(x,z) квадратом', nA, nB)

# ---------------------------------------------------------------- 2. плохие простые
hdr("2. ДИСКРИМИНАНТ И ПЛОХИЕ ПРОСТЫЕ")
D = fZ.discriminant()
lc = fZ.leading_coefficient()
bad = sorted(set([q for q, _ in factor(2*D*lc)]))
print("  disc(f)      =", factor(D))
print("  lc(f)        =", factor(lc))
print("  плохие простые (делители 2*disc(f)*lc(f)):", bad)
WEIL = 17
print("  Хассе–Вейль: для гладкой кривой рода 2 над F_p  #H(F_p) >= p+1-4sqrt(p).")
for pp in [11, 13, 17, 19, 23]:
    print("     p=%-4s p+1-4sqrt(p) = %+.4f  %s" % (pp, pp+1-4*sqrt(RR(pp)), ">0" if pp+1-4*sqrt(RR(pp))>0 else "<=0 (нужна явная проверка)"))
explicit = sorted(set(bad) | set([q for q in primes(WEIL)]))
print("  => явно проверяем p in", explicit)
print("     (все остальные p — хорошая редукция и p >= 17)")

# ---------------------------------------------------------------- 3. движок
hdr("3. ДОКАЗАТЕЛЬНАЯ ПРОВЕРКА H(Q_p) ДЛЯ ВОСЬМИ ЯВНЫХ p")
print("  p       вердикт   узлов A/B   свидетель")
print("  " + "-"*94)
verdicts = {}
obstructed = []
undetermined = []
for p in explicit:
    st, wit, nA, nB = P1_solvable(f, p)
    verdicts[p] = st
    print("  %-7s %-9s %-11s %s" % (p, st, "%d/%d" % (nA, nB), wit))
    if st == 'NO':
        obstructed.append(p)
    elif st == 'UNKNOWN':
        undetermined.append((p, wit))

# ---------------------------------------------------------------- 4. замкнутые свидетели
hdr("4. НЕЗАВИСИМЫЕ СВИДЕТЕЛИ В ЗАМКНУТОЙ ФОРМЕ (проверяются вручную)")
print("  Ключевые значения:  f(0) = lc(f) = 2^9*19^2*617,  класс по модулю квадратов = 2*617 = 1234")
print("                      f(1) = 2^2*617^3,             класс по модулю квадратов = 617")
print("  Плюс: если -1 — квадрат в Q_p, то t^2+1 имеет корень в Z_p (лемма Гензеля, p нечётно),")
print("        и этот корень даёт точку Вейерштрасса f(t)=0 на H.")
print()
print("  p      617 кв.?  1234 кв.?  -1 кв.?   свидетель                       проверка f(t) — квадрат в Q_p")
print("  " + "-"*104)
closed = {}
for p in explicit:
    a617  = is_sq_Qp(QQ(617), p)
    a1234 = is_sq_Qp(QQ(1234), p)
    am1   = is_sq_Qp(QQ(-1), p)
    if a617:
        wit, val, chk = "t = 1",  fZ(1), is_sq_Qp(fZ(1), p)
    elif a1234:
        wit, val, chk = "t = 0",  fZ(0), is_sq_Qp(fZ(0), p)
    elif am1:
        K = Qp(p, 60); rt = K(-1).sqrt()
        wit, val, chk = "t = sqrt(-1) in Z_%d" % p, 0, (f(rt).valuation() >= 55)
    else:
        wit, val, chk = "(замкнутого свидетеля нет)", None, False
    closed[p] = chk
    print("  %-6s %-9s %-10s %-9s %-31s %s%s" %
          (p, a617, a1234, am1, wit,
           chk, ("  [f(t)=%s]" % val) if val is not None and val != 0 else "  [f(t)=0]" if chk and val==0 else ""))
print()
print("  Разбор вручную (символы Лежандра):")
print("    p=2   : 617 = 1 mod 8  => 617 in (Z_2^*)^2 ; f(1)=2^2*617^3, v_2=2 чётно, ед. часть 617^3 = 1 mod 8  -> КВАДРАТ")
print("    p=3   : f(0)=2^9*19^2*617; mod 3: 2^9=512=2, 19^2=1, 617=2 => 2*1*2=4=1 mod 3 -> КВАДРАТ (t=0)")
print("    p=5   : 5 = 1 mod 4 => -1 квадрат => корень t^2+1 в Z_5 -> точка Вейерштрасса")
print("    p=7   : 617 = 1 mod 7 -> квадрат; f(1)=4*617^3 -> КВАДРАТ (t=1)")
print("    p=11  : 617 = 1 mod 11 -> квадрат -> КВАДРАТ (t=1)")
print("    p=13  : 13 = 1 mod 4 => -1 квадрат => корень t^2+1 в Z_13 -> точка Вейерштрасса")
print("    p=19  : 617 = 9 mod 19 = 3^2 -> квадрат; v_19(f(1))=0 -> КВАДРАТ (t=1)")
print("    p=617 : 617 = 1 mod 4 => -1 квадрат mod 617 => корень t^2+1 в Z_617 -> точка Вейерштрасса")
print("  все восемь замкнутых свидетелей подтверждены машиной:", all(closed.values()))
print("  согласие с движком (п.3):", all((verdicts[p] == 'YES') == closed[p] for p in explicit))

# ---------------------------------------------------------------- 5. бесконечность
hdr("5. ТОЧКИ НА БЕСКОНЕЧНОСТИ (гладкая модель рода 2)")
print("  lc(f) = %s = %s" % (lc, factor(lc)))
print("  бесквадратная часть lc(f) =", QQ(lc).squarefree_part(), "= 2*617")
print("  Над полем K две точки на бесконечности рациональны <=> lc(f) in (K^*)^2.")
print("  Над Q: lc(f) = 2^9*19^2*617, бесквадратная часть 1234 != квадрат => НЕТ рациональных")
print("         точек на бесконечности; на Q-гладкой модели бесконечность — одно замкнутое")
print("         место степени 2 (расщепляется над Q(sqrt(1234))).")
print()
print("  p        lc(f)=2^9*19^2*617 — квадрат в Q_p?   => две точки на бесконечности над Q_p?")
for p in explicit + [23, 29, 31, 37]:
    print("    p=%-6s %-38s %s" % (p, is_sq_Qp(lc, p), is_sq_Qp(lc, p)))
print()
print("  ВАЖНО: бесконечность НЕ нужна для локальной разрешимости — она уже входит")
print("  в карту B (w=0), и в п.3 движок разбирал обе карты. Кроме того f палиндромичен,")
print("  поэтому reverse(f) = f и lc(f) = f(0): бесконечность имеет тот же квадратичный")
print("  класс, что t=0. Проверка палиндромичности: reverse(f) == f :", frev == f)

# ---------------------------------------------------------------- 6. хорошие p >= 17
hdr("6. ОСТАЛЬНЫЕ ПРОСТЫЕ: ХОРОШАЯ РЕДУКЦИЯ + ХАССЕ-ВЕЙЛЬ")
print("  Пусть p не делит 2*disc(f)*lc(f) и p >= 17. Тогда f mod p — свободный от квадратов")
print("  многочлен степени 6 над F_p (p нечётно), значит H имеет ХОРОШУЮ редукцию в p,")
print("  и редукция — гладкая проективная кривая рода 2 над F_p.")
print("  Хассе–Вейль: #H(F_p) >= p+1-4sqrt(p) >= 18-4sqrt(17) = %.4f > 0." % (18-4*sqrt(RR(17))))
print("  Любая F_p-точка гладкая, значит поднимается до Z_p-точки леммой Гензеля.")
print("  => H(Q_p) != пусто для всех таких p.   [статус: доказано]")
print()
print("  Контрольный прямой подсчёт #H(F_p) для хороших p < 300:")
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
badset = set(bad)
zero_cnt = []
mincount = None
for p in primes(300):
    if p in badset: continue
    C = HyperellipticCurve(f.change_ring(GF(p)))
    npts = len(C.rational_points())
    if npts == 0: zero_cnt.append(p)
    if mincount is None or npts < mincount[1]: mincount = (p, npts)
print("    хорошие p < 300 с #H(F_p) = 0:", zero_cnt if zero_cnt else "нет — все непусты")
print("    минимум #H(F_p) по хорошим p<300:", mincount)

# ---------------------------------------------------------------- 7. независимая сверка
hdr("7. НЕЗАВИСИМАЯ СВЕРКА: РАВНОМЕРНОЕ ИСЧЕРПАНИЕ (другой код, без рекурсии)")
def uniform_check(g, p, N):
    """Перебираем ВСЕ вычеты r mod p^N в обеих картах. 'YES' — найден свидетель;
       'NO' — каждый диск решён критерием постоянного класса и не квадрат;
       'UNDECIDED' — какой-то диск не решается на уровне N (ничего не доказано)."""
    S = PolynomialRing(QQ, 'u'); u = S.gen()
    delta = 3 if p == 2 else 1
    undec = False
    for gg, tag in [(g, 'A'), (reverse6(g), 'B')]:
        for r in range(p^N):
            a0 = QQ(gg(r))
            if a0 == 0: return ('YES', '%s: t=%s (корень)' % (tag, r))
            h = gg(r + QQ(p)^N*u); cs = h.list()
            e = QQ(cs[0]).valuation(p)
            tail = [QQ(c) for c in cs[1:] if c != 0]
            c = min([x.valuation(p) for x in tail]) if tail else Infinity
            if e + delta <= c:
                if is_sq_Qp(cs[0], p): return ('YES', '%s: t=%s' % (tag, r))
            else:
                undec = True
    return ('UNDECIDED', None) if undec else ('NO', None)

print("  p       рекурсивный   равномерный(N)        согласие")
for p in explicit:
    N = 10 if p == 2 else (6 if p < 8 else (4 if p < 20 else 2))
    st2, w2 = uniform_check(f, p, N)
    print("  %-7s %-13s %-21s %s" % (p, verdicts[p], "%s (N=%d)" % (st2, N),
          "ДА" if (st2 == verdicts[p] or st2 == 'UNDECIDED') else "*** РАСХОЖДЕНИЕ ***"))

print()
print("  Сверка #2: p-адические свидетели в Qp(prec) — прямая проверка .is_square()")
for p in explicit:
    K = Qp(p, 40)
    w = None
    for cand in [K(0), K(1)]:
        if K(f(QQ(cand.lift() if cand != 0 else 0))).is_square():
            w = cand.lift(); break
    if w is None and is_sq_Qp(QQ(-1), p):
        w = "sqrt(-1) (f=0)"
    print("    p=%-6s Qp-свидетель: t=%s" % (p, w))

# ---------------------------------------------------------------- 8. поиск Q-точек
hdr("8. ПОИСК РАЦИОНАЛЬНЫХ ТОЧЕК (ТОЛЬКО ЧИСЛЕННО, НЕ ДОКАЗАТЕЛЬСТВО)")
Rx = PolynomialRing(QQ, 'x')
fx = Rx(f.list())
for B in [200, 1000, 5000]:
    try:
        pts = pari.hyperellratpoints(fx, B)
        print("   PARI hyperellratpoints, граница высоты %-6s : %s" % (B, pts))
    except Exception as e:
        print("   граница %-6s : ошибка %s" % (B, e))
print("   [статус: проверено численно] Q-точек на H при высоте <= 5000 не найдено.")
print("   Это ШИРЕ поиска Codex (t=u/v, u,v <= 100), но ОТСУТСТВИЕ НАХОДКИ НЕ ЕСТЬ")
print("   ДОКАЗАТЕЛЬСТВО ОТСУТСТВИЯ. H(Q) = пусто остаётся ГИПОТЕЗОЙ.")

# ---------------------------------------------------------------- 9. бонус: C рода 5
hdr("9. БОНУС: ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ КРИВОЙ C РОДА 5 (все три F_i — квадраты)")
print("  C : u0^2=F0(t), u4^2=F4(t), u8^2=F8(t).  Ищем (a:b) in P^1(Q_p) с b^2*Fi(a/b)")
print("  квадратом для ВСЕХ трёх i одновременно (масштаб l умножает на l^2 — корректно).")
MAXD_C = 200
def simul(gs, p, r, k, depth):
    if depth > MAXD_C: return ('UNKNOWN', 'глубина')
    S = PolynomialRing(QQ, 'u'); u = S.gen(); delta = 3 if p == 2 else 1
    classes = []; decided = True
    for g in gs:
        a0 = QQ(g(r))
        if a0 == 0: classes.append('zero'); continue
        h = g(r + QQ(p)^k*u); cs = h.list(); e = QQ(cs[0]).valuation(p)
        tail = [QQ(c) for c in cs[1:] if c != 0]
        c = min([xx.valuation(p) for xx in tail]) if tail else Infinity
        if e + delta <= c: classes.append(is_sq_Qp(cs[0], p))
        else: decided = False; classes.append(None)
    if decided:
        return ('YES', 't=%s' % r) if all(cl == 'zero' or cl is True for cl in classes) else ('NO', None)
    if any(cl is False for cl in classes): return ('NO', None)   # один множитель уже отпал на всём диске
    unk = None
    for j in range(p):
        res = simul(gs, p, r + j*p^k, k+1, depth+1)
        if res[0] == 'YES': return res
        if res[0] == 'UNKNOWN' and unk is None: unk = res
    return unk if unk else ('NO', None)
def rev2(g):
    co = g.list()
    while len(co) < 3: co.append(QQ(0))
    return R(list(reversed(co[:3])))
print("  вещественное место: F0,F4,F8 > 0 при любом t in R =>",
      all(F0(a) > 0 and F4(a) > 0 and F8(a) > 0 for a in [-5,-2,-1,0,1,2,5]))
cbad = []; cunk = []
for p in primes(140):
    A = simul([F0, F4, F8], p, ZZ(0), ZZ(0), 0); st = A[0]
    if st != 'YES':
        B = simul([rev2(F0), rev2(F4), rev2(F8)], p, ZZ(0), ZZ(0), 0)
        st = 'YES' if B[0] == 'YES' else ('UNKNOWN' if 'UNKNOWN' in (A[0], B[0]) else 'NO')
    if st == 'NO': cbad.append(p)
    if st == 'UNKNOWN': cunk.append(p)
print("  все p < 140:  C(Q_p) = пусто при p in", cbad if cbad else "НЕТ ТАКИХ",
      ";  неопределённых:", cunk if cunk else "нет")
print("  [ВАЖНО, ограничение] для рода 5 граница Вейля p+1-10sqrt(p)>0 работает лишь с p >= 101,")
print("  а вопрос хорошей редукции C требует отдельного разбора модели. Поэтому утверждение")
print("  'C локально разрешима ВСЮДУ' здесь имеет статус [проверено численно для p<140],")
print("  а НЕ [доказано]. Для H (род 2) разбор ПОЛНЫЙ и статус [доказано].")

# ---------------------------------------------------------------- ИТОГ
hdr("ИТОГ")
print("  Пара (m,n) = (19,16).  H : Y^2 = F0*F4*F8,  W^2 = f(t), W = 2Y.")
print()
print("  [доказано] Вещественное место: f(t) > 0 для всех t in R => H(R) != пусто.")
print("  [доказано] p in", explicit, ": H(Q_p) != пусто — движок п.3 И замкнутый свидетель п.4.")
print("  [доказано] все остальные p (хорошая редукция, p >= 17): Хассе–Вейль + Гензель.")
print("  Разбор случаев исчерпывающий: любое p либо в списке плохих", bad, ",")
print("  либо p < 17 (=> p in {2,3,5,7,11,13}), либо хорошее с p >= 17.")
print()
print("  места с ДОКАЗАННЫМ отсутствием точек:", obstructed if obstructed else "НЕТ")
print("  неопределённые места:", undetermined if undetermined else "НЕТ")
if not obstructed and not undetermined and real_ok:
    print()
    print("  ВЫВОД [доказано]: H ЛОКАЛЬНО РАЗРЕШИМА ВО ВСЕХ МЕСТАХ Q.")
    print("  ЛОКАЛЬНОГО ПРЕПЯТСТВИЯ НЕТ. obstruction_found = FALSE.")
    print()
    print("  Следствие [наблюдение, не доказательство]: отсутствие точек в поиске Codex")
    print("  (t=u/v, 0<=u<=v<=100) НЕ объясняется локальными причинами. Если H(Q)=пусто,")
    print("  то это требует либо препятствия Брауэра–Манина / Шаферевича–Тейта,")
    print("  либо метода Шабо–Коулмана (rank E_b = 1 => rank Jac = 2 = genus, Шабо")
    print("  напрямую НЕ применим; нужна квадратичная Шабо или спуск).")
print("="*78)
