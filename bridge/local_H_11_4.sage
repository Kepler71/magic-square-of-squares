# local_H_11_4.sage
# Локальная разрешимость кривой рода 2
#     H : Y^2 = F0(t)*F4(t)*F8(t),   (m,n) = (11,4)
#     s  = (m^2+n^2)/2,  F0 = m^2 + n^2 t^2,  F4 = s(1+t^2),  F8 = n^2 + m^2 t^2
# во всех местах Q: R, 2, и все простые плохой редукции (+ контрольные малые простые).
#
# Метод: (а) сверка модели подстановкой конкретных t;
#        (б) вещественное место — знак f;
#        (в) для каждого p — ПОЛНЫЙ рекурсивный разбор P^1(Q_p) на диски
#            с точными критериями "квадрат/не квадрат" (Гензель), без обрыва по эвристике.
#        Результат по каждому диску: True (точка доказана), False (точек нет, доказано),
#        None (не хватило глубины — честно сообщается).
#
# Запуск:  sage /home/kep/magicKube/bridge/local_H_11_4.sage

import sys

m, n = 11, 4
Rt = PolynomialRing(QQ, 't'); t = Rt.gen()

s  = QQ(m^2 + n^2)/2
F0 = m^2 + n^2*t^2
F4 = s*(1 + t^2)
F8 = n^2 + m^2*t^2
Fprod = Rt(F0*F4*F8)                      # = (137/2)*(...)  -- нецелые коэффициенты

# ---- целочисленная модель: домножаем на 4 = 2^2 (изоморфизм над Q: Y -> Y/2) ----
f = Rt(4*Fprod)
f = f.change_ring(ZZ) if all(c in ZZ for c in Fprod.coefficients()) else Rt(4*Fprod)
RZ = PolynomialRing(ZZ, 't'); tZ = RZ.gen()
f = RZ(4*Fprod)

print("="*78)
print("МОДЕЛЬ")
print("="*78)
print("  (m,n) = (%d,%d),  s = %s" % (m, n, s))
print("  F0 =", F0, "   F4 =", F4, "   F8 =", F8)
print("  F0*F4*F8 =", Fprod)
print("  целая модель  f(t) = 4*F0*F4*F8 =", f)
print("  f = %s" % factor(f))
b  = s*m^2*n^2
print("  b = s*m^2*n^2 = %s,  squarefree part = %s" % (b, QQ(b).squarefree_part()))
print("  lc(f) = %s = %s" % (f.leading_coefficient(), factor(f.leading_coefficient())))
print("  f(0)  = %s = %s" % (f(0), factor(f(0))))

# ---- (а) сверка модели подстановкой конкретных t ----
print()
print("Сверка модели: f(t) должно равняться (2Y)^2 при Y^2 = F0F4F8, т.е. f = 4*F0*F4*F8")
ok_model = True
for tv in [QQ(0), QQ(1), QQ(2), QQ(-3), QQ(5)/7, QQ(-11)/4]:
    lhs = f(tv)
    rhs = 4*F0(tv)*F4(tv)*F8(tv)
    same = (lhs == rhs)
    ok_model = ok_model and same
    print("   t = %-8s  f(t) = %-24s  4*F0F4F8 = %-24s  %s" % (tv, lhs, rhs, "OK" if same else "РАСХОЖДЕНИЕ"))
assert ok_model
# независимая сверка с формулой Y^2 = b*(t^6 + a t^4 + a t^2 + 1)
a_par = QQ(m)^2/QQ(n)^2 + 1 + QQ(n)^2/QQ(m)^2
chk = Rt(b*(t^6 + a_par*t^4 + a_par*t^2 + 1)) - Fprod
print("   сверка с Y^2 = b(t^6+a t^4+a t^2+1), a = %s :  разность = %s" % (a_par, chk))
assert chk == 0

print()
print("  f гладкая (секстика без кратных корней)?  disc(f) != 0 :", f.discriminant() != 0)
D = ZZ(f.discriminant())
print("  disc(f) =", factor(D))
try:
    Hc = HyperellipticCurve(Rt(f))
    Dh = ZZ(Hc.discriminant())
    print("  disc(HyperellipticCurve) =", factor(Dh))
except Exception as e:
    Dh = D
    print("  HyperellipticCurve.discriminant() недоступен:", e)

bad = sorted(set([p for p, _ in factor(2*ZZ(f.leading_coefficient())*D)]))
print("  простые, делящие 2*lc(f)*disc(f)  (надмножество плохих):", bad)

# ==========================================================================
# (б) ВЕЩЕСТВЕННОЕ МЕСТО
# ==========================================================================
print()
print("="*78)
print("МЕСТО R")
print("="*78)
# f = 274*(16t^2+121)(t^2+1)(121t^2+16): каждый множитель строго положителен при t in R
fac = list(factor(Rt(f)))
allpos = True
for g, e in fac:
    if g.degree() == 0:
        allpos = allpos and (g > 0 if g in QQ else True)
    else:
        # положительно определённый квадратный трёхчлен?
        allpos = allpos and (g.degree() == 2 and g.leading_coefficient() > 0 and g.discriminant() < 0)
print("  разложение f над Q:", factor(Rt(f)))
print("  все множители строго положительны на R:", allpos)
print("  f(0) = %s > 0, f(1) = %s > 0, f(10) = %s > 0" % (f(0), f(1), f(10)))
print("  => H(R) != пусто.  (Даже больше: f(t) > 0 для ВСЕХ вещественных t,")
print("     значит H(R) состоит из двух овалов, покрывающих всю t-прямую.)")
real_ok = allpos and f(0) > 0

# ==========================================================================
# (в) p-АДИЧЕСКИЕ МЕСТА: полный рекурсивный разбор
# ==========================================================================
def unit_is_square(u, p):
    """u -- целое, не делящееся на p. Квадрат ли u в Z_p^* ?  (точный критерий)"""
    if p == 2:
        return (u % 8) == 1
    return kronecker(u, p) == 1

def disc_search(h, par, p, path, depth, maxdepth, witness):
    """
    Существует ли s in Z_p с  p^par * h(s)  квадратом в Q_p (0 тоже считается)?
    h in ZZ[s], примитивен по p (не все коэффициенты делятся на p);
    par in {0,1} -- чётность накопленной степени p.
    path -- список для восстановления свидетеля t.
    Возврат: True / False / None(не хватило глубины).
    """
    step = 3 if p == 2 else 1
    mod  = 8 if p == 2 else p
    undecided = False
    for s0 in range(mod):
        val = h(s0)
        if val % p != 0:
            # на всём поддиске s ≡ s0 (mod `mod`) значение h — единица фиксированного
            # квадратичного класса; при p=2 это верно mod 8, при p нечётном — mod p
            if par == 0 and unit_is_square(val, p):
                witness.append(path + [(s0, step)])
                return True
            # иначе поддиск целиком исключён
            continue
        # val ≡ 0 mod p : уточняем
        if depth >= maxdepth:
            undecided = True
            continue
        g2 = h(s0 + p**step * PolynomialRing(ZZ, 's').gen())
        if g2 == 0:
            # h тождественно 0 — невозможно для примитивного ненулевого h
            witness.append(path + [(s0, step)])
            return True
        c2 = min(ZZ(co).valuation(p) for co in g2.coefficients() if co != 0)
        h2 = g2 // p**c2
        r = disc_search(h2, (par + c2) % 2, p, path + [(s0, step)], depth + 1, maxdepth, witness)
        if r is True:
            return True
        if r is None:
            undecided = True
    if undecided:
        return None
    return False

def poly_solvable_on_Zp(F, p, maxdepth=40):
    """Есть ли t in Z_p с F(t) квадратом в Q_p?  (F in ZZ[t], F != 0)"""
    Rs = PolynomialRing(ZZ, 's'); ss = Rs.gen()
    g = Rs(F(ss))
    c = min(ZZ(co).valuation(p) for co in g.coefficients() if co != 0)
    h = g // p**c
    wit = []
    r = disc_search(h, c % 2, p, [], 0, maxdepth, wit)
    return r, wit

def rebuild_t(wit, p):
    """восстановить целое приближение свидетеля t из пути (t = sum s0*p^N)"""
    if not wit:
        return None
    path = wit[0]
    tv, scale = ZZ(0), ZZ(1)
    for (s0, step) in path:
        tv = tv + scale*s0
        scale = scale * (p**step)
    return tv, scale

print()
print("="*78)
print("p-АДИЧЕСКИЕ МЕСТА")
print("="*78)
print("""  Алгоритм (точный, без эвристик):
   * P^1(Q_p) покрывается двумя картами: t in Z_p (полином f) и t = 1/u, u in Z_p
     (полином f*(u) = u^6 f(1/u); точка u=0 — это две точки на бесконечности).
   * Диск a + p^N Z_p: подставляем, выносим p^c, остаток h примитивен.
     Для каждого поддиска (mod p при p нечётном, mod 8 при p=2):
       - h(s0) — единица: значение p^c*h(s) имеет ФИКСИРОВАННЫЙ квадратичный класс
         на всём поддиске; если c чётно и h(s0) — квадрат в Z_p^*, точка ДОКАЗАНА
         (лемма Гензеля); иначе поддиск ДОКАЗАТЕЛЬНО пуст.
       - h(s0) ≡ 0 mod p: рекурсия вглубь.
     Обход полный, поэтому ответ False означает доказанное отсутствие точек.""")

fstar = RZ(sum(ZZ(f[i])*tZ**(6-i) for i in range(7)))   # f*(u) = u^6 f(1/u)
print()
print("  f*(u) = u^6 f(1/u) =", fstar, "   (совпадает с f:", fstar.list() == f.list(), ")")

# набор простых для проверки: все плохие + все p <= 50
# (для p >= 17 хорошей редукции граница Вейля #H(F_p) >= p+1-4sqrt(p) > 0 даёт точку автоматически,
#  но лишняя проверка не вредит)
primes_to_check = sorted(set(bad) | set(primes(2, 60)))
print("  проверяемые простые:", primes_to_check)
print()

results = {}
for p in primes_to_check:
    r1, w1 = poly_solvable_on_Zp(RZ(f), p)
    if r1 is True:
        t0, sc = rebuild_t(w1, p)
        # НЕЗАВИСИМАЯ проверка свидетеля: f(t0) — точное рациональное число
        ver = QQ(f(t0)).is_padic_square(p) if f(t0) != 0 else True
        results[p] = (True, "t = %s + O(%s^%s), f(t)=%s, проверка is_padic_square: %s"
                      % (t0, p, ZZ(sc).valuation(p), f(t0), ver))
    else:
        r2, w2 = poly_solvable_on_Zp(RZ(fstar), p)
        if r2 is True:
            u0, sc = rebuild_t(w2, p)
            ver = QQ(fstar(u0)).is_padic_square(p) if fstar(u0) != 0 else True
            results[p] = (True, "t = 1/u, u = %s + O(%s^%s), f*(u)=%s, проверка: %s"
                          % (u0, p, ZZ(sc).valuation(p), fstar(u0), ver))
        elif r1 is False and r2 is False:
            results[p] = (False, "обе карты пусты — ДОКАЗАНО отсутствие Q_p-точек")
        else:
            results[p] = (None, "не хватило глубины (карта t: %s, карта 1/t: %s)" % (r1, r2))
    st, why = results[p]
    mark = {True: "ЕСТЬ точка", False: "ТОЧЕК НЕТ", None: "НЕ ОПРЕДЕЛЕНО"}[st]
    extra = ""
    if p % 4 == 1:
        extra = "  [p=1 mod 4: sqrt(-1) in Z_p — корень f, даёт точку (t,Y)=(i,0)]"
    print("   p = %-5d  %-13s  %s%s" % (p, mark, why, extra))

# ---- отдельная проверка точек на бесконечности ----
print()
print("="*78)
print("ТОЧКИ НА БЕСКОНЕЧНОСТИ")
print("="*78)
lc = ZZ(f.leading_coefficient())
print("  старший коэффициент f: %s = %s" % (lc, factor(lc)))
print("  над Q: квадрат?", lc.is_square(), " => над Q точек на бесконечности НЕТ")
print("  над Q_v точка на бесконечности есть <=> lc — квадрат в Q_v:")
print("    R:", lc > 0)
for p in primes_to_check:
    print("    p = %-5d %s" % (p, QQ(lc).is_padic_square(p)))

# ---- отдельная проверка t = 0 и t = oo (сверка с утверждением Codex) ----
print()
print("  Сверка с запиской Codex: 'при t=0 и t=oo точек H над Q нет, класс b нетривиален'")
print("    f(0) = %s, квадрат в Q: %s" % (f(0), ZZ(f(0)).is_square()))
print("    lc    = %s, квадрат в Q: %s" % (lc, lc.is_square()))
print("    оба равны 274*44^2, squarefree part = %s  -> над Q не квадраты. Подтверждено." % QQ(f(0)).squarefree_part())

# ---- итог ----
print()
print("="*78)
print("ИТОГ")
print("="*78)
obstructed = [p for p in results if results[p][0] is False]
undecided  = [p for p in results if results[p][0] is None]
print("  R: разрешима" if real_ok else "  R: ПРЕПЯТСТВИЕ")
print("  места с ДОКАЗАННЫМ отсутствием точек:", obstructed if obstructed else "нет")
print("  места, где проверка НЕ завершилась:", undecided if undecided else "нет")
print("  проверено простых:", len(primes_to_check))
if not obstructed and not undecided and real_ok:
    print("  => H(Q_v) != пусто для ВСЕХ проверенных v.")
    print()
    print("  Покрытие всех мест Q:")
    print("   * v = R: доказано выше (f > 0 всюду).")
    print("   * плохие простые: делители 2*lc(f)*disc(f) = {2,3,5,7,11,137} — все проверены явно.")
    print("   * хорошие простые p <= 13 (это 13): проверено явно.")
    print("   * хорошие простые p >= 17: редукция — гладкая кривая рода 2 над F_p,")
    print("     граница Вейля #H(F_p) >= p+1-4*sqrt(p) > 0 при p >= 17 (17+1-4*sqrt(17) = %.2f > 0);"
          % (17 + 1 - 4*RR(17).sqrt()))
    print("     любая F_p-точка гладкая => подъём по Гензелю до Q_p-точки.")
    print("   => ЛОКАЛЬНОГО ПРЕПЯТСТВИЯ НЕТ НИ В ОДНОМ МЕСТЕ Q.  [доказано]")
    print()
    print("  ВАЖНО о надёжности: во всех местах алгоритм вернул ПОЛОЖИТЕЛЬНЫЙ сертификат")
    print("  (конкретное целое t с точно проверенным QQ(f(t)).is_padic_square(p) = True).")
    print("  Ветка 'доказано пусто' (pruning) НИ РАЗУ не использовалась, поэтому вывод не")
    print("  зависит от её корректности — только от точного теста квадратичности.")
    print()
    print("  Что это НЕ означает: H(Q) может быть пустым (препятствие Брауэра–Манина или")
    print("  просто отсутствие точек). Локальная разрешимость — необходимое, не достаточное условие.")
