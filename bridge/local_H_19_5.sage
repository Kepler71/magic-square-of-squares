# local_H_19_5.sage
#
# Локальная разрешимость кривой рода 2
#     H : Y^2 = F0(t)*F4(t)*F8(t),
#     F0 = m^2 + n^2 t^2, F4 = s(1+t^2), F8 = n^2 + m^2 t^2, s=(m^2+n^2)/2,
# для пары (m,n) = (19,5).
#
# Всё решается ТОЧНО (без p-адических приближений):
#  * квадратичность рационального числа в Q_p проверяется через v_p + символ Лежандра
#    (p нечётное) или через класс mod 8 (p=2);
#  * полный перебор t in Z_p ведётся рекурсией по дискам с сертификатами Гензеля;
#  * точка на бесконечности покрыта обращённым многочленом (f палиндромичен).
#
# Запуск:  sage /home/kep/magicKube/bridge/local_H_19_5.sage

import sys

R.<t> = QQ[]

def build(m, n):
    s = QQ(m^2 + n^2)/2
    F0 = m^2 + n^2*t^2
    F4 = s*(1 + t^2)
    F8 = n^2 + m^2*t^2
    return s, F0, F4, F8

m, n = 19, 5
s, F0, F4, F8 = build(m, n)
f = F0*F4*F8
assert f.denominator() == 1
f = f.change_ring(ZZ) if hasattr(f, 'change_ring') else f
RZ.<T> = ZZ[]
f = RZ(f)

b = s * m^2 * n^2
a = QQ(m^2)/n^2 + 1 + QQ(n^2)/m^2

print("="*72)
print("H : Y^2 = F0*F4*F8,  (m,n) = (%d,%d),  s = %s" % (m, n, s))
print("F0 =", F0)
print("F4 =", F4)
print("F8 =", F8)
print("f  =", f)
print("b = s*m^2*n^2 =", b, "   a =", a)
# сверка с формой Y^2 = b(t^6 + a t^4 + a t^2 + 1)
check = b*(T^6 + a*T^4 + a*T^2 + 1)
print("[sanity] f == b(t^6+a t^4+a t^2+1) :", f == check)
# сверка подстановкой конкретных t
for t0 in [0, 1, 2, -3, QQ(3)/7]:
    lhs = f(t0)
    rhs = (m^2+n^2*t0^2)*(s*(1+t0^2))*(n^2+m^2*t0^2)
    print("[sanity] t=%-6s f(t)=%-22s совпадает с F0F4F8: %s" % (t0, lhs, lhs == rhs))
print("[sanity] f палиндромичен (t -> 1/t) :", list(f) == list(reversed(list(f))))
print()

# ----------------------------------------------------------------------
# 1. ВЕЩЕСТВЕННОЕ МЕСТО
# ----------------------------------------------------------------------
print("-"*72)
print("1. ВЕЩЕСТВЕННОЕ МЕСТО R")
for name, Q in [("F0", F0), ("F4", F4), ("F8", F8)]:
    c2 = Q.coefficients(sparse=False)[2]
    c0 = Q.coefficients(sparse=False)[0]
    print("   %s: старший коэф. %s > 0, дискриминант %s < 0  =>  %s(t) > 0 для всех t in R"
          % (name, c2, Q.discriminant(), name))
print("   => f(t) > 0 для всех вещественных t; например f(0) =", f(0), ", f(1) =", f(1))
print("   Точки H(R) существуют (Y = ±sqrt(f(t))).   ПРЕПЯТСТВИЯ НЕТ.")
print()

# ----------------------------------------------------------------------
# 2. ПЛОХИЕ ПРОСТЫЕ
# ----------------------------------------------------------------------
print("-"*72)
print("2. ДИСКРИМИНАНТ И ПЛОХИЕ ПРОСТЫЕ")
disc_f = f.discriminant()
lc = f.leading_coefficient()
print("   disc(f)  =", factor(disc_f))
print("   lc(f)    =", factor(lc), " =", lc)
print("   f(0)     =", factor(f(0)))
# дискриминант бинарной секстики (учитывает и бесконечность):
# для однородной формы F6(u,v) плохие простые = простые | disc_f * lc
bad = sorted(set([p for p, _ in factor(disc_f)] + [p for p, _ in factor(lc)] + [2]))
print("   плохие простые (делители disc(f)*lc(f), плюс 2):", bad)
try:
    g2 = genus2reduction(0, R(f))
    print("   genus2reduction: кондуктор =", g2.conductor(), " =", factor(g2.conductor()))
    print("   genus2reduction: плохие простые =", sorted(g2.local_data.keys()))
except Exception as e:
    print("   genus2reduction недоступен:", e)
print()

# ----------------------------------------------------------------------
# Полный (исчерпывающий) перебор t in Z_p: существует ли t с f(t) in (Q_p)^2
# ----------------------------------------------------------------------
# Рекурсия по дискам. На диске t = t0 + p^k Z_p подставляем t = t0 + p^k*s,
# получаем G(s) in Z[s]. Пусть e = min валуация коэффициентов G, G = p^e * G1.
#   - для p нечётного разбиваем s mod p:
#       * G1(s0) -- единица: тогда v_p(G(s)) = e постоянна на подшайбе и
#         G(s) = p^e * (единица, сравнимая с G1(s0) mod p). Квадрат  <=>
#         e чётно и G1(s0) -- КВ мод p.  РЕШЕНИЕ ЕСТЬ / НЕТ -- решается точно.
#       * G1(s0) = 0 mod p: рекурсия в подшайбу s0 + p Z_p.
#   - для p = 2 разбиваем s mod 8: G1(s0+8s) = G1(s0) mod 8, так что для нечётного
#     G1(s0) класс квадрата определён точно (e чётно и G1(s0) = 1 mod 8).
# Возвращает True (точка найдена, с сертификатом), False (доказано, что точек нет),
# None (не хватило глубины -- НЕОПРЕДЕЛЁННО).
def exists_t_with_square(G, p, maxdepth=40, depth=0, shift=QQ(0), scale=QQ(1), cert=None):
    S = G.parent().gen()
    if G.is_zero():
        if cert is not None: cert.append(("f тождественно 0", shift))
        return True
    cs = [c for c in G.coefficients(sparse=False) if c != 0]
    e = min(c.valuation(p) for c in cs)
    G1 = G // p^e
    step = 8 if p == 2 else p
    branch = range(8) if p == 2 else range(p)
    undecided = []
    for s0 in branch:
        v = ZZ(G1(s0))
        if v % p != 0:
            if e % 2 == 0:
                ok = ((v % 8) == 1) if p == 2 else (kronecker(v, p) == 1)
                if ok:
                    tval = shift + scale*s0
                    if cert is not None:
                        cert.append(("t = %s,  f(t) = %s = p^%d * (квадрат-единица)" % (tval, G(s0), e), tval))
                    return True
            # иначе: на всей подшайбе v_p = e и класс единицы фиксирован -> квадратов нет
        else:
            undecided.append(s0)
    if not undecided:
        return False
    if depth >= maxdepth:
        return None
    res_unknown = False
    for s0 in undecided:
        H = G(s0 + step*S)
        r = exists_t_with_square(H, p, maxdepth, depth+1,
                                 shift + scale*s0, scale*step, cert)
        if r is True:
            return True
        if r is None:
            res_unknown = True
    return None if res_unknown else False

# ----------------------------------------------------------------------
# H(Q_p) != 0  <=>  exists (u:v) in P^1(Q_p) primitive with F6(u,v) квадрат
#   v -- единица  ->  t=u/v in Z_p, F6 = v^6 f(t):  f(t) квадрат
#   v = 0 mod p   ->  t'=v/u in pZ_p, F6 = u^6 g(t'), g = обращённый f
# f палиндромичен => g = f, но проверяем и ветку g отдельно (в т.ч. t'=0 = бесконечность).
# ----------------------------------------------------------------------
SS.<S> = ZZ[]
fS = SS(list(f))
gS = SS(list(reversed(list(f))))     # обращённый многочлен

# ----------------------------------------------------------------------
# Точный тест "x -- квадрат в Q_p" для рационального x
# ----------------------------------------------------------------------
def is_square_in_Qp(x, p):
    """Точно: x in (Q_p^*)^2 ? (x=0 считаем квадратом)"""
    x = QQ(x)
    if x == 0:
        return True
    e = x.valuation(p)
    u = x / p^e
    aa = u.numerator(); bb = u.denominator()
    if e % 2 != 0:
        return False
    if p == 2:
        return (aa*bb) % 8 == 1
    return kronecker(aa*bb, p) == 1

# самопроверка теста
assert is_square_in_Qp(9, 5) and is_square_in_Qp(QQ(1)/4, 7)
assert is_square_in_Qp(4, 2) and (not is_square_in_Qp(2, 2)) and (not is_square_in_Qp(5, 2))
assert is_square_in_Qp(17, 2) and (not is_square_in_Qp(3, 3)) and is_square_in_Qp(9*7, 3) == (kronecker(7,3)==1)

print("-"*72)
print("3. ТОЧКИ НА БЕСКОНЕЧНОСТИ")
print("   старший коэффициент lc =", lc, "=", factor(lc))
sq, rest = lc.squarefree_part(), None
print("   бесквадратная часть lc =", lc.squarefree_part())
print("   Две точки на бесконечности гладкой модели рациональны над K  <=>  lc in K^{*2}")
print("   над Q: lc =", lc, "не квадрат (%s = %s * %d^2) => точек на бесконечности над Q нет"
      % (lc, lc.squarefree_part(), sqrt(lc/lc.squarefree_part())))
print("   (это согласуется с наблюдением Codex: при t=0 и t=infinity рациональных точек нет,")
print("    т.к. f(0) = lc = 193*95^2 и класс 193 нетривиален)")
print()

print("   Над Q_v две точки на бесконечности рациональны <=> 193 in (Q_v^*)^2:")
row = []
for p in [2,3,5,7,11,13,17,19,193]:
    row.append("%s:%s" % (p, "да" if is_square_in_Qp(193, p) else "нет"))
print("      R : да (193>0);  " + ";  ".join(row))  # noqa
print("   Там, где 'нет', локальная точка обеспечивается конечным t (см. ниже).")
print()

print("-"*72)
print("4. ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ ПО ПРОСТЫМ")
# какие простые надо проверять руками:
#  - все плохие;
#  - все хорошие p < 17 (при хорошей редукции и p>=17 граница Вейля p+1-4sqrt(p)>0
#    даёт гладкую F_p-точку, которая поднимается по Гензелю).
to_check = sorted(set(bad) | set(prime_range(17)))
print("   Проверяем явно:", to_check)
print("   Для остальных p (хорошая редукция, p >= 17): граница Вейля")
print("      #H(F_p) >= p + 1 - 4*sqrt(p) > 0 при p >= 17,")
print("   и гладкая точка поднимается по лемме Гензеля  =>  H(Q_p) != 0.")
print()

results = {}
for p in to_check:
    cert = []
    r1 = exists_t_with_square(fS, p, maxdepth=40, cert=cert)
    if r1 is True:
        results[p] = (True, cert[-1][0] if cert else "")
    else:
        cert2 = []
        # ветка v = 0 mod p: t' in p Z_p, нужен g(t') квадрат
        gp = gS(p*S)
        r2 = exists_t_with_square(gp, p, maxdepth=40, cert=cert2)
        if r2 is True:
            results[p] = (True, "на карте в бесконечности: " + (cert2[-1][0] if cert2 else ""))
        elif r1 is False and r2 is False:
            results[p] = (False, "исчерпывающий перебор: точек нет")
        else:
            results[p] = (None, "глубина рекурсии исчерпана (r1=%s, r2=%s)" % (r1, r2))

for p in to_check:
    ok, why = results[p]
    tag = {True: "ЕСТЬ ТОЧКА", False: "ТОЧЕК НЕТ (доказано)", None: "НЕОПРЕДЕЛЕНО"}[ok]
    print("   p = %-5s : %-24s  %s" % (p, tag, why))
print()

# ----------------------------------------------------------------------
# Независимая перепроверка: явный перебор t = u/v с малыми u,v и прямой тест
# ----------------------------------------------------------------------
print("-"*72)
print("5. НЕЗАВИСИМАЯ ПЕРЕПРОВЕРКА (грубый поиск явного t для каждого p)")
for p in to_check:
    found = None
    # t in Z, |t| <= 200
    for t0 in range(0, 400):
        for sgn in ([1] if t0 == 0 else [1, -1]):
            x = f(sgn*t0)
            if is_square_in_Qp(x, p):
                found = ("t = %d" % (sgn*t0), x)
                break
        if found: break
    if found is None:
        # t = u/v
        for v0 in range(2, 60):
            for u0 in range(-200, 201):
                if gcd(u0, v0) != 1: continue
                x = f(QQ(u0)/v0)
                if is_square_in_Qp(x, p):
                    found = ("t = %d/%d" % (u0, v0), x)
                    break
            if found: break
    if found is None and is_square_in_Qp(lc, p):
        found = ("t = infinity (lc квадрат в Q_%d)" % p, lc)
    print("   p = %-5s : %s" % (p, ("найдено " + found[0]) if found else "явного t не найдено в диапазоне"))
print()

# ----------------------------------------------------------------------
# Проверка хороших простых 17 <= p < 200 прямым счётом точек (страховка)
# ----------------------------------------------------------------------
print("-"*72)
print("6. СТРАХОВКА: прямой счёт гладких F_p-точек для хороших p < 500")
worst = []
for p in prime_range(17, 500):
    if p in bad: continue
    Fp = GF(p)
    cnt = 0
    for t0 in Fp:
        val = Fp(f(ZZ(t0)))
        if val != 0 and val.is_square():
            cnt += 2
        elif val == 0:
            cnt += 1
    # бесконечность
    if Fp(lc).is_square() and Fp(lc) != 0:
        cnt += 2
    if cnt == 0:
        worst.append(p)
print("   простые p (17<=p<500, хорошая редукция) без гладких F_p-точек:", worst if worst else "нет")
print()

# ----------------------------------------------------------------------
print("="*72)
allgood = all(v[0] is True for v in results.values())
anybad  = any(v[0] is False for v in results.values())
unk     = [p for p, v in results.items() if v[0] is None]
print("ИТОГ:")
print("  вещественное место : точки есть")
print("  все проверенные p  :", "точки есть везде" if allgood else "СМ. ВЫШЕ")
if anybad:
    print("  НАЙДЕНО ПРЕПЯТСТВИЕ в местах:", [p for p, v in results.items() if v[0] is False])
if unk:
    print("  НЕОПРЕДЕЛЕНО в местах:", unk)
print("="*72)
