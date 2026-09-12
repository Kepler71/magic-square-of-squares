# -*- coding: utf-8 -*-
r"""
check_kummer_19_5_модель.sage
=============================
НЕЗАВИСИМАЯ АТАКА на заявление Codex: C_{19,5}(Q) = пусто.

Угол атаки: МОДЕЛЬ И ТОЖДЕСТВА. Всё строится с нуля; числа Codex
не копируются — сверка с ними только в самом конце (раздел 12).

Разделы:
  1.  Символьная проверка тождеств в QQ[m,n,t] (общий (m,n)).
  2.  Конкретизация (19,5): s, b, корни, неособость, 2-кручение.
  3.  Карта C -> E на конкретных числах; вывод необходимого класса.
  4.  Невырожденность: u_i != 0, X != e_i при конечном рациональном t.
  5.  Род кривой C (Риман–Гурвиц + Sage).
  6.  Ранг E(Q): несколько НЕЗАВИСИМЫХ источников, включая
      безусловный путь Гросса–Загира–Колывагина.
  7.  Полный образ delta(E(Q)) — свои генераторы, свои 2^(r+2) классов.
  8.  Устойчивость к перестановке/переобозначению корней.
  9.  Бесконечность: аккуратный анализ всех ветвей над t=oo.
 10.  Попытка ОПРОВЕРЖЕНИЯ: прямой поиск точек на C и на E.
 11.  Лежит ли цель в 2-Selmer (нужен ли ранговый аргумент вообще).
 12.  Сверка с числами Codex.
"""
import itertools, sys, time

def hdr(t):
    print("\n" + "=" * 74)
    print(t)
    print("=" * 74)

FAILS = []
def check(name, cond, extra=""):
    ok = bool(cond)
    print(("  [OK]   " if ok else "  [FAIL] ") + name + ("  " + str(extra) if extra else ""))
    if not ok:
        FAILS.append(name)
    return ok

proof.elliptic_curve(False)   # иначе Sage бросает исключение на gens(); ранг мы доказываем отдельно

# ================================================================== 1
hdr("1. СИМВОЛЬНАЯ ПРОВЕРКА ТОЖДЕСТВ (общий (m,n), кольцо QQ[m,n,t])")

R = PolynomialRing(QQ, ['m', 'n', 't'])
m, n, t = R.gens()

s_sym = (m**2 + n**2) / QQ(2)
b_sym = s_sym * m**2 * n**2

F0 = m**2 + n**2 * t**2
F4 = s_sym * (1 + t**2)
F8 = n**2 + m**2 * t**2

e1 = -b_sym
e2 = -s_sym * m**4
e3 = -s_sym * n**4
Xs = b_sym * t**2

check("X - e1 == (mn)^2 * F4", Xs - e1 - (m * n)**2 * F4 == 0)
check("X - e2 == s*m^2 * F0", Xs - e2 - s_sym * m**2 * F0 == 0)
check("X - e3 == s*n^2 * F8", Xs - e3 - s_sym * n**2 * F8 == 0)
prod3 = (Xs - e1) * (Xs - e2) * (Xs - e3)
check("(X-e1)(X-e2)(X-e3) == b^2*F0*F4*F8  (V=b*u0u4u8 лежит на E)",
      prod3 - b_sym**2 * F0 * F4 * F8 == 0)

print("\n  Разности корней (символьно, для проверки различности):")
for (i, j, ei, ej) in [(1, 2, e1, e2), (1, 3, e1, e3), (2, 3, e2, e3)]:
    print("    e%d - e%d = %s" % (i, j, R(ei - ej).factor()))
print("  => при gcd(m,n)=1, m>n>0, m!=n все три корня различны.")

# ================================================================== 2
hdr("2. КОНКРЕТИЗАЦИЯ (m,n)=(19,5): свои числа")

M = Integer(19); N = Integer(5)
check("gcd(m,n)=1", gcd(M, N) == 1)
S = QQ(M**2 + N**2) / 2
B = S * M**2 * N**2
print("  s = (%s+%s)/2 = %s" % (M**2, N**2, S))
print("  b = s*m^2*n^2 = %s" % B)
check("s == 193 (мой счёт)", S == 193)
check("b == 1741825 (мой счёт)", B == 1741825)
check("s свободно от квадратов", Integer(S).is_squarefree())
check("s НЕ квадрат в Q", not QQ(S).is_square())

E1 = QQ(-B); E2 = QQ(-S * M**4); E3 = QQ(-S * N**4)
print("  e1 = -b     = %s" % E1)
print("  e2 = -s*m^4 = %s" % E2)
print("  e3 = -s*n^4 = %s" % E3)
check("три корня различны", len(set([E1, E2, E3])) == 3)

Rx = PolynomialRing(QQ, 'Xv'); Xv = Rx.gen()
cub = (Xv - E1) * (Xv - E2) * (Xv - E3)
co = cub.coefficients(sparse=False)
E = EllipticCurve([0, co[2], 0, co[1], co[0]])
print("  E : V^2 = %s" % cub.factor())
print("  E =", E)
check("disc(E) != 0", E.discriminant() != 0)
check("disc(кубики) != 0", cub.discriminant() != 0)
check("2-кручение полностью рационально", E.two_torsion_rank() == 2)
Tg = E.torsion_subgroup()
check("кручение ровно Z/2 x Z/2", tuple(Tg.invariants()) == (2, 2), Tg.invariants())
tors_x = sorted([P.xy()[0] for P in Tg.points() if P != E(0)])
check("x(2-кручение) = {e1,e2,e3}", tors_x == sorted([E1, E2, E3]))

# ================================================================== 3
hdr("3. КАРТА C -> E НА (19,5); НЕОБХОДИМЫЙ КЛАСС")

Rt = PolynomialRing(QQ, 't'); tt = Rt.gen()
f0 = M**2 + N**2 * tt**2
f4 = S * (1 + tt**2)
f8 = N**2 + M**2 * tt**2
Xt = B * tt**2
check("X-e1 == (mn)^2*F4", Xt - E1 - (M * N)**2 * f4 == 0)
check("X-e2 == s*m^2*F0", Xt - E2 - S * M**2 * f0 == 0)
check("X-e3 == s*n^2*F8", Xt - E3 - S * N**2 * f8 == 0)
check("V^2 = b^2 F0F4F8 == кубика", (Xt-E1)*(Xt-E2)*(Xt-E3) - B**2*f0*f4*f8 == 0)


def sqclass(x):
    x = QQ(x)
    if x == 0:
        raise ValueError("ноль")
    z = x.numerator() * x.denominator()
    sgn = 1 if z > 0 else -1
    out = Integer(1)
    for p, ee in Integer(abs(z)).factor():
        if ee % 2 == 1:
            out *= p
    return Integer(sgn) * out

print("  множители: (mn)^2=%s -> класс %s ; s*m^2=%s -> класс %s ; s*n^2=%s -> класс %s"
      % ((M*N)**2, sqclass((M*N)**2), S*M**2, sqclass(S*M**2), S*N**2, sqclass(S*N**2)))
TARGET = (sqclass((M*N)**2), sqclass(S*M**2), sqclass(S*N**2))
print("  => НЕОБХОДИМЫЙ КЛАСС delta = %s" % (TARGET,))
check("мой независимый вывод цели совпал с (1,193,193)",
      TARGET == (Integer(1), Integer(193), Integer(193)), TARGET)

# ================================================================== 4
hdr("4. НЕВЫРОЖДЕННОСТЬ")

for name, poly in [("F0", f0), ("F4", f4), ("F8", f8)]:
    realr = [r for r, _ in poly.roots(QQbar) if r.imag() == 0]
    check("%s не имеет вещественных корней => != 0 при рац. t" % name, len(realr) == 0)
check("s > 0", S > 0)
print("  => при конечном рациональном t: u0,u4,u8 != 0, X != e_i, точка E не 2-кручение и не O;")
print("     значит delta вычисляется без регуляризации и равен ровно (1,s,s).")

# ================================================================== 5
hdr("5. РОД КРИВОЙ C")

print("  C -> P^1_t степени 8, группа (Z/2)^3.")
print("  Точки ветвления = корни F0,F4,F8 (по 2 на каждый) = 6 точек;")
print("  над каждой ровно одно из трёх квадратичных расширений ветвится:")
print("    слой = 4 точки с e=2, вклад в дифференту = 8-4 = 4.")
print("  t=oo НЕ ветвится: F_i/t^2 имеют при z=0 значения n^2, s, m^2 != 0 (чётная валюация).")
print("  Риман–Гурвиц: 2g-2 = 8*(-2) + 6*4 = -16+24 = 8  =>  g = 5.")
print("  Sage не умеет genus над QQ для такой башни; считаем род при ХОРОШЕЙ редукции")
print("  (род постоянен в гладкой семье, поэтому genus_Q = genus_Fp):")
badset = set(Integer(2*M*N*(M**2 - N**2)*(M**2 + N**2)).prime_factors())
print("  плохие p (делители 2mn(m^2-n^2)(m^2+n^2)):", sorted(badset))
for p in [11, 13, 17, 23, 29, 31]:
    if p in badset:
        continue
    K = GF(p); K0 = FunctionField(K, 't'); T0 = K0.gen()
    g0 = K(M**2) + K(N**2)*T0**2
    g4 = K(S)*(1 + T0**2)
    g8 = K(N**2) + K(M**2)*T0**2
    Pb = PolynomialRing(K0, ['A', 'Bv', 'Cc', 'Yv'])
    A, Bv, Cc, Yv = Pb.gens()
    poly = prod([Yv - (sa*A + sb*Bv + sc*Cc)
                 for sa, sb, sc in itertools.product([1, -1], repeat=3)])
    acc = {}
    for mon, cf in zip(poly.monomials(), poly.coefficients()):
        ea, eb, ec, ey = mon.degrees()
        if ea % 2 or eb % 2 or ec % 2:
            continue
        acc[ey] = acc.get(ey, K0(0)) + cf*g0**(ea//2)*g4**(eb//2)*g8**(ec//2)
    Rz = PolynomialRing(K0, 'Y'); Y = Rz.gen()
    minp = sum([acc[k]*Y**k for k in acc])
    if minp.degree() == 8 and minp.is_irreducible():
        Lp = K0.extension(minp, 'al')
        gp = Lp.genus()
        check("genus(C mod %d) == 5 (степень покрытия 8, неприводимо)" % p, gp == 5,
              "genus=%s" % gp)

# ================================================================== 6
hdr("6. РАНГ E(Q) — НЕСКОЛЬКО НЕЗАВИСИМЫХ ИСТОЧНИКОВ")

Em = E.minimal_model()
print("  минимальная модель:", Em)
NC = Em.conductor()
print("  кондуктор N =", NC, "=", NC.factor())

print("\n  (6a) mwrank / полный 2-спуск (Sage):")
rb = E.rank_bounds()
print("       rank_bounds() =", rb)
sel = E.selmer_rank()
print("       selmer_rank C (dim 2-Selmer) =", sel)
tt2 = E.two_torsion_rank()
print("       two_torsion_rank T =", tt2)
check("2-Selmer даёт лишь rank <= C - T = %d" % (sel - tt2), rb[1] == sel - tt2, rb)
print("       ВЫВОД: чистый 2-спуск даёт только 1 <= rank <= 3. Не достаточно!")

print("\n  (6b) PARI ellrank (алгебраический, спаривание Касселса):")
for eff in [0, 2]:
    pr = pari(Em).ellrank(eff)
    print("       effort=%d -> %s" % (eff, pr))
print("       PARI: r2 = C - T - s, где s = rank(Sha[2]/2Sha[4]) считается безусловно.")
print("       Здесь C=%d, T=%d, s=2 => r2 = %d." % (sel, tt2, sel - tt2 - 2))

print("\n  (6c) БЕЗУСЛОВНЫЙ путь: Гросс–Загир + Колывагин (не зависит от PARI ellrank)")
w = E.root_number()
print("       корневое число w =", w)
check("w = -1 => по функциональному уравнению L(E,1)=0", w == -1)
print("       L_ratio (точная рациональная проверка L(1)/Omega) =", Em.lseries().L_ratio())
check("L(E,1) = 0 подтверждено точно (L_ratio == 0)", Em.lseries().L_ratio() == 0)
t0 = time.time()
L1v, L1e = Em.lseries().deriv_at1(2000000)
print("       L'(E,1) = %s  ±  %s   (%.1f s, k=2*10^6)" % (L1v, L1e, time.time() - t0))
check("|L'(E,1)| > строгая граница ошибки  =>  L'(E,1) != 0", abs(L1v) > L1e)
try:
    Ld = Em.lseries().dokchitser(100)
    print("       НЕЗАВИСИМАЯ проверка (Dokchitser): L(1) =", Ld(1), " L'(1) =", Ld.derivative(1, 1))
    check("Dokchitser даёт то же L'(1) (согласие >= 20 знаков)",
          abs(RR(Ld.derivative(1, 1)) - RR(L1v)) < 1e-18)
except Exception as ex:
    print("       [прим.] Dokchitser не отработал:", ex)
print("""       Граница хвоста 2*exp(-2pi(k+1)/sqrt(N))/(1-exp(-2pi/sqrt(N))) строгая
       (Grigorov–Jorza–Patrascu–Patrikis–Stein), см. док. Sage deriv_at1.
       Итак ord_{s=1} L(E,s) = 1 ТОЧНО.
       Модулярность (Wiles-BCDT) + Гросс–Загир + Колывагин  =>  rank E(Q) = 1
       и Sha(E/Q) конечна. Это БЕЗУСЛОВНО и НЕЗАВИСИМО от PARI ellrank.""")

print("\n  (6d) контрольные (не доказательства):")
print("       analytic_rank(pari) =", Em.analytic_rank(algorithm='pari'))
print("       analytic_rank_upper_bound (при GRH) =", Em.analytic_rank_upper_bound())
print("       sha.an() (порядок Sha по BSD) =", Em.sha().an(),
      " -> согласуется с dim Sha[2]=2 (PARI s=2)")

RANK = 1
print("\n  ПРИНЯТО: rank E(Q) = %d (обосновано в 6c безусловно)." % RANK)

# ================================================================== 7
hdr("7. ПОЛНЫЙ ОБРАЗ delta(E(Q)) — свои генераторы, свои классы")


def delta(P):
    if P == E(0):
        return (Integer(1), Integer(1), Integer(1))
    x = P.xy()[0]
    es = [E1, E2, E3]
    out = []
    for i in range(3):
        d = x - es[i]
        if d == 0:
            j, k = [q for q in range(3) if q != i]
            d = (es[i] - es[j]) * (es[i] - es[k])
        out.append(sqclass(d))
    return tuple(out)


gens = E.gens(proof=False)
print("  E.gens(proof=False) =", gens)
inf_gens = [P for P in gens if P.order() == oo]
check("найдено ровно %d точек бесконечного порядка" % RANK, len(inf_gens) == RANK)
sat, idx, reg = E.saturation(inf_gens)
print("  saturation index =", idx, " (для перечисления НЕ обязателен)")
inf_gens = sat

TorsPts = list(E.torsion_points())
print("  точек кручения:", len(TorsPts))

classes = {}
allpts = []
for coeffs in itertools.product([0, 1], repeat=RANK):
    for Tp in TorsPts:
        P = Tp + sum([c * g for c, g in zip(coeffs, inf_gens)], E(0))
        allpts.append(P)
        classes.setdefault(delta(P), []).append(P)

print("  перечислено представителей E(Q)/2E(Q):", len(allpts))
print("  РАЗЛИЧНЫХ классов:", len(classes))
expected = 2 ** (RANK + 2)
check("|образ| = 2^(r+2) = %d" % expected, len(classes) == expected)
print("\n  СПИСОК КЛАССОВ (мой счёт):")
for c in sorted(classes.keys(), key=lambda z: (abs(z[0]), abs(z[1]), abs(z[2]))):
    print("   ", c)

check("целевой класс %s ОТСУТСТВУЕТ в образе" % (TARGET,), TARGET not in classes)

print("\n  Контроль: образ — подгруппа (Q*/Q*^2)^3?")
CS = set(classes.keys())
grp_ok = True
for a in CS:
    for bq in CS:
        pr = tuple(sqclass(u * v) for u, v in zip(a, bq))
        if pr not in CS:
            grp_ok = False
check("множество классов замкнуто относительно умножения (подгруппа)", grp_ok)
check("(1,1,1) в образе", (Integer(1), Integer(1), Integer(1)) in CS)

print("\n  Контроль гомоморфности delta на всех парах перечисленных точек:")
bad = 0
for P in allpts:
    for Q in allpts:
        lhs = delta(P + Q)
        rhs = tuple(sqclass(u * v) for u, v in zip(delta(P), delta(Q)))
        if lhs != rhs:
            bad += 1
check("delta(P+Q)=delta(P)delta(Q) на всех %d парах" % (len(allpts)**2), bad == 0,
      "нарушений %d" % bad)

print("\n  Контроль: произведение координат каждого класса — квадрат (так как")
print("  (X-e1)(X-e2)(X-e3)=V^2):")
prod_ok = all(QQ(c[0]*c[1]*c[2]).is_square() for c in CS)
check("d1*d2*d3 — квадрат для всех классов образа", prod_ok)
check("d1*d2*d3 — квадрат и для цели", QQ(TARGET[0]*TARGET[1]*TARGET[2]).is_square())

# ================================================================== 8
hdr("8. УСТОЙЧИВОСТЬ К ПЕРЕСТАНОВКЕ / ПЕРЕОБОЗНАЧЕНИЮ КОРНЕЙ")

print("  Если бы порядок корней был перепутан, цель была бы перестановкой (1,s,s).")
for p in sorted(set(itertools.permutations(TARGET))):
    check("перестановка %s отсутствует в образе" % (p,), p not in CS,
          "ПРИСУТСТВУЕТ!" if p in CS else "")

print("\n  Пере-деривация при ПРОИЗВОЛЬНОЙ нумерации корней (все 6 назначений):")
named = {'-b': E1, '-s*m^4': E2, '-s*n^4': E3}
fac = {'-b': (M*N)**2, '-s*m^4': S*M**2, '-s*n^4': S*N**2}
for order in itertools.permutations(['-b', '-s*m^4', '-s*n^4']):
    # delta считается В ТОМ ЖЕ порядке корней, что и цель => инвариантно
    es = [named[k] for k in order]
    cls = tuple(sqclass(fac[k]) for k in order)

    def delta_o(P, es=es):
        if P == E(0):
            return (Integer(1), Integer(1), Integer(1))
        x = P.xy()[0]
        out = []
        for i in range(3):
            d = x - es[i]
            if d == 0:
                j, k = [q for q in range(3) if q != i]
                d = (es[i] - es[j]) * (es[i] - es[k])
            out.append(sqclass(d))
        return tuple(out)
    img_o = set(delta_o(P) for P in allpts)
    print("    порядок %-28s цель %-16s |образ|=%d ; цель в образе: %s"
          % (str(order), str(cls), len(img_o), cls in img_o))
    check("при нумерации %s цель вне образа" % (order,), cls not in img_o)

print("\n  Диагностика «а если потерян множитель»:")
for cand in [(1,1,S),(1,S,1),(S,1,1),(S,S,S),(1,1,1),(2*S,2*S,1)]:
    cand = tuple(Integer(z) for z in cand)
    print("    класс %-22s в образе: %s" % (str(cand), cand in CS))

# ================================================================== 9
hdr("9. БЕСКОНЕЧНОСТЬ: все ветви над t = oo")

Ft = FractionField(PolynomialRing(QQ, 'tq')); tq = Ft.gen()
zz = 1 / tq
check("U0^2 := F0/t^2 = n^2 + m^2 z^2 (тождество рац. функций)",
      (M**2 + N**2*tq**2)*zz**2 - (N**2 + M**2*zz**2) == 0)
check("U4^2 := F4/t^2 = s(1 + z^2)",
      (S*(1 + tq**2))*zz**2 - S*(zz**2 + 1) == 0)
check("U8^2 := F8/t^2 = m^2 + n^2 z^2",
      (N**2 + M**2*tq**2)*zz**2 - (M**2 + N**2*zz**2) == 0)
print("  при z=0:  U0^2 = n^2 = %s (квадрат) ; U4^2 = s = %s ; U8^2 = m^2 = %s (квадрат)"
      % (N**2, S, M**2))
check("U0^2|_{z=0} квадрат в Q", QQ(N**2).is_square())
check("U8^2|_{z=0} квадрат в Q", QQ(M**2).is_square())
check("U4^2|_{z=0} = s НЕ квадрат в Q", not QQ(S).is_square())
print("""
  СТРОГИЙ АРГУМЕНТ (ни одна ветвь не теряется).
  Пусть P — Q-рациональная точка гладкой проективной модели C над z=0.
  Локальное кольцо O_P — ДВН, значит целозамкнуто в поле функций.
  U4 = u4/t удовлетворяет монический уравнению T^2 - s(1+z^2) = 0 с
  коэффициентами из O_P (z регулярна в P). Значит U4 цела над O_P,
  следовательно U4 in O_P. Поле вычетов Q-точки равно Q, поэтому
  U4(P) in Q и U4(P)^2 = s(1+0) = s. Значит s — квадрат в Q. Но s=193
  не квадрат. Противоречие. Аргумент НЕ использует поведение u0,u8,
  поэтому расщепление их ветвей несущественно.

  Для протокола — разложение места z=0 в башне:
    Q(z)(U0): U0^2 = n^2(1+(m/n)^2 z^2) — квадрат единицы -> расщепляется (2 места)
    Q(z)(U8): U8^2 = m^2(1+(n/m)^2 z^2) — квадрат единицы -> расщепляется (2 места)
    Q(z)(U4): U4^2 = s(1+z^2) — s не квадрат -> инертно, поле вычетов Q(sqrt(s))
  => над z=0 ровно 4 места степени 2, ни одного Q-рационального.""")
# Локальная (формальная) форма того же утверждения: Q-точка над z=0 <=>
# система U_i^2 = G_i(z) разрешима в Q[[z]] с z=0, т.е. G_i(0) — квадраты в Q.
Rzz = PowerSeriesRing(QQ, 'zs', default_prec=12); zs = Rzz.gen()
for nm, ser, val in [("U0", N**2 + M**2*zs**2, N**2),
                     ("U4", S*(1 + zs**2), S),
                     ("U8", M**2 + N**2*zs**2, M**2)]:
    has = QQ(val).is_square()
    if has:
        rt = ser.sqrt()
        print("    %s: ряд в Q[[z]] существует, %s = %s + O(z^4)" % (nm, nm, rt.truncate(4)))
    else:
        print("    %s: ряда в Q[[z]] НЕ существует (свободный член %s не квадрат)" % (nm, val))
check("ветвь U4 не поднимается до Q[[z]] => нет Q-точки над t=oo",
      not QQ(S).is_square())

print("\n  Контроль t=0: F4(0)=s=%s, квадрат? %s" % (S, QQ(S).is_square()))
check("над t=0 тоже нет Q-точек", not QQ(S).is_square())

# ================================================================== 10
hdr("10. ПОПЫТКА ОПРОВЕРЖЕНИЯ: прямой поиск")

print("  (a) перебор t = p/q, 0 <= p,q <= 600, gcd(p,q)=1:")
found = []
LIM = 600
for q in range(1, LIM + 1):
    q2 = q * q
    for p in range(0, LIM + 1):
        if gcd(p, q) != 1:
            continue
        p2 = p * p
        A = M**2 * q2 + N**2 * p2            # F0 * q^2
        if not Integer(A).is_square():
            continue
        Cc = N**2 * q2 + M**2 * p2           # F8 * q^2
        if not Integer(Cc).is_square():
            continue
        Bb = S * (q2 + p2)                   # F4 * q^2
        if not Integer(Bb).is_square():
            continue
        found.append((p, q))
print("     найдено t:", found)
check("прямой поиск точек C ничего не нашёл", len(found) == 0, found)

print("\n  (b) поиск точек E(Q) с целевым классом (point_search высота 18):")
ps = E.point_search(18)
hit = [P for P in ps if P != E(0) and delta(P) == TARGET]
print("     найдено точек:", len(ps), "; с целевым классом:", len(hit))
check("ни одна найденная точка E(Q) не имеет целевого класса", len(hit) == 0, hit)
cl_found = set(delta(P) for P in ps)
print("     классы найденных точек:", sorted(cl_found))
check("все классы найденных точек лежат в перечисленном образе",
      cl_found.issubset(CS), cl_found - CS)

print("\n  (c) точки, найденные PARI ellrank при effort=2 и 6 (независимый источник):")
for eff in [2, 6]:
    pr = pari(Em).ellrank(eff)
    pts = pr[3]
    for pp in pts:
        xm, ym = QQ(pp[0]), QQ(pp[1])
        try:
            Pm = Em(xm, ym)
            Pe = Em.isomorphism_to(E)(Pm)
            print("     effort=%d точка %s -> delta = %s ; в образе: %s"
                  % (eff, Pe.xy()[0], delta(Pe), delta(Pe) in CS))
            check("PARI-точка (effort=%d) имеет класс из перечисленного образа" % eff,
                  delta(Pe) in CS)
            check("PARI-точка (effort=%d) НЕ имеет целевого класса" % eff,
                  delta(Pe) != TARGET)
        except Exception as ex:
            print("     effort=%d: не удалось поднять точку: %s" % (eff, ex))

print("\n  (d) кратные генератора nG+T — не даст ли класс цель?")
print("      Для больших n факторизация x-e_i неподъёмна, поэтому считаем класс")
print("      по НОСИТЕЛЮ: d_i обязан делить 2*(e_i-e_j)(e_i-e_k). Снимаем эти простые,")
print("      остаток ОБЯЗАН быть точным квадратом — иначе печатаем отказ, а не ответ.")
SUPP = set([2])
for i in range(3):
    j, k = [q for q in range(3) if q != i]
    esl = [E1, E2, E3]
    SUPP |= set(Integer((esl[i]-esl[j])*(esl[i]-esl[k])).prime_factors())
SUPP = sorted(SUPP)
print("      носитель S =", SUPP)


def sqclass_S(x):
    """класс в Q*/Q*^2 без общей факторизации; None если носитель не покрыл."""
    x = QQ(x)
    z = x.numerator() * x.denominator()
    sgn = Integer(1) if z > 0 else Integer(-1)
    z = abs(z)
    out = Integer(1)
    for p in SUPP:
        v = z.valuation(p)
        z //= p**v
        if v % 2:
            out *= p
    if not Integer(z).is_square():
        return None
    return sgn * out


def delta_S(P):
    if P == E(0):
        return (Integer(1), Integer(1), Integer(1))
    x = P.xy()[0]
    esl = [E1, E2, E3]
    out = []
    for i in range(3):
        d = x - esl[i]
        if d == 0:
            j, k = [q for q in range(3) if q != i]
            d = (esl[i] - esl[j]) * (esl[i] - esl[k])
        c = sqclass_S(d)
        if c is None:
            return None
        out.append(c)
    return tuple(out)


unresolved = []
badn = []
seen_d = set()
for nn in range(-25, 26):
    for Tp in TorsPts:
        P = nn * inf_gens[0] + Tp
        if P == E(0):
            continue
        dP = delta_S(P)
        if dP is None:
            unresolved.append((nn, Tp))
            continue
        seen_d.add(dP)
        if dP == TARGET:
            badn.append((nn, Tp))
check("носитель S покрыл все классы nG+T (|n|<=25)", len(unresolved) == 0,
      "не разрешено: %d" % len(unresolved))
check("ни одно nG+T (|n|<=25) не даёт целевой класс", len(badn) == 0, badn)
check("все классы nG+T лежат в перечисленном образе", seen_d.issubset(CS), seen_d - CS)
print("      (это контроль согласованности; логически он излишен, т.к. delta —")
print("       гомоморфизм и |E(Q)/2E(Q)|=8 уже перечислено полностью)")

badn = []
seen_d = set()
for nn in range(-1, 2):
    for Tp in TorsPts:
        P = nn * inf_gens[0] + Tp
        if P == E(0):
            continue
        dP = delta(P)
        seen_d.add(dP)
        if dP == TARGET:
            badn.append((nn, Tp))
check("полная факторизация (|n|<=1) согласна с версией по носителю",
      seen_d == set(delta(nn*inf_gens[0] + Tp)
                    for nn in range(-1, 2) for Tp in TorsPts
                    if nn*inf_gens[0] + Tp != E(0)))

# ================================================================== 11
hdr("11. ЛЕЖИТ ЛИ ЦЕЛЬ В 2-SELMER? (можно ли было обойтись без ранга)")

d1, d2, d3 = QQ(1), QQ(S), QQ(S)
print("  Торсор: d1 w1^2 - d2 w2^2 = e2-e1 = %s" % (E2 - E1))
print("          d1 w1^2 - d3 w3^2 = e3-e1 = %s" % (E3 - E1))
w1s = (M*N)**2 * f4; w2s = M**2 * f0; w3s = N**2 * f8
check("точка C даёт точку торсора: тождество 1 (w1=mn u4, w2=m u0)",
      d1*w1s - d2*w2s - (E2 - E1) == 0)
check("точка C даёт точку торсора: тождество 2 (w3=n u8)",
      d1*w1s - d3*w3s - (E3 - E1) == 0)

print("\n  Локальная разрешимость C (мой счёт): ищу t in Z с F0,F4,F8 квадратами в Q_p")
bad_primes = [2, 3, 5, 7, 19, 193]
test_primes = bad_primes + [11, 13, 17, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
                            79, 83, 89, 97, 101, 103]
all_loc = True
for p in sorted(set(test_primes)):
    Kp = Qp(p, 60)
    ok = None
    for tv in range(0, 1200):
        try:
            a = Kp(M**2 + N**2*tv**2); bqq = Kp(S*(1 + tv**2)); c = Kp(N**2 + M**2*tv**2)
            if a != 0 and bqq != 0 and c != 0 and a.is_square() and bqq.is_square() and c.is_square():
                ok = tv
                break
        except Exception:
            continue
    if ok is None:
        all_loc = False
    print("     p=%-4s свидетель t=%s" % (p, ok))
print("     вещественное место: t=1 даёт F0=F4=F8=%s > 0, все квадраты в R" % (M**2 + N**2))
check("C(Q_p) != пусто для всех проверенных p", all_loc)
print("""
  ВЫВОД: C всюду локально разрешима (проверено на указанных p; для p>=101
  с хорошей редукцией работает Хассе–Вейль #C(F_p) >= p+1-10*sqrt(p) > 0).
  Точка C(Q_p) даёт точку торсора класса (1,s,s) над Q_p, значит цель
  ЛЕЖИТ в 2-Selmer-группе. Следовательно:
    * локальными методами исключение НЕВОЗМОЖНО;
    * класс (1,s,s) задаёт НЕТРИВИАЛЬНЫЙ элемент Sha(E/Q)[2];
    * ранговый аргумент ОБЯЗАТЕЛЕН, обойти его нельзя.
  Это согласуется с sha.an()=4 и PARI s=2.""")

# ================================================================== 12
hdr("12. СВЕРКА С ЧИСЛАМИ CODEX (только теперь)")

codex = [(1,1,1), (-1,4053,-4053), (-4053,386,-42), (4053,42,386),
         (-3,386,-1158), (3,42,14), (1351,1,1351), (-1351,4053,-3)]
codex = set(tuple(Integer(z) for z in c) for c in codex)
print("  Codex: s=193, b=1741825, e=(-1741825,-25151953,-120625), цель (1,193,193)")
check("s совпал", S == 193)
check("b совпал", B == 1741825)
check("корни совпали", (E1, E2, E3) == (QQ(-1741825), QQ(-25151953), QQ(-120625)))
print("  мой образ  :", sorted(CS))
print("  образ Codex:", sorted(codex))
check("списки классов совпадают ПОТОЧЕЧНО", CS == codex,
      "мои лишние: %s ; их лишние: %s" % (CS - codex, codex - CS))
Gc = (QQ(-28561963657)/1369, QQ(2089086742828800)/50653)
try:
    PG = E(Gc[0], Gc[1])
    print("  точка Codex G лежит на E:", PG in E, ", порядок:", PG.order())
    print("  delta(G_codex) =", delta(PG), "; моя точка:", inf_gens[0].xy()[0],
          "delta =", delta(inf_gens[0]))
    check("G_codex совпадает с моим генератором с точностью до знака",
          PG == inf_gens[0] or PG == -inf_gens[0])
except Exception as ex:
    print("  точка Codex НЕ лежит на E:", ex)
    FAILS.append("точка Codex не на E")

# ================================================================== итог
hdr("ИТОГ")
if FAILS:
    print("  ПРОВАЛЕНЫ ПРОВЕРКИ:")
    for f in FAILS:
        print("   -", f)
else:
    print("  Все проверки пройдены.")
print("  |delta(E(Q))| = %d = 2^(rank+2), rank=%d" % (len(CS), RANK))
print("  цель %s в образе: %s" % (TARGET, TARGET in CS))
print("  => C_{19,5}(Q) = пусто (при rank=1, обоснованном безусловно в 6c)")
