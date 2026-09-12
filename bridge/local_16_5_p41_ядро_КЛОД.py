# -*- coding: utf-8 -*-
# ЯДРО прямой проверки локального исключения (m,n) = (16,5) при p = 41.
# Чистый Python (никакого препарсера Sage), точная арифметика Fraction/isqrt.
# Запускается либо напрямую (python3), либо из обёртки local_16_5_p41_прямая.sage.
#
# Проверяемое утверждение (Codex, PARALLEL_RESULT_FULL_NINE_126_2026-09-12.md):
#   полная система девяти клеток G1 для (16,5) не имеет решений уже над F_41.

from fractions import Fraction as Fr
from math import isqrt

M, N, P = 16, 5, 41

def head(s):
    print()
    print("=" * 100)
    print(s)
    print("=" * 100)

# ==============================================================================================
def section0_sympy():
    """§0. Вывод девяти клеток из АКСИОМ магического квадрата, а не из чужих файлов."""
    head("§0. СИМВОЛИЧЕСКИЙ ВЫВОД ДЕВЯТИ КЛЕТОК ИЗ АКСИОМ (sympy)")
    import sympy as sp
    a, b, m, n = sp.symbols('a b m n')
    # Аксиомы ассоциативного магического квадрата: центр C, все линии = 3C, противоположные = 2C.
    # c0=C+al, c8=C-al, c1=C+be, c7=C-be, c2=C+ga, c6=C-ga, c3=C+de, c5=C-de,
    # строка 0 даёт al+be+ga=0, столбец 0 даёт de = ga-al. Два свободных параметра.
    # Семейство G1 = «четыре рёберные клетки суть квадраты линейных форм»:
    c1 = sp.expand((a*m + b*n)**2)
    c3 = sp.expand((a*m - b*n)**2)
    c5 = sp.expand((a*n + b*m)**2)
    c7 = sp.expand((a*n - b*m)**2)
    same_center = sp.expand((c1 + c7) - (c3 + c5)) == 0
    C = sp.expand((c1 + c7) / 2)
    be = sp.expand(c1 - C)
    de = sp.expand(c3 - C)
    ga = sp.expand((de - be) / 2)           # из al+be+ga=0 и de=ga-al следует 2ga = de-be
    al = sp.expand(-be - ga)
    c0, c8 = sp.expand(C + al), sp.expand(C - al)
    c2, c6 = sp.expand(C + ga), sp.expand(C - ga)
    c4 = C
    cells = [c0, c1, c2, c3, c4, c5, c6, c7, c8]
    LINES = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
    ok_lines = all(sp.expand(cells[i]+cells[j]+cells[k] - 3*c4) == 0 for i,j,k in LINES)
    ok_pairs = all(sp.expand(cells[i]+cells[8-i] - 2*c4) == 0 for i in range(4))

    print("  обе рёберные пары дают ОДИН центр (c1+c7 = c3+c5):", same_center)
    print("  c4 = C =", sp.factor(C), "  == (a^2+b^2)(m^2+n^2)/2 :",
          sp.expand(C - (a*a+b*b)*(m*m+n*n)/2) == 0)
    print("  ВЫВЕДЕНО из аксиом:")
    print("    c0 =", c0, "  == n^2 a^2 + m^2 b^2 :", sp.expand(c0 - (n*n*a*a + m*m*b*b)) == 0)
    print("    c2 =", c2, "  == C - 2mn*ab :", sp.expand(c2 - (C - 2*m*n*a*b)) == 0)
    print("    c6 =", c6, "  == C + 2mn*ab :", sp.expand(c6 - (C + 2*m*n*a*b)) == 0)
    print("    c8 =", c8, "  == m^2 a^2 + n^2 b^2 :", sp.expand(c8 - (m*m*a*a + n*n*b*b)) == 0)
    print("  все восемь линий = 3*c4:", ok_lines, "; все четыре противоположные пары = 2*c4:", ok_pairs)

    # какие клетки — квадраты ТОЖДЕСТВЕННО. Корректный тест: 2*c_i должен быть полным квадратом
    # в кольце Q[a,b,m,n]. (Множитель 2 у c2,c4,c6 не мешает: проверяем квадратность самого c_i.)
    def poly_square(e):
        e = sp.expand(e)
        if e == 0:
            return True
        f = sp.factor_list(sp.together(e))
        # e — полный квадрат <=> рациональный коэффициент квадрат и все кратности чётны
        coeff, facs = f
        if not (sp.nsimplify(coeff).is_rational and sp.sqrt(coeff).is_rational):
            return False
        return all(k % 2 == 0 for _, k in facs)
    auto = [i for i in range(9) if poly_square(cells[i])]
    nonauto = [i for i in range(9) if i not in auto]
    print("  ТОЖДЕСТВЕННЫЕ (автоматические) квадраты — клетки", auto, "; ожидались рёбра [1, 3, 5, 7]:",
          auto == [1, 3, 5, 7])
    print("  НЕАВТОМАТИЧЕСКИЕ клетки (их и надо проверять) —", nonauto,
          "; ожидались [0, 2, 4, 6, 8]:", nonauto == [0, 2, 4, 6, 8])
    print("  (контроль к чужому логу: c0 = a^2 n^2 + b^2 m^2 НЕ является тождественным квадратом:",
          not poly_square(c0), ")")
    return auto == [1, 3, 5, 7] and ok_lines and ok_pairs and same_center

# ==============================================================================================
def cells_Q(mm, nn, aa, bb):
    """девять клеток как Fraction для проективного параметра t = aa/bb (однородно, x bb^2)"""
    mm, nn, aa, bb = int(mm), int(nn), int(aa), int(bb)
    F4 = Fr((mm*mm + nn*nn) * (aa*aa + bb*bb), 2)
    return [Fr(nn*nn*aa*aa + mm*mm*bb*bb),
            Fr((aa*mm + bb*nn)**2),
            F4 - 2*mm*nn*aa*bb,
            Fr((aa*mm - bb*nn)**2),
            F4,
            Fr((aa*nn + bb*mm)**2),
            F4 + 2*mm*nn*aa*bb,
            Fr((aa*nn - bb*mm)**2),
            Fr(mm*mm*aa*aa + nn*nn*bb*bb)]

def five_Q(mm, nn, aa, bb):
    c = cells_Q(mm, nn, aa, bb)
    return [c[i] for i in (0, 2, 4, 6, 8)]

def is_rat_square(x):
    """точный тест «x — квадрат рационального»; isqrt, никаких float"""
    if x < 0:
        return False
    nu, de = x.numerator, x.denominator
    return isqrt(nu)**2 == nu and isqrt(de)**2 == de

def sq_or_zero_mod(x, p):
    """редукция Fraction в F_p: квадрат или ноль?  None — если знаменатель делится на p"""
    num, den = x.numerator, x.denominator
    if den % p == 0:
        return None
    v = (num % p) * pow(den % p, p - 2, p) % p
    return True if v == 0 else pow(v, (p - 1)//2, p) == 1

# ==============================================================================================
def section1_controls():
    head("§1. КОНТРОЛИ ПРОЦЕДУРЫ. Без них пустой перебор ничего не доказывает")
    print("  Логика: если у параметра есть ГЛОБАЛЬНОЕ решение пяти квадратностей, тест обязан")
    print("  пропускать его при КАЖДОМ p. Провал контроля означал бы, что тест режет настоящие решения.")
    PR = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113]
    ctrl = [("(7,1),  t=0    [вырожденная точка из REVIEW §1: (u0,u4,u8)=(7,5,1)]", 7, 1, 0, 1),
            ("(15,8), t=1    [вырожденная точка]",                                 15, 8, 1, 1),
            ("(15,8), t=-1   [вырожденная точка]",                                 15, 8, -1, 1),
            ("(7,1),  t=inf  [точка на бесконечности (a:b)=(1:0)]",                 7, 1, 1, 0),
            ("(1,7),  t=inf  [та же точка после симметрии m<->n]",                  1, 7, 1, 0)]
    ok = True
    for name, mm, nn, aa, bb in ctrl:
        vals = five_Q(mm, nn, aa, bb)
        allsq = all(is_rat_square(v) for v in vals)
        bad = [p for p in PR if not all(sq_or_zero_mod(v, p) for v in vals)]
        print(f"  {name}")
        print(f"     пять клеток над Q: {[str(v) for v in vals]}  -> все квадраты в Q: {allsq}")
        print(f"     редукция «квадрат-или-ноль» при всех p из списка: {not bad}"
              + (f"   ПРОВАЛ при p={bad}" if bad else ""))
        if allsq and bad:
            ok = False
    print("  ИТОГ §1:", "контроли пройдены" if ok else "!!! КОНТРОЛЬ ПРОВАЛЕН, всё ниже недействительно")
    return ok

# ==============================================================================================
def five_mod(mm, nn, p, aa, bb):
    m2, n2 = mm*mm % p, nn*nn % p
    mn2 = 2*mm*nn % p
    s = (m2 + n2) * pow(2, p-2, p) % p
    F0 = (n2*aa*aa + m2*bb*bb) % p
    F4 = s*(aa*aa + bb*bb) % p
    F8 = (m2*aa*aa + n2*bb*bb) % p
    return F0, F4, F8, (F4 - mn2*aa*bb) % p, (F4 + mn2*aa*bb) % p

def section2_main():
    head(f"§2. ГЛАВНОЕ: (m,n)=({M},{N}), p={P}. ПОЛНЫЙ перебор всех {P+1} точек P^1(F_{P})")
    m2, n2, mn2 = M*M % P, N*N % P, 2*M*N % P
    s = (m2 + n2) * pow(2, P-2, P) % P
    QR = sorted({x*x % P for x in range(1, P)})
    SQ0 = set(QR) | {0}
    print(f"  константы mod {P}: m^2={m2}, n^2={n2}, 2mn={mn2}, m^2+n^2={(m2+n2)%P}, s={s}")
    print(f"  ненулевые квадраты F_{P} ({len(QR)} шт.): {QR}")
    print(f"  сверка списка QR с Codex: "
          f"{QR == [1,2,4,5,8,9,10,16,18,20,21,23,25,31,32,33,36,37,39,40]}")
    print(f"  ВАЖНО: нуль РАЗРЕШЁН — ненулевой рациональный квадрат может редуцироваться в 0 mod p.")

    pts = [(t, 1) for t in range(P)] + [(1, 0)]
    assert len(pts) == P + 1
    rows, s3, s5 = [], [], []
    for (aa, bb) in pts:
        F0, F4, F8, L, U = five_mod(M, N, P, aa, bb)
        ok3 = all(v in SQ0 for v in (F0, F4, F8))
        ok5 = ok3 and L in SQ0 and U in SQ0
        rows.append(((aa, bb), (F0, F4, F8), (L, U), ok3, ok5))
        if ok3: s3.append((aa, bb))
        if ok5: s5.append((aa, bb))
    nm = lambda aa, bb: "inf" if bb == 0 else str(aa)
    print(f"\n  (а) после ТРЁХ квадратностей F0,F4,F8 (это кривая C7) выживает {len(s3)} точек: "
          f"{[nm(*q) for q in s3]}")
    print(f"      Codex: «остаются ровно четыре конечных параметра 4,10,31,37» -> совпадает: "
          f"{sorted(nm(*q) for q in s3) == sorted(['4','10','31','37'])}")
    print("\n  (б) выжившие тройки и что с ними делают красные клетки L=c2, U=c6:")
    print("          t |   F0   F4   F8 |    L    U | L кв? | U кв? | вердикт")
    for q, (F0, F4, F8), (L, U), ok3, ok5 in rows:
        if not ok3:
            continue
        print(f"      {nm(*q):>5} | {F0:4} {F4:4} {F8:4} | {L:4} {U:4} | "
              f"{str(L in SQ0):>5} | {str(U in SQ0):>5} | {'ПРОХОД' if ok5 else 'отказ'}")
    inf = five_mod(M, N, P, 1, 0)
    print(f"\n  (в) бесконечность (1:0): (F0,F4,F8,L,U) = {inf}; центр F4 = s = {s}, "
          f"квадрат-или-ноль: {s in SQ0}")
    print(f"\n  ИТОГ §2: точек P^1(F_{P}) со ВСЕМИ ПЯТЬЮ квадратностями: {len(s5)} из {len(pts)}")
    return len(s5), rows, SQ0, QR

# ==============================================================================================
def section3_independent():
    head("§3. ВТОРОЕ НЕЗАВИСИМОЕ КОДИРОВАНИЕ: все 1680 пар (a,b)!=(0,0) в F_41^2, ДЕВЯТЬ клеток,"
         "\n    квадратность — по таблице значений x^2 (символ Лежандра не используется)")
    sq_set = {x*x % P for x in range(P)}
    inv2 = pow(2, P-2, P)
    good = []
    auto_ok = True
    for aa in range(P):
        for bb in range(P):
            if aa == 0 and bb == 0:
                continue
            F4 = (M*M + N*N) % P * inv2 % P * ((aa*aa + bb*bb) % P) % P
            cs = [(N*N*aa*aa + M*M*bb*bb) % P,
                  (aa*M + bb*N)**2 % P,
                  (F4 - 2*M*N*aa*bb) % P,
                  (aa*M - bb*N)**2 % P,
                  F4,
                  (aa*N + bb*M)**2 % P,
                  (F4 + 2*M*N*aa*bb) % P,
                  (aa*N - bb*M)**2 % P,
                  (M*M*aa*aa + N*N*bb*bb) % P]
            if any(cs[i] not in sq_set for i in (1, 3, 5, 7)):
                auto_ok = False
            if all(c in sq_set for c in cs):
                good.append((aa, bb))
    print(f"  проверено пар (a,b): {P*P - 1}")
    print(f"  рёберные c1,c3,c5,c7 квадраты ВСЕГДА (контроль вывода §0): {auto_ok}")
    print(f"  пар, где ВСЕ ДЕВЯТЬ клеток — квадрат-или-ноль: {len(good)}")
    if good:
        print(f"  примеры: {good[:10]}")
    return len(good)

# ==============================================================================================
def section4_codex():
    head("§4. ПОСТРОЧНАЯ СВЕРКА С ТАБЛИЦЕЙ CODEX (PARALLEL_RESULT_FULL_NINE_126, «(16,5), p=41»)")
    codex = {4: ((0,31,21), (6,15)), 10: ((9,25,0), (24,26)),
             31: ((9,25,0), (26,24)), 37: ((0,31,21), (15,6))}
    SQ0 = {x*x % P for x in range(P)}
    allm = True
    for t in sorted(codex):
        F0, F4, F8, L, U = five_mod(M, N, P, t, 1)
        ok = (F0, F4, F8) == codex[t][0] and (L, U) == codex[t][1]
        allm &= ok
        print(f"  t={t:>2}: мой расчёт (F0,F4,F8)=({F0},{F4},{F8}), (L,U)=({L},{U}); "
              f"Codex {codex[t][0]}, {codex[t][1]} -> {'совпадает' if ok else 'РАСХОЖДЕНИЕ'}")
    s = (M*M + N*N) * pow(2, P-2, P) % P
    print(f"  Codex «на бесконечности F4 = s = 38 не квадрат»: s={s}, не квадрат: {s not in SQ0}")
    both = all(five_mod(M,N,P,t,1)[3] not in SQ0 and five_mod(M,N,P,t,1)[4] not in SQ0 for t in codex)
    print(f"  Codex «оба не квадраты» во всех четырёх строках: {both}")
    print(f"  ИТОГ §4: таблица Codex воспроизведена полностью: {allm and both}")
    return allm and both

# ==============================================================================================
def section5_degeneracy():
    head("§5. ВЫРОЖДЕНИЯ: опирается ли вывод на положительность, различность или запрет нулей")
    SQ0 = {x*x % P for x in range(P)}
    inv2 = pow(2, P-2, P)
    print("  (а) девять клеток mod 41 у четырёх точек, выживших после F0,F4,F8:")
    for t in (4, 10, 31, 37):
        F4 = (M*M + N*N) % P * inv2 % P * ((t*t + 1) % P) % P
        cs = [(N*N*t*t + M*M) % P, (t*M + N)**2 % P, (F4 - 2*M*N*t) % P, (t*M - N)**2 % P, F4,
              (t*N + M)**2 % P, (F4 + 2*M*N*t) % P, (t*N - M)**2 % P, (M*M*t*t + N*N) % P]
        zeros = [i for i in range(9) if cs[i] == 0]
        dups = [(i, j) for i in range(9) for j in range(i+1, 9) if cs[i] == cs[j]]
        print(f"      t={t:>2}: {cs}; нулевые клетки: {zeros or 'нет'}; совпадающих пар: {len(dups)}")
    print("  (б) в решающем переборе §2/§3 НЕТ ни одного сравнения c_i > 0 и ни одной проверки")
    print("      c_i != c_j. Тест — только «квадрат или ноль в F_p». Значит вывод от положительности")
    print("      и попарной различности НЕ зависит: они его могли бы лишь усилить.")
    strict = [q for q in [(t,1) for t in range(P)] + [(1,0)]
              if all(v in SQ0 and v != 0 for v in five_mod(M, N, P, *q))]
    perm = [q for q in [(t,1) for t in range(P)] + [(1,0)]
            if all(v in SQ0 for v in five_mod(M, N, P, *q))]
    print(f"  (в) счёт в обе стороны: при ЗАПРЕЩЁННЫХ нулях выживает {len(strict)} точек, "
          f"при РАЗРЕШЁННЫХ — {len(perm)}.")
    print(f"      Используется более слабый (разрешающий) вариант, поэтому пустота — более сильный факт.")
    print("  (г) явная проверка нескольких вырожденных t над Q (тест обязан их пропускать, если они есть):")
    for t in [Fr(0), Fr(1), Fr(-1), Fr(M, N), Fr(N, M), Fr(3, 2), Fr(-M, N)]:
        vals = five_Q(M, N, t.numerator, t.denominator)
        print(f"      t={str(t):>6}: {[str(v) for v in vals]} -> все квадраты в Q: "
              f"{all(is_rat_square(v) for v in vals)}")
    return len(perm)

# ==============================================================================================
def vp(x, p):
    if x == 0:
        return None
    nu, de = x.numerator, x.denominator
    v = 0
    while nu % p == 0: nu //= p; v += 1
    while de % p == 0: de //= p; v -= 1
    return v

def is_sq_Qp(x, p):
    v = vp(x, p)
    if v is None:
        return True
    if v % 2:
        return False
    y = x / Fr(p)**v
    u = (y.numerator % p) * pow(y.denominator % p, p-2, p) % p
    return pow(u, (p-1)//2, p) == 1

def solve_disk(mm, nn, p, r, k, KMAX, at_inf, budget):
    """есть ли t в диске r + p^k Z_p (или t = 1/u, u в диске) со всеми пятью клетками-квадратами в Q_p"""
    if budget[0] <= 0:
        return None
    budget[0] -= 1
    vals = five_Q(mm, nn, 1, r) if at_inf else five_Q(mm, nn, r, 1)
    undec = False
    for x in vals:
        v = vp(x, p)
        if v is not None and v < k:          # значение определяет класс квадратов на всём диске
            if not is_sq_Qp(x, p):
                return False
        else:
            undec = True
    if not undec:
        return True
    if k >= KMAX:
        return None
    res = False
    for z in range(p):
        out = solve_disk(mm, nn, p, r + z*p**k, k+1, KMAX, at_inf, budget)
        if out is True:
            return True
        if out is None:
            res = None
    return res

def section6_qp():
    head("§6. НЕЗАВИСИМЫЙ ПОДЪЁМ В Q_41 (рекурсия по дискам + Гензель), без ссылки на §2")
    bud = [300000]
    fin = solve_disk(M, N, P, 0, 0, 4, False, bud)
    inf = solve_disk(M, N, P, 0, 1, 4, True, bud)
    print(f"  ветка t in Z_{P}:            {fin}    (False = точек нет, None = НЕ РЕШЕНО)")
    print(f"  ветка t = 1/u, u in {P}Z_{P}:  {inf}")
    proved = (fin is False and inf is False)
    print(f"  ИТОГ §6: отсутствие Q_{P}-точки ДОКАЗАНО подъёмом: {proved}")
    # контроль самой процедуры: там, где точка есть, она обязана находиться
    ctrl = solve_disk(15, 8, P, 0, 0, 4, False, [300000])
    print(f"  контроль процедуры: (15,8) при p={P} (есть глобальная точка t=1) -> {ctrl} (ожидалось True)")
    return proved, ctrl

# ==============================================================================================
def empty_at(mm, nn, p, k=5):
    """пуста ли система из k клеток над F_p для пары (mm,nn); перебор всех p+1 точек P^1"""
    m2, n2, mn2 = mm*mm % p, nn*nn % p, 2*mm*nn % p
    s = (m2 + n2) * pow(2, p-2, p) % p
    S = {x*x % p for x in range(p)}
    for (aa, bb) in [(t, 1) for t in range(p)] + [(1, 0)]:
        F0 = (n2*aa*aa + m2*bb*bb) % p
        F4 = s*(aa*aa + bb*bb) % p
        F8 = (m2*aa*aa + n2*bb*bb) % p
        vals = (F0, F4, F8) if k == 3 else (F0, F4, F8, (F4 - mn2*aa*bb) % p, (F4 + mn2*aa*bb) % p)
        if all(v in S for v in vals):
            return False
    return True

def section7_extra():
    head("§7. ПОБОЧНОЕ: минимальное исключающее p; запрещённые классы λ = m/n; статус C7")
    PR = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
    k5 = [p for p in PR if empty_at(M, N, p, 5)]
    k3 = [p for p in PR if empty_at(M, N, p, 3)]
    print(f"  (16,5), p<=97: ПЯТЬ клеток локально невозможны при p = {k5 or 'нет'}")
    print(f"  (16,5), p<=97: ТРИ клетки (кривая C7) локально невозможны при p = {k3 or 'нет'}")
    print(f"  => минимальное исключающее p для девяти клеток: {k5[0] if k5 else 'нет'}"
          f"  (Codex заявлял 41: {bool(k5) and k5[0] == 41})")
    print(f"  => C7 локально НЕ исключается ни одним p<=97: {not k3}. Значит C7(Q)=пусто отсюда")
    print(f"     НЕ следует и не заявляется (различение из REVIEW §1 соблюдено).")
    for p, codex in ((41, [5,8,13,19,22,28,33,36]), (17, [2,3,6,8,9,11,14,15])):
        lam = [l for l in range(p) if empty_at(l, 1, p, 5)]
        ninf = empty_at(1, 0, p, 5)
        print(f"  p={p}: запрещённые λ=m/n (n!=0): {lam}; совпадает с Codex: {lam == codex}; "
              f"класс n=0 mod p исключён: {ninf}")
    lam = M * pow(N, P-2, P) % P
    print(f"  λ для нашей пары: {M}/{N} = {lam} mod {P}; в запрещённом списке: "
          f"{lam in [l for l in range(P) if empty_at(l, 1, P, 5)]}")
    return k5, k3

# ==============================================================================================
def section8_lambda():
    head("§8. ЗАПАС ПРОЧНОСТИ: ослабленный тест с общим множителем λ")
    print("  Над Q общий множитель λ обязан быть квадратом (рёберные клетки — квадраты линейных форм),")
    print("  поэтому λ-свободы на самом деле нет. Но mod p ребро может занулиться, и тогда λ формально")
    print("  не обязан быть квадратом. Проверяем заведомо БОЛЕЕ СЛАБЫЙ тест: существуют ли (a,b) и")
    print("  λ in F_41^*, при которых все девять λ*c_i — квадрат-или-ноль. Пустота здесь — запас прочности.")
    SQ = {x*x % P for x in range(P)}
    inv2 = pow(2, P-2, P)
    found9, found5 = [], []
    for aa in range(P):
        for bb in range(P):
            if aa == 0 and bb == 0:
                continue
            F4 = (M*M + N*N) % P * inv2 % P * ((aa*aa + bb*bb) % P) % P
            cs = [(N*N*aa*aa + M*M*bb*bb) % P, (aa*M + bb*N)**2 % P, (F4 - 2*M*N*aa*bb) % P,
                  (aa*M - bb*N)**2 % P, F4, (aa*N + bb*M)**2 % P, (F4 + 2*M*N*aa*bb) % P,
                  (aa*N - bb*M)**2 % P, (M*M*aa*aa + N*N*bb*bb) % P]
            for lam in range(1, P):
                if all((lam*c) % P in SQ for c in cs):
                    found9.append((aa, bb, lam)); break
            for lam in range(1, P):
                if all((lam*cs[i]) % P in SQ for i in (0, 2, 4, 6, 8)):
                    found5.append((aa, bb, lam)); break
    print(f"  пар (a,b,λ) с девятью λ*c_i квадратами-или-нулями: {len(found9)}")
    print(f"  пар (a,b,λ) с пятью неавтоматическими:             {len(found5)}")
    print(f"  ИТОГ §8: препятствие держится даже в ослабленной постановке: "
          f"{len(found9) == 0 and len(found5) == 0}")
    return len(found9), len(found5)

# ==============================================================================================
if __name__ == "__main__":
    ok0 = section0_sympy()
    ok1 = section1_controls()
    n5, rows, SQ0, QR = section2_main()
    n9 = section3_independent()
    ok4 = section4_codex()
    nperm = section5_degeneracy()
    proved6, ctrl6 = section6_qp()
    k5, k3 = section7_extra()
    l9, l5 = section8_lambda()

    head("ФИНАЛЬНЫЙ ВЕРДИКТ")
    print(f"  §0 вывод клеток из аксиом сошёлся с формулами проекта:      {ok0}")
    print(f"  §1 контроли процедуры пройдены:                             {ok1}")
    print(f"  §2 точек P^1(F_41) со всеми пятью квадратностями:           {n5}")
    print(f"  §3 пар (a,b) с девятью квадратами-или-нулями (2-е кодир.):  {n9}")
    print(f"  §4 таблица Codex воспроизведена построчно:                  {ok4}")
    print(f"  §6 отсутствие Q_41-точки доказано подъёмом:                 {proved6} "
          f"(контроль процедуры: {ctrl6})")
    print(f"  §8 ослабленный λ-тест тоже пуст (запас прочности):          {l9 == 0 and l5 == 0}")
    verdict = ok0 and ok1 and n5 == 0 and n9 == 0 and proved6 and ctrl6 is True
    print()
    print(f"  ЛОКАЛЬНОЕ ПРЕПЯТСТВИЕ (16,5) ПРИ p=41 ПОДТВЕРЖДЕНО: {verdict}")
    print("  Точная область: не существует (a:b) in P^1(Q_41), при котором все пять неавтоматических")
    print("  клеток — квадраты в Q_41; следовательно нет и рациональной точки полной девятки.")
    print("  Это НЕ утверждение «C7(Q) пусто» и НЕ утверждение о задаче магического квадрата вообще.")
