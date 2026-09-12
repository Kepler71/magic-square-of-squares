# -*- coding: utf-8 -*-
# СТАТУС: проверка чужого утверждения — прямой перебор, без теории
# ДЛЯ: Claude/Codex/пользователь
# ИТОГ: независимая проверка локального исключения (m,n)=(16,5) по модулю 41
# ОТМЕНЯЕТ: ничего
# ПРОВЕРЕНО: девять клеток выведены символически заново; QR-таблицы построены возведением в квадрат
# ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: нужен воспроизводимый счёт всех p+1 проективных случаев
#
# Запуск: python3 /home/kep/magicKube/bridge/local_16_5_p41_прямая.py
#
# ЧТО ЗДЕСЬ ЕСТЬ И ЧЕГО НЕТ.
# Есть: (1) самостоятельный символический вывод девяти клеток G1 и проверка всех восьми магических
#           линий = 3*центр (sympy, тождества в Q(m,n,t));
#       (2) независимое построение множества квадратичных вычетов возведением в квадрат всех
#           residues (не через символ Лежандра), затем сверка с pow(x,(p-1)//2);
#       (3) ПОЛНЫЙ перебор всех p+1 точек P^1(F_p), включая (1:0);
#       (4) контроли: пары/параметры с ЗАВЕДОМО существующим глобальным решением пяти квадратностей
#           обязаны выживать при каждом p — если контроль падает, падает вся процедура.
# Нет: ни положительности, ни попарной различности клеток нигде не используется.
#      Нули при редукции разрешены (ненулевой рациональный квадрат может редуцироваться в 0).

from fractions import Fraction as Fr
from math import gcd
import sys

# ----------------------------------------------------------------------------------
# 0. Девять клеток семейства G1 — выводим и проверяем заново
# ----------------------------------------------------------------------------------
# Однородная (a:b) форма, t = a/b; каждая клетка домножена на b^2 (квадрат => квадратность не меняется)
#
#   F0 = m^2 b^2 + n^2 a^2        c1 = (n b + m a)^2        L = F4 - 2 m n a b
#   c3 = (m a - n b)^2            F4 = s (a^2 + b^2)        c5 = (m b + n a)^2
#   U  = F4 + 2 m n a b           c7 = (m b - n a)^2        F8 = n^2 b^2 + m^2 a^2
#   s = (m^2 + n^2)/2
#
# Раскладка по позициям магического квадрата:
#      c0 c1 c2          F0  c1  L
#      c3 c4 c5    =     c3  F4  c5
#      c6 c7 c8          U   c7  F8
# Неавтоматические квадратности ровно пять: F0, F4, F8, L, U.

def cells_homog(m, n, a, b, half=Fr(1, 2)):
    """девять клеток как однородные квадратичные формы от (a,b). half = 1/2 в нужном кольце."""
    s = half * (m * m + n * n)
    F4 = s * (a * a + b * b)
    c0 = m * m * b * b + n * n * a * a
    c1 = (n * b + m * a) ** 2
    c2 = F4 - 2 * m * n * a * b
    c3 = (m * a - n * b) ** 2
    c4 = F4
    c5 = (m * b + n * a) ** 2
    c6 = F4 + 2 * m * n * a * b
    c7 = (m * b - n * a) ** 2
    c8 = n * n * b * b + m * m * a * a
    return [c0, c1, c2, c3, c4, c5, c6, c7, c8]

NONAUTO = [0, 4, 8, 2, 6]          # F0, F4, F8, L, U
NAMES = {0: "F0", 4: "F4", 8: "F8", 2: "L", 6: "U"}
LINES = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]


def symbolic_check():
    print("=" * 100)
    print("0. САМОСТОЯТЕЛЬНЫЙ СИМВОЛИЧЕСКИЙ ВЫВОД: магичность и автоматические квадраты")
    print("=" * 100)
    try:
        import sympy as sp
    except ImportError:
        print("  sympy нет — символическая часть пропущена (НЕ считать проверенной)")
        return False
    m, n, a, b = sp.symbols("m n a b")
    C = cells_homog(m, n, a, b, half=sp.Rational(1, 2))
    ok = True
    for L in LINES:
        d = sp.simplify(C[L[0]] + C[L[1]] + C[L[2]] - 3 * C[4])
        good = (d == 0)
        ok &= good
        print(f"  линия {L}: сумма - 3*центр = {d}   {'OK' if good else 'ОШИБКА'}")
    print("  --- какие клетки — квадраты ТОЖДЕСТВЕННО (как многочлены от a,b,m,n):")

    def is_perfect_square_poly(expr):
        """честный тест: многочлен — полный квадрат в Q[a,b,m,n] <=> в разложении
           все неконстантные множители входят в чётной степени, а константа — квадрат в Q."""
        expr = sp.together(expr)
        num, den = sp.fraction(expr)
        const, facs = sp.factor_list(sp.expand(num))
        for f, e in facs:
            if e % 2:
                return False
        c = sp.Rational(const) / sp.Rational(den)
        return sp.sqrt(c).is_rational is True

    for i in range(9):
        auto = is_perfect_square_poly(C[i])
        print(f"    c{i} = {sp.factor(C[i])}    автоматический квадрат: {auto}")
        if i in NONAUTO and auto:
            print("      !!! неожиданно: неавтоматическая клетка оказалась квадратом")
            ok = False
        if i not in NONAUTO and not auto:
            print("      !!! неожиданно: клетка, считавшаяся автоматической, не квадрат")
            ok = False
    # тождества-контроли
    for name, expr in [("c2+c6-2*c4", C[2] + C[6] - 2 * C[4]),
                       ("2*c2-c3-c7", 2 * C[2] - C[3] - C[7]),
                       ("2*c6-c1-c5", 2 * C[6] - C[1] - C[5])]:
        d = sp.simplify(expr)
        print(f"  тождество {name} = {d}  {'OK' if d == 0 else 'ОШИБКА'}")
        ok &= (d == 0)
    print(f"  ИТОГ раздела 0: {'все проверки пройдены' if ok else 'ЕСТЬ РАСХОЖДЕНИЯ'}")
    return ok


# ----------------------------------------------------------------------------------
# 1. Квадратичные вычеты — строим перебором, сверяем с критерием Эйлера
# ----------------------------------------------------------------------------------
def qr_set(p):
    """множество НЕНУЛЕВЫХ квадратов в F_p, построенное возведением в квадрат"""
    S = set((x * x) % p for x in range(1, p))
    # независимая сверка критерием Эйлера
    T = set(x for x in range(1, p) if pow(x, (p - 1) // 2, p) == 1)
    assert S == T, f"QR-таблицы разошлись при p={p}"
    return S


def is_sq_mod(x, p, QR, allow_zero=True):
    x %= p
    if x == 0:
        return allow_zero
    return x in QR


# ----------------------------------------------------------------------------------
# 2. Перебор всех p+1 точек P^1(F_p)
# ----------------------------------------------------------------------------------
def projective_points(p):
    """все p+1 точек P^1(F_p): (t:1) для t=0..p-1 и (1:0)"""
    return [(t, 1) for t in range(p)] + [(1, 0)]


def cells_mod(m, n, a, b, p):
    inv2 = pow(2, p - 2, p)
    s = ((m * m + n * n) % p) * inv2 % p
    F4 = s * ((a * a + b * b) % p) % p
    c0 = (m * m * b * b + n * n * a * a) % p
    c1 = pow(n * b + m * a, 2, p)
    c2 = (F4 - 2 * m * n * a * b) % p
    c3 = pow(m * a - n * b, 2, p)
    c4 = F4
    c5 = pow(m * b + n * a, 2, p)
    c6 = (F4 + 2 * m * n * a * b) % p
    c7 = pow(m * b - n * a, 2, p)
    c8 = (n * n * b * b + m * m * a * a) % p
    return [c0, c1, c2, c3, c4, c5, c6, c7, c8]


def scan(m, n, p, idx=NONAUTO, allow_zero=True, verbose=False):
    """возвращает список выживших проективных точек и подробную таблицу"""
    QR = qr_set(p)
    survivors, table = [], []
    for (a, b) in projective_points(p):
        C = cells_mod(m, n, a, b, p)
        flags = [is_sq_mod(C[i], p, QR, allow_zero) for i in idx]
        row = {"pt": (a, b), "cells": C,
               "vals": {NAMES.get(i, f"c{i}"): C[i] for i in idx},
               "flags": dict(zip([NAMES.get(i, f"c{i}") for i in idx], flags)),
               "ok": all(flags)}
        table.append(row)
        if row["ok"]:
            survivors.append((a, b))
    if verbose:
        for row in table:
            if row["ok"]:
                print("   ВЫЖИЛ", row["pt"], row["vals"])
    return survivors, table


def label(pt, p):
    a, b = pt
    return "inf" if b == 0 else str(a % p)


# ----------------------------------------------------------------------------------
# 3. Контроли: заведомые глобальные решения пяти квадратностей
# ----------------------------------------------------------------------------------
# (7,1), t=0   : c0=49, c4=25, c8=1, c2=25, c6=25 — все пять ЦЕЛЫЕ КВАДРАТЫ (вырожденно)
# (15,8), t=1  : c0=c4=c8=289=17^2, c2=49, c6=529 — все пять квадраты
# (15,8), t=-1 : то же зеркально
# (1,7),  t=inf: проверка ветки бесконечности (симметрия m<->n, t<->1/t)
CONTROLS = [((7, 1), (0, 1), "t=0"), ((15, 8), (1, 1), "t=1"), ((15, 8), (-1, 1), "t=-1"),
            ((1, 7), (1, 0), "t=inf")]


def check_controls(primes):
    print("=" * 100)
    print("1. КОНТРОЛИ ПРОЦЕДУРЫ. Параметры с ГЛОБАЛЬНЫМ решением пяти квадратностей")
    print("   обязаны выживать при КАЖДОМ простом. Падение контроля = процедура негодна.")
    print("=" * 100)
    all_ok = True
    for (m, n), (a, b), tag in CONTROLS:
        # сначала честно проверим над Q, что все пять — квадраты рациональных
        if b != 0:
            C = cells_homog(m, n, a, b)
            vals = [C[i] for i in NONAUTO]
            def issq(x):
                if x < 0: return False
                r = Fr(x).limit_denominator() if isinstance(x, Fr) else x
                num, den = Fr(x).numerator, Fr(x).denominator
                return int(num ** 0.5 + 0.5) ** 2 == num and int(den ** 0.5 + 0.5) ** 2 == den
            glob = all(issq(v) for v in vals)
            print(f"  ({m},{n}) {tag}: пять клеток над Q = {[str(v) for v in vals]}; все квадраты: {glob}")
            if not glob:
                print("    !!! контрольная точка не является глобальным решением — контроль недействителен")
                all_ok = False
                continue
        else:
            print(f"  ({m},{n}) {tag}: ветка бесконечности, F0=n^2={n*n}, F8=m^2={m*m}, "
                  f"F4=L=U=s={Fr(m*m+n*n,2)}")
        bad = []
        for p in primes:
            if p == 2:
                continue
            QR = qr_set(p)
            C = cells_mod(m, n, a % p, b % p, p)
            if not all(is_sq_mod(C[i], p, QR) for i in NONAUTO):
                bad.append(p)
        print(f"     редукция даёт квадрат-или-ноль при всех p: {'ДА' if not bad else 'НЕТ, сбой при ' + str(bad)}")
        all_ok &= (not bad)
    print(f"  ИТОГ раздела 1: {'контроли пройдены' if all_ok else 'КОНТРОЛЬ УПАЛ'}")
    return all_ok


# ----------------------------------------------------------------------------------
# 4. Главное: (16,5) по модулю 41
# ----------------------------------------------------------------------------------
def main_case():
    m, n, p = 16, 5, 41
    QR = qr_set(p)
    inv2 = pow(2, p - 2, p)
    print("=" * 100)
    print(f"2. ГЛАВНОЕ: (m,n)=({m},{n}), p={p}. ПОЛНЫЙ перебор всех {p+1} точек P^1(F_{p})")
    print("=" * 100)
    print(f"  константы mod {p}: m^2={m*m%p}, n^2={n*n%p}, 2mn={2*m*n%p}, "
          f"s=(m^2+n^2)/2={((m*m+n*n)%p)*inv2%p}, m^2+n^2={(m*m+n*n)%p}")
    print(f"  QR({p}) без нуля, {len(QR)} штук: {sorted(QR)}")
    print(f"  сверка со списком Codex: {sorted(QR | {0}) == sorted([0,1,2,4,5,8,9,10,16,18,20,21,23,25,31,32,33,36,37,39,40])}")
    print()

    # 2a. сколько выживает после ТРЁХ белых квадратностей (кривая C7)
    surv3, tab3 = scan(m, n, p, idx=[0, 4, 8])
    print(f"  (a) после трёх квадратностей F0,F4,F8 выживает {len(surv3)} точек: "
          f"{[label(x,p) for x in surv3]}")
    print(f"      сверка с Codex «остаются ровно четыре конечных параметра {{4,10,31,37}}»: "
          f"{sorted(label(x,p) for x in surv3) == sorted(['4','10','31','37'])}")
    print()
    print("  (b) подробная таблица выживших троек и что с ними делают красные L,U:")
    print(f"      {'t':>5} | {'F0':>4} {'F4':>4} {'F8':>4} | {'L':>4} {'U':>4} | L кв? U кв? | вердикт")
    for (a, b) in surv3:
        C = cells_mod(m, n, a, b, p)
        L, U = C[2], C[6]
        lo, uo = is_sq_mod(L, p, QR), is_sq_mod(U, p, QR)
        print(f"      {label((a,b),p):>5} | {C[0]:>4} {C[4]:>4} {C[8]:>4} | {L:>4} {U:>4} | "
              f"{str(lo):>5} {str(uo):>5} | {'ПРОХОДИТ' if lo and uo else 'отказ'}")
    print()

    # 2b. полный перебор пяти квадратностей
    surv5, tab5 = scan(m, n, p, idx=NONAUTO)
    print(f"  (c) ПОЛНЫЙ перебор пяти квадратностей по всем {p+1} проективным точкам:")
    print(f"      выживших точек: {len(surv5)}  -> {[label(x,p) for x in surv5] if surv5 else 'НИ ОДНОЙ'}")
    # явный пересчёт, сколько отказов на каждой клетке
    from collections import Counter
    cnt = Counter()
    for row in tab5:
        first = [k for k, v in row["flags"].items() if not v]
        cnt[tuple(first)] += 1
    print("      разбиение всех точек по множеству НЕвыполненных квадратностей:")
    for k, v in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"        {('все выполнены' if not k else 'не квадраты: ' + ','.join(k)):<40} {v:>3} точек")
    print()

    # 2c. точка на бесконечности отдельно
    C = cells_mod(m, n, 1, 0, p)
    print(f"  (d) точка (1:0) (t=infinity) отдельно: F0={C[0]} (кв:{is_sq_mod(C[0],p,QR)}), "
          f"F4={C[4]} (кв:{is_sq_mod(C[4],p,QR)}), F8={C[8]} (кв:{is_sq_mod(C[8],p,QR)}), "
          f"L={C[2]} (кв:{is_sq_mod(C[2],p,QR)}), U={C[6]} (кв:{is_sq_mod(C[6],p,QR)})")
    print()

    # 2d. зависит ли вывод от разрешения нулей / от различности / положительности
    surv5_nz, _ = scan(m, n, p, idx=NONAUTO, allow_zero=False)
    print("  (e) чувствительность вывода к вырожденностям:")
    print(f"      нули РАЗРЕШЕНЫ (самый мягкий критерий, честный):   выживших {len(surv5)}")
    print(f"      нули ЗАПРЕЩЕНЫ (более жёсткий, необоснованный):    выживших {len(surv5_nz)}")
    print("      => вывод получен при САМОМ МЯГКОМ критерии, значит запрет нулей на него не влияет")
    # проверим, встречаются ли нули и совпадения среди тех, кто дошёл до красных
    print("      вырожденности среди четырёх дошедших до красных клеток:")
    for (a, b) in surv3:
        C = cells_mod(m, n, a, b, p)
        zeros = [f"c{i}" for i in range(9) if C[i] % p == 0]
        eqs = [f"c{i}=c{j}" for i in range(9) for j in range(i + 1, 9) if C[i] % p == C[j] % p]
        print(f"        t={label((a,b),p):>3}: нулевые клетки mod p: {zeros or 'нет'}; "
              f"совпадения mod p: {len(eqs)} пар")
    print("      (совпадения КЛЕТОК mod p ничего не запрещают и в критерии не используются)")
    print()

    # 2e. независимый маршрут: перебор по ДЕВЯТИ клеткам целиком
    surv9, _ = scan(m, n, p, idx=list(range(9)))
    print(f"  (f) независимый маршрут — требуем квадратности ВСЕХ ДЕВЯТИ клеток (а не пяти): "
          f"выживших {len(surv9)}. Совпадает с пятиклеточным: {len(surv9) == len(surv5)}")
    return len(surv5) == 0


# ----------------------------------------------------------------------------------
# 5. Повторение таблицы Codex для (13,8) mod 17 — как второй контроль его метода
# ----------------------------------------------------------------------------------
def case_13_8():
    m, n, p = 13, 8, 17
    QR = qr_set(p)
    print("=" * 100)
    print("3. КОНТРОЛЬ ЧУЖОЙ ТАБЛИЦЫ: (13,8) mod 17 — воспроизводим независимо")
    print("=" * 100)
    surv3, _ = scan(m, n, p, idx=[0, 4, 8])
    print(f"  после F0,F4,F8 выживает {len(surv3)}: {[label(x,p) for x in surv3]}; "
          f"ожидалось {{2,8,9,15}}: {sorted(label(x,p) for x in surv3) == sorted(['2','8','9','15'])}")
    for (a, b) in surv3:
        C = cells_mod(m, n, a, b, p)
        print(f"    t={label((a,b),p):>3}: (F0,F4,F8)=({C[0]},{C[4]},{C[8]}) (L,U)=({C[2]},{C[6]}) "
              f"L кв:{is_sq_mod(C[2],p,QR)} U кв:{is_sq_mod(C[6],p,QR)}")
    surv5, _ = scan(m, n, p, idx=NONAUTO)
    print(f"  полный перебор 18 точек: выживших {len(surv5)}")
    return len(surv5) == 0


# ----------------------------------------------------------------------------------
# 6. Насколько 41 особое: какие простые ловят (16,5), и что с остальными парами
# ----------------------------------------------------------------------------------
PRIMES = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113]


def context():
    print("=" * 100)
    print("4. КОНТЕКСТ")
    print("=" * 100)
    print("  (a) какие простые дают локальное препятствие для (16,5) на пяти клетках:")
    hit = [p for p in PRIMES if not scan(16, 5, p, idx=NONAUTO)[0]]
    print(f"      {hit}")
    hit3 = [p for p in PRIMES if not scan(16, 5, p, idx=[0, 4, 8])[0]]
    print(f"  (b) те же простые на ТРЁХ клетках (кривая C7): {hit3 or 'нет ни одного'}")
    print("      => если тут пусто, то препятствие живёт именно на красных клетках L,U")
    print()
    print("  (c) три пары, остававшиеся открытыми, на пяти клетках по всем p <= 113:")
    for (a, b) in [(13, 8), (15, 8), (16, 5)]:
        mm, nn = (b, a) if a % 2 == 0 else (a, b)
        h = [p for p in PRIMES if not scan(mm, nn, p, idx=NONAUTO)[0]]
        print(f"      ({a},{b}) -> (m,n)=({mm},{nn}): препятствие при p = {h or 'НЕТ НИ ПРИ ОДНОМ'}")
    print()
    print("  (d) запрещённые классы lambda = m/n mod p (сверка обобщения Codex):")
    for p, expect in [(17, [2,3,6,8,9,11,14,15]), (41, [5,8,13,19,22,28,33,36])]:
        bad = []
        for lam in range(p):
            if not scan(lam, 1, p, idx=NONAUTO)[0]:
                bad.append(lam)
        # плюс класс n = 0 mod p, то есть (m:n) = (1:0)
        n0 = not scan(1, 0, p, idx=NONAUTO)[0]
        print(f"      p={p}: запрещены lambda = {bad}; ожидалось {expect}; совпало: {bad == expect}")
        print(f"             класс n=0 mod {p} запрещён: {n0} (Codex: не входит в список)")
        lam_165 = (16 * pow(5, p - 2, p)) % p if p == 41 else None
        if lam_165 is not None:
            print(f"             lambda для (16,5) при p=41: {lam_165}; в запрещённом списке: {lam_165 in bad}")


# ----------------------------------------------------------------------------------
# 7. ДВА НЕЗАВИСИМЫХ МАРШРУТА (другой код, другая постановка)
# ----------------------------------------------------------------------------------
def independent_bruteforce(m=16, n=5, p=41):
    """Маршрут 1: не пользуемся ни P^1, ни таблицей QR. Перебираем ВСЕ (a,b) != (0,0)
       в F_p^2 и ищем u0,u4,u8,uL,uU в F_p с u^2 = соответствующая форма.
       Это буквально вопрос «есть ли точка аффинной схемы над F_p», решённый грубой силой."""
    print("=" * 100)
    print("5. НЕЗАВИСИМЫЙ МАРШРУТ 1: грубый перебор всех (a,b) в F_p^2 и всех корней u_i")
    print("=" * 100)
    sq = {}
    for u in range(p):
        sq.setdefault((u * u) % p, []).append(u)          # значение -> список корней
    found = []
    total = 0
    for a in range(p):
        for b in range(p):
            if a == 0 and b == 0:
                continue
            total += 1
            C = cells_mod(m, n, a, b, p)
            vals = [C[i] for i in NONAUTO]
            if all(v in sq for v in vals):
                found.append(((a, b), vals, [sq[v][0] for v in vals]))
    print(f"  перебрано пар (a,b): {total} (= p^2-1 = {p*p-1})")
    print(f"  наборов (a,b,u0,u4,u8,uL,uU) над F_{p} с ВСЕМИ пятью уравнениями: {len(found)}")
    if found:
        for f in found[:10]:
            print("    ", f)
    print(f"  => множество F_{p}-точек пятиквадратной схемы (кроме a=b=0) {'ПУСТО' if not found else 'НЕ ПУСТО'}")
    # контроль этого же кода: (15,8) при p=41 обязан дать непустой ответ (есть глобальная точка t=1)
    ctrl = []
    for a in range(p):
        for b in range(p):
            if a == 0 and b == 0:
                continue
            C = cells_mod(15, 8, a, b, p)
            if all(C[i] in sq for i in NONAUTO):
                ctrl.append((a, b))
    print(f"  КОНТРОЛЬ тем же кодом: (15,8) mod 41 -> найдено {len(ctrl)} пар (должно быть > 0, "
          f"есть глобальная точка t=1): {'OK' if ctrl else 'СБОЙ КОНТРОЛЯ'}")
    return len(found) == 0


def lemma_empirical_test(trials=4000, seed=20260912):
    """Маршрут 2: эмпирическая проверка самой ЛЕММЫ, на которой всё держится:
       если клетка c_i(t) — квадрат рационального числа, то её однородное значение
       F_i(a,b) при t = a/b в несократимом виде редуцируется в QR_p u {0}.
       Берём случайные (m,n,t), смотрим на те клетки, что реально оказались квадратами."""
    import random
    print("=" * 100)
    print("6. НЕЗАВИСИМЫЙ МАРШРУТ 2: эмпирическая проверка леммы о редукции")
    print("=" * 100)
    rnd = random.Random(seed)
    import math
    checked = 0
    viol = 0
    for _ in range(trials):
        m = rnd.randint(1, 60); n = rnd.randint(1, 60)
        if gcd(m, n) != 1:
            continue
        a = rnd.randint(-80, 80); b = rnd.randint(1, 80)
        g = gcd(abs(a), b)
        if g: a //= g; b //= g
        C = cells_homog(m, n, a, b)          # однородные значения, тип Fraction
        for i in NONAUTO:
            v = Fr(C[i])
            if v <= 0:
                continue
            num, den = v.numerator, v.denominator
            rn, rd = math.isqrt(num), math.isqrt(den)
            if rn * rn != num or rd * rd != den:
                continue                      # клетка не квадрат — лемма про неё молчит
            for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
                if den % p == 0:
                    continue                  # знаменатель делится на p -> значение не p-целое
                QR = qr_set(p)
                x = (num % p) * pow(den % p, p - 2, p) % p
                checked += 1
                if not is_sq_mod(x, p, QR):
                    viol += 1
                    print(f"    НАРУШЕНИЕ: (m,n)=({m},{n}) t={a}/{b} клетка {NAMES[i]} = {v}, p={p}, red={x}")
    print(f"  проверено пар (клетка-квадрат, простое): {checked}; нарушений леммы: {viol}")
    print(f"  => лемма {'подтверждена численно' if viol == 0 else 'ОПРОВЕРГНУТА'}")
    return viol == 0


def p47_extra():
    print("=" * 100)
    print("7. ВТОРОЕ, НЕЗАВИСИМОЕ ОТ 41 ПРЕПЯТСТВИЕ ДЛЯ (16,5): p = 47")
    print("=" * 100)
    for (m, n) in [(16, 5), (5, 16)]:
        p = 47
        QR = qr_set(p)
        s3, _ = scan(m, n, p, idx=[0, 4, 8])
        s5, _ = scan(m, n, p, idx=NONAUTO)
        print(f"  (m,n)=({m},{n}), p=47: после F0,F4,F8 выживает {len(s3)} точек "
              f"{[label(x,p) for x in s3]}; после всех пяти — {len(s5)}")
        for pt in s3:
            C = cells_mod(m, n, pt[0], pt[1], p)
            print(f"     t={label(pt,p):>4}: (F0,F4,F8)=({C[0]},{C[4]},{C[8]})  (L,U)=({C[2]},{C[6]})"
                  f"  L кв:{is_sq_mod(C[2],p,QR)}  U кв:{is_sq_mod(C[6],p,QR)}")


if __name__ == "__main__":
    ok0 = symbolic_check(); print()
    ok1 = check_controls(PRIMES); print()
    ok2 = main_case(); print()
    ok3 = case_13_8(); print()
    context(); print()
    print("=" * 100)
    print("ИТОГ")
    print("=" * 100)
    print(f"  символический вывод девяти клеток и магичность .......... {'OK' if ok0 else 'СБОЙ'}")
    print(f"  контроли (глобальные решения выживают везде) ............ {'OK' if ok1 else 'СБОЙ'}")
    print(f"  (16,5) mod 41: пяти­клеточная система ПУСТА ............. {'ДА' if ok2 else 'НЕТ'}")
    print(f"  (13,8) mod 17: пяти­клеточная система ПУСТА ............. {'ДА' if ok3 else 'НЕТ'}")
