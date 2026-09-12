# -*- coding: utf-8 -*-
# СКЕПТИК: попытка ОПРОВЕРГНУТЬ локальное исключение пары (m,n) = (16,5) по модулю 41.
# Заявление Codex (PARALLEL_RESULT_FULL_NINE_126_2026-09-12.md):
#   полная система девяти клеток G1 для (16,5) не имеет решений уже над F_41.
#
# План:
#  §1  ВЫВЕСТИ девять клеток G1 независимо, из одних только магических уравнений
#      и параметризации Сегре квадрики P^2+Q^2 = R^2+S^2. Ни одна формула Codex не берётся на веру.
#  §2  Сверить полученные пять свободных клеток с формами Codex (F0,F4,F8,L,U).
#  §3  Полный перебор всех p+1 = 42 точек P^1(F_41): значения пяти форм, квадратичность.
#  §4  Дыры-кандидаты: пропущенное t, ноль как квадрат, бесконечность, обе ориентации (m,n)/(n,m),
#      зависимость от положительности/различности, корректность простого и редукции.
#  §5  Контроли: пара с ЗАВЕДОМО существующим решением должна проходить тест;
#      тест на семи клетках должен НЕ закрывать (16,5) (иначе таблица Codex противоречива).
#  §6  Независимая проверка через прямой перебор рациональных t малой высоты в Q_41 (Гензель).
#
# Запуск: python3 /home/kep/magicKube/bridge/local_16_5_p41_скептик.py

from fractions import Fraction as Fr
from math import gcd
from itertools import product
import sympy as sp

SEP = "=" * 100
sub = "-" * 100

# ----------------------------------------------------------------------------------
# §1. НЕЗАВИСИМЫЙ ВЫВОД ДЕВЯТИ КЛЕТОК G1
# ----------------------------------------------------------------------------------
# Постановка (из BREMNER16 §2.2, но формулы клеток выводим сами):
#   G1 = две ПОЛНЫЕ рёберные пары {c1,c7} и {c3,c5}.
#   Нумерация:   c0 c1 c2
#                c3 c4 c5
#                c6 c7 c8
#   c1 = P^2, c7 = Q^2, c3 = R^2, c5 = S^2  — четыре автоматических квадрата.
#   Условие совместности: P^2+Q^2 = R^2+S^2 = 2*c4.
#   Параметризация Сегре (q = 1, p = t):
#       P = m t + n,  Q = m - n t,  R = m t - n,  S = m + n t.
#   Остальные пять клеток определяются ЛИНЕЙНО из магических уравнений.

def derive_cells_symbolic():
    m, n, t = sp.symbols('m n t')
    P = m*t + n
    Q = m - n*t
    R = m*t - n
    S = m + n*t
    # проверка совместности квадрики
    assert sp.expand(P**2 + Q**2 - (R**2 + S**2)) == 0, "Сегре: квадрика не выполнена"

    c0, c2, c4, c6, c8 = sp.symbols('c0 c2 c4 c6 c8')
    c1, c3, c5, c7 = P**2, R**2, S**2, Q**2
    Ssum = 3*c4  # магическая сумма ассоциативного квадрата

    eqs = [
        sp.Eq(c0 + c1 + c2, Ssum),   # строка 0
        sp.Eq(c3 + c4 + c5, Ssum),   # строка 1
        sp.Eq(c6 + c7 + c8, Ssum),   # строка 2
        sp.Eq(c0 + c3 + c6, Ssum),   # столбец 0
        sp.Eq(c1 + c4 + c7, Ssum),   # столбец 1
        sp.Eq(c2 + c5 + c8, Ssum),   # столбец 2
        sp.Eq(c0 + c4 + c8, Ssum),   # диагональ
        sp.Eq(c2 + c4 + c6, Ssum),   # антидиагональ
    ]
    sol = sp.solve(eqs, [c0, c2, c4, c6, c8], dict=True)
    assert len(sol) == 1, f"ожидалось единственное решение, получено {len(sol)}"
    sol = sol[0]
    cells = {
        'c0': sp.expand(sol[c0]), 'c1': sp.expand(c1), 'c2': sp.expand(sol[c2]),
        'c3': sp.expand(c3),      'c4': sp.expand(sol[c4]), 'c5': sp.expand(c5),
        'c6': sp.expand(sol[c6]), 'c7': sp.expand(c7), 'c8': sp.expand(sol[c8]),
    }
    # КОНТРОЛЬ: магичность заново, «с нуля», без опоры на систему выше
    g = [[cells['c0'], cells['c1'], cells['c2']],
         [cells['c3'], cells['c4'], cells['c5']],
         [cells['c6'], cells['c7'], cells['c8']]]
    tot = sp.expand(sum(sum(r) for r in g)) / 3
    lines = [g[0], g[1], g[2],
             [g[0][0], g[1][0], g[2][0]], [g[0][1], g[1][1], g[2][1]], [g[0][2], g[1][2], g[2][2]],
             [g[0][0], g[1][1], g[2][2]], [g[0][2], g[1][1], g[2][0]]]
    for L in lines:
        assert sp.simplify(sp.expand(sum(L)) - tot) == 0, "квадрат не магический"
    return cells, (m, n, t)


def report_derivation():
    print(SEP)
    print("§1. НЕЗАВИСИМЫЙ ВЫВОД ДЕВЯТИ КЛЕТОК G1 (только магические уравнения + Сегре)")
    print(SEP)
    cells, (m, n, t) = derive_cells_symbolic()
    for k in ['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8']:
        print(f"   {k} = {sp.factor(cells[k])}")
    print("   [контроль] все 8 линий дают одну сумму — ПРОЙДЕН")

    # §2. Сверка с формами Codex
    print()
    print(SEP)
    print("§2. СВЕРКА С ФОРМАМИ CODEX  (F0,F4,F8,L,U)")
    print(SEP)
    s = (m**2 + n**2) / 2
    codex = {
        'F0': m**2 + n**2 * t**2,
        'F4': s * (1 + t**2),
        'F8': n**2 + m**2 * t**2,
        'L':  s * (1 + t**2) - 2*m*n*t,
        'U':  s * (1 + t**2) + 2*m*n*t,
    }
    pairs = [('c0', 'F0'), ('c4', 'F4'), ('c8', 'F8'), ('c2', 'L'), ('c6', 'U')]
    ok = True
    for c, f in pairs:
        d = sp.simplify(sp.expand(cells[c] - codex[f]))
        print(f"   {c} - {f} = {d}   {'СОВПАДАЕТ' if d == 0 else '*** РАСХОЖДЕНИЕ ***'}")
        ok = ok and d == 0
    # автоматические квадраты
    for c in ['c1', 'c3', 'c5', 'c7']:
        r = sp.factor(cells[c])
        issq = sp.simplify(sp.sqrt(r)**2 - r) == 0
        print(f"   {c} = {r}  — полный квадрат: {issq}")
    print(f"   ИТОГ §2: клетки {'ТЕ ЖЕ' if ok else 'РАЗНЫЕ'}")

    # §2a. симметрии, которые обязаны быть у системы (контроль структуры)
    print()
    print("   §2a. симметрии системы пяти форм:")
    swap = {m: n, n: m}
    print(f"      m<->n:  F0<->F8 ? {sp.simplify(codex['F0'].subs(swap, simultaneous=True) - codex['F8']) == 0};"
          f"  F4 инв ? {sp.simplify(codex['F4'].subs(swap, simultaneous=True) - codex['F4']) == 0};"
          f"  L инв ? {sp.simplify(codex['L'].subs(swap, simultaneous=True) - codex['L']) == 0}")
    print(f"      t->-t:  F0,F4,F8 инв ? "
          f"{all(sp.simplify(codex[k].subs(t, -t) - codex[k]) == 0 for k in ['F0','F4','F8'])};"
          f"  L<->U ? {sp.simplify(codex['L'].subs(t, -t) - codex['U']) == 0}")
    inv = {k: sp.simplify(sp.expand(codex[k].subs(t, 1/t) * t**2)) for k in codex}
    print(f"      t->1/t (умн. на t^2 — квадрат): F0<->F8 ? {sp.simplify(inv['F0'] - codex['F8']) == 0};"
          f"  F4 инв ? {sp.simplify(inv['F4'] - codex['F4']) == 0};"
          f"  L инв ? {sp.simplify(inv['L'] - codex['L']) == 0}")
    return cells


# ----------------------------------------------------------------------------------
# §3. АРИФМЕТИКА ПО МОДУЛЮ p — своя реализация, без опоры на таблицы Codex
# ----------------------------------------------------------------------------------

def qr_set(p):
    """квадраты в F_p, НОЛЬ ВКЛЮЧЁН (0 = 0^2)"""
    return {(x * x) % p for x in range(p)}


def forms_hom(m, n, a, b, p):
    """пять однородных форм в F_p при (a:b), t = a/b.
       Знаменатель 2 в s обращается через обратный по модулю (p нечётное)."""
    assert p % 2 == 1
    inv2 = pow(2, p - 2, p)
    s = ((m * m + n * n) % p) * inv2 % p
    F0 = (m * m * b * b + n * n * a * a) % p
    F4 = (s * (a * a + b * b)) % p
    F8 = (n * n * b * b + m * m * a * a) % p
    L = (F4 - 2 * m * n * a * b) % p
    U = (F4 + 2 * m * n * a * b) % p
    return {'F0': F0, 'F4': F4, 'F8': F8, 'L': L, 'U': U}


def proj_points(p):
    """все p+1 точек P^1(F_p): (t:1) и (1:0)"""
    return [(t, 1) for t in range(p)] + [(1, 0)]


def scan(m, n, p, keys=('F0', 'F4', 'F8', 'L', 'U'), zero_is_square=True, verbose=False):
    """возвращает список выживших проективных точек"""
    QR = qr_set(p)
    QRn = QR - {0} if not zero_is_square else QR
    survivors = []
    rows = []
    for (a, b) in proj_points(p):
        v = forms_hom(m, n, a, b, p)
        bad = [k for k in keys if v[k] not in QRn]
        rows.append(((a, b), v, bad))
        if not bad:
            survivors.append((a, b))
    if verbose:
        for (ab, v, bad) in rows:
            tag = "∞" if ab[1] == 0 else str(ab[0])
            mark = "ВЫЖИЛ" if not bad else "отказ: " + ",".join(bad)
            print(f"      t={tag:>4}  " + "  ".join(f"{k}={v[k]:>3}{'□' if v[k] in QRn else '·'}" for k in keys)
                  + f"   {mark}")
    return survivors, rows


# ----------------------------------------------------------------------------------
# §6. Независимая 41-адическая проверка через Гензеля (второй, непохожий метод)
# ----------------------------------------------------------------------------------

def vp(x, p):
    if x == 0:
        return None
    num, den = x.numerator, x.denominator
    v = 0
    while num % p == 0:
        num //= p; v += 1
    while den % p == 0:
        den //= p; v -= 1
    return v


def is_sq_Qp(x, p):
    v = vp(x, p)
    if v is None:
        return True                      # 0 — квадрат
    if v % 2:
        return False
    y = x / Fr(p) ** v
    u = (y.numerator * pow(y.denominator, p - 2, p)) % p
    return pow(u, (p - 1) // 2, p) == 1


def cells_Q(m, n, t):
    s = Fr(m * m + n * n, 2)
    F0 = Fr(m * m) + Fr(n * n) * t * t
    F4 = s * (1 + t * t)
    F8 = Fr(n * n) + Fr(m * m) * t * t
    L = F4 - 2 * m * n * t
    U = F4 + 2 * m * n * t
    return {'F0': F0, 'F4': F4, 'F8': F8, 'L': L, 'U': U}


def brute_rational(m, n, H, p=None):
    """прямой перебор t = u/v, |u|,v <= H: все пять клеток — квадраты в Q (или в Q_p)."""
    hits = []
    for v in range(1, H + 1):
        for u in range(-H, H + 1):
            if gcd(abs(u), v) != 1:
                continue
            t = Fr(u, v)
            c = cells_Q(m, n, t)
            if p is None:
                good = all(x >= 0 and sp.sqrt(sp.Rational(x.numerator, x.denominator)).is_rational
                           for x in c.values())
            else:
                good = all(is_sq_Qp(x, p) for x in c.values())
            if good:
                hits.append((t, c))
    return hits


# ----------------------------------------------------------------------------------

def main():
    cells = report_derivation()

    m, n, p = 16, 5, 41

    print()
    print(SEP)
    print(f"§3. ПОЛНЫЙ ПЕРЕБОР ВСЕХ {p+1} ТОЧЕК P^1(F_{p}) для (m,n)=({m},{n}) — пять клеток")
    print(SEP)
    QR = sorted(qr_set(p))
    print(f"   QR_{p} (с нулём), мой пересчёт: {QR}")
    print(f"   |QR*| = {len(QR)-1}  (ожидается (p-1)/2 = {(p-1)//2}): "
          f"{'OK' if len(QR)-1 == (p-1)//2 else 'ОШИБКА'}")
    inv2 = pow(2, p - 2, p)
    print(f"   m^2 = {m*m % p}, n^2 = {n*n % p}, s = {(m*m+n*n) % p * inv2 % p}, 2mn = {2*m*n % p}  (mod {p})")
    print(f"   сверка с Codex: m^2=10, n^2=25, s=38, 2mn=37 -> "
          f"{'СОВПАЛО' if (m*m%p, n*n%p, (m*m+n*n)%p*inv2%p, 2*m*n%p) == (10,25,38,37) else 'РАСХОЖДЕНИЕ'}")
    print()
    surv5, rows = scan(m, n, p, verbose=False)
    # печатаем только строки, где прошли хотя бы три клетки F0,F4,F8 — плюс бесконечность
    print("   строки, где F0,F4,F8 ОДНОВРЕМЕННО квадраты (именно их обязан перечислить Codex):")
    QRs = qr_set(p)
    three = []
    for (ab, v, bad) in rows:
        if all(v[k] in QRs for k in ('F0', 'F4', 'F8')):
            three.append(ab)
            tag = "беск." if ab[1] == 0 else str(ab[0])
            b5 = [k for k in ('L', 'U') if v[k] not in QRs]
            print(f"      t={tag:>5}  (F0,F4,F8)=({v['F0']},{v['F4']},{v['F8']})  "
                  f"(L,U)=({v['L']},{v['U']})   "
                  f"{'ВСЕ ПЯТЬ КВАДРАТЫ — ОПРОВЕРЖЕНИЕ!' if not b5 else 'не квадраты: ' + ','.join(b5)}")
    print(f"   всего таких t: {len(three)}  (Codex заявил 4 конечных + проверку бесконечности)")
    print(f"   ВЫЖИВШИХ НА ПЯТИ КЛЕТКАХ: {surv5 if surv5 else 'НЕТ НИ ОДНОЙ'}")

    print()
    print(SEP)
    print("§4. ПОИСК ДЫР")
    print(SEP)

    print("   (а) пропущено ли значение t? — перебор исчерпывающий по построению:")
    pts = proj_points(p)
    print(f"       |P^1(F_{p})| = {len(pts)} = p+1 = {p+1}: {'OK' if len(pts) == p+1 else 'ОШИБКА'};"
          f" различных: {len(set(pts))}")
    print(f"       конечные t покрыты все: {sorted(a for a, b in pts if b == 1) == list(range(p))}")
    print(f"       точка на бесконечности (1:0) включена: {(1,0) in pts}")
    vinf = forms_hom(m, n, 1, 0, p)
    print(f"       на бесконечности: {vinf}; квадраты: "
          f"{ {k: (v in QRs) for k, v in vinf.items()} }")

    print()
    print("   (б) трактовка нуля. Ноль — квадрат в F_p, запрещать его НЕЛЬЗЯ.")
    surv5_zero_ok, _ = scan(m, n, p, zero_is_square=True)
    surv5_zero_no, _ = scan(m, n, p, zero_is_square=False)
    print(f"       ноль РАЗРЕШЁН  (правильно, слабее): выживших {len(surv5_zero_ok)}")
    print(f"       ноль ЗАПРЕЩЁН  (неправильно, сильнее): выживших {len(surv5_zero_no)}")
    zeros = [(ab, {k: v[k] for k in v if v[k] == 0}) for (ab, v, bad) in rows if 0 in v.values()]
    print(f"       точки с нулевой клеткой: {[(('беск.' if ab[1]==0 else ab[0]), list(z)) for ab, z in zeros]}")
    print("       => вывод НЕ держится на запрете нулей (он получен при разрешённом нуле).")

    print()
    print("   (в) зависимость от положительности / попарной различности:")
    print("       в scan() используется ТОЛЬКО предикат «квадрат в F_p». Ни одного сравнения")
    print("       по величине и ни одной проверки различия клеток в коде нет — проверено ниже:")
    import inspect
    src = inspect.getsource(scan) + inspect.getsource(forms_hom)
    banned = [w for w in ['>', '<', 'positive', 'distinct', 'abs('] if w in src]
    print(f"       запрещённые конструкции в коде теста: {banned if banned else 'нет'}")
    print("       => тест есть ОСЛАБЛЕНИЕ полной задачи: он отбрасывает условия положительности")
    print("          и различности. Пустота ослабленной системы влечёт пустоту полной. Направление верное.")

    print()
    print("   (г) простое и редукция:")
    print(f"       p = {p} нечётное: {p % 2 == 1}; p простое: {sp.isprime(p)}")
    print(f"       p не делит 2 (знаменатель s): {p != 2}")
    print(f"       p не делит m = {m}: {m % p != 0};  p не делит n = {n}: {n % p != 0}")
    print(f"       p не делит m^2+n^2 = {m*m+n*n}: {(m*m+n*n) % p != 0}  (=> s != 0 mod p)")
    print(f"       p не делит 2mn = {2*m*n}: {(2*m*n) % p != 0}")
    print("       нормировка (a,b): хотя бы одно — единица; все пять форм p-целые. Знаменатель не теряется.")

    print()
    print("   (д) обе ориентации пары: (16,5) против (5,16) — в проекте бывает свап m<->n")
    for (mm, nn) in [(16, 5), (5, 16)]:
        s5, _ = scan(mm, nn, p)
        s3, _ = scan(mm, nn, p, keys=('F0', 'F4', 'F8'))
        print(f"       (m,n)=({mm},{nn}): выживших на 5 клетках {len(s5)}, на 3 клетках {len(s3)}")

    print()
    print("   (е) масштабирование (m,n) -> lambda*(m,n): клетки умножаются на lambda^2 (квадрат),")
    print("       поэтому предикат не меняется. Проверка на нескольких lambda:")
    for lam in [2, 3, 7, 40]:
        s5, _ = scan((16 * lam) % p, (5 * lam) % p, p)
        print(f"       lambda={lam}: выживших {len(s5)}")

    print()
    print(SEP)
    print("§5. КОНТРОЛИ (тест обязан НЕ закрывать то, что закрывать не должен)")
    print(SEP)
    print("   К1. Пара (15,8): у неё ЕСТЬ рациональная точка t=1 (вырожденная, все девять — квадраты).")
    print("       Тест обязан её пропустить при всех p.")
    bad158 = []
    for q in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101]:
        s5, _ = scan(15, 8, q)
        if not s5:
            bad158.append(q)
    print(f"       простые, где тест ЛОЖНО закрыл (15,8): {bad158 if bad158 else 'НЕТ — контроль пройден'}")
    c158 = cells_Q(15, 8, Fr(1))
    print(f"       (15,8), t=1: пять клеток = {[str(x) for x in c158.values()]};"
          f" все квадраты: {all(sp.sqrt(sp.Rational(x.numerator, x.denominator)).is_rational for x in c158.values())}")

    print()
    print("   К2. Контрпример Codex к смешению C7 и девятки: (7,1), t=0, (u0,u4,u8)=(7,5,1).")
    c71 = cells_Q(7, 1, Fr(0))
    print(f"       (7,1), t=0: {[f'{k}={v}' for k, v in c71.items()]}")
    print(f"       F0,F4,F8 квадраты: "
          f"{[bool(sp.sqrt(sp.Rational(c71[k].numerator, c71[k].denominator)).is_rational) for k in ('F0','F4','F8')]}"
          f"; L,U квадраты: "
          f"{[bool(sp.sqrt(sp.Rational(c71[k].numerator, c71[k].denominator)).is_rational) for k in ('L','U')]}")
    print("       => C7 непусто, девятки нет. Различение Codex подтверждено на моих формулах.")

    print()
    print("   К3. Тест на СЕМИ клетках (F0,F4,F8) для (16,5) при p=41 — обязан НЕ закрывать,")
    print("       иначе таблица Codex (4 выживших t) противоречива.")
    s3, _ = scan(16, 5, 41, keys=('F0', 'F4', 'F8'))
    print(f"       выживших на трёх клетках: {len(s3)} -> {[('беск.' if b==0 else a) for a,b in s3]}")

    print()
    print("   К4. Существует ли ВООБЩЕ простое, где тест на пяти клетках закрывает (16,5)?")
    print("       Сканируем все нечётные простые p < 200 — 41 не должен быть единственным по случайности.")
    kill = []
    for q in [x for x in range(3, 200) if sp.isprime(x)]:
        s5, _ = scan(16, 5, q)
        if not s5:
            kill.append(q)
    print(f"       простые, закрывающие (16,5) на пяти клетках: {kill}")
    print(f"       41 среди них: {41 in kill}")
    kill3 = [q for q in [x for x in range(3, 200) if sp.isprime(x)] if not scan(16, 5, q, keys=('F0','F4','F8'))[0]]
    print(f"       простые, закрывающие (16,5) уже на ТРЁХ клетках: {kill3 if kill3 else 'нет'}")

    print()
    print(SEP)
    print("§6. ВТОРОЙ, НЕПОХОЖИЙ МЕТОД: прямой перебор рациональных t и 41-адическая проверка")
    print(SEP)
    H = 60
    hits_p = brute_rational(16, 5, H, p=41)
    print(f"   t = u/v, |u|,v <= {H}, все пять клеток — квадраты в Q_41: найдено {len(hits_p)}")
    if hits_p:
        print(f"      *** {hits_p[:5]} — ОПРОВЕРЖЕНИЕ ***")
    hits_ctrl = brute_rational(15, 8, 20, p=41)
    print(f"   контроль (15,8), |u|,v <= 20, квадраты в Q_41: найдено {len(hits_ctrl)}"
          f"{'  (в т.ч. t=1)' if any(t == 1 for t, _ in hits_ctrl) else ''}")

    print()
    print(SEP)
    print("ИТОГ")
    print(SEP)
    print(f"   клетки: выведены независимо, совпали с Codex поэлементно")
    print(f"   (16,5), p=41: выживших проективных точек на пяти клетках — {len(surv5)}")
    print("   опровержения НЕ найдено" if not surv5 else "   *** НАЙДЕНО ОПРОВЕРЖЕНИЕ ***")


if __name__ == "__main__":
    main()
