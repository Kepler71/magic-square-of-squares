# -*- coding: utf-8 -*-
# ============================================================================
#  Проверка гипотезы «при чётном n красные клетки семейства G1 не квадраты»
#  (зацепка Grok, GROK_IDEAS_2026-09-12.md п.2 — со ссылкой на Boyer/Wesolowski).
#
#  Запуск:  sage /home/kep/magicKube/bridge/even_n_численно.sage
#
#  Разделы:
#   A. Символьный вывод девяти клеток G1; идентификация «красных» клеток.
#   B. Целочисленная модель (t = u/v) и мод-8 леммы.
#   C. ТОЧНОЕ решение вопроса «есть ли точка над Q_2» (BFS по P^1(Z/2^k)) —
#      это доказательство, а не таблица: процедура завершается либо
#      «все классы убиты» (препятствие ДОКАЗАНО), либо предъявляет
#      открытое множество решений (препятствия НЕТ, тоже доказано).
#   D. Табличная разведка: вычеты красных клеток mod 8/16/32/64 и v_2.
#   E. Реальный поиск рациональных t, где обе красные клетки — квадраты.
# ============================================================================

from sage.all import *

Rp = PolynomialRing(QQ, 'p'); pv = Rp.gen()
COEF = [(1,1,0),(1,-1,-1),(1,0,1),(1,-1,1),(1,0,0),(1,1,-1),(1,0,-1),(1,1,1),(1,-1,0)]

def cells_G1(m, n):
    """Девять клеток семейства G1 как многочлены от параметра p = t."""
    m, n = QQ(m), QQ(n)
    P = m*pv + n; Q = m - n*pv; R = m*pv - n; S = m + n*pv
    c = (m^2 + n^2)*(pv^2 + 1)/2
    M = matrix(QQ, [[COEF[1][1], COEF[1][2]], [COEF[3][1], COEF[3][2]]]); Mi = M.inverse()
    ab = Mi*vector(Rp, [Rp(P^2 - c), Rp(R^2 - c)])
    a, b = Rp(ab[0]), Rp(ab[1])
    return [Rp(cc*c + ca*a + cb*b) for cc, ca, cb in COEF]

# ---------------------------------------------------------------- A ---------
print("="*78)
print("A. СИМВОЛЬНАЯ ИДЕНТИФИКАЦИЯ КЛЕТОК G1")
print("="*78)
Rmn = PolynomialRing(QQ, ['M','N','T']); (Mv, Nv, Tv) = Rmn.gens()
Fr = Rmn.fraction_field()
Pp = Mv*Tv + Nv; Qq = Mv - Nv*Tv; Rr = Mv*Tv - Nv; Ss = Mv + Nv*Tv
cc = (Mv^2 + Nv^2)*(Tv^2 + 1)/2
bb = (Rr^2 - Pp^2)/2
aa = cc - (Pp^2 + Rr^2)/2
C = [cc + aa, cc - aa - bb, cc + bb, cc - aa + bb, cc, cc + aa - bb, cc - bb, cc + aa + bb, cc - aa]
names = ["c0","c1","c2","c3","c4","c5","c6","c7","c8"]
for i, x in enumerate(C):
    print(f"  {names[i]} = {Fr(x).numerator()} / {Fr(x).denominator()}")
rows = [C[0]+C[1]+C[2], C[3]+C[4]+C[5], C[6]+C[7]+C[8]]
cols = [C[0]+C[3]+C[6], C[1]+C[4]+C[7], C[2]+C[5]+C[8]]
dgs  = [C[0]+C[4]+C[8], C[2]+C[4]+C[6]]
allsums = rows+cols+dgs
print("  магический (все 8 сумм равны 3*центр):", all(x == 3*C[4] for x in allsums))
print("  c1 = (mt+n)^2 :", C[1] == Pp^2, "  c3 = (mt-n)^2 :", C[3] == Rr^2,
      "  c5 = (m+nt)^2 :", C[5] == Ss^2, "  c7 = (m-nt)^2 :", C[7] == Qq^2)
print("  c0 = m^2+n^2 t^2 = F0 :", C[0] == Mv^2 + Nv^2*Tv^2)
print("  c4 = s(1+t^2)    = F4 :", C[4] == (Mv^2+Nv^2)*(1+Tv^2)/2)
print("  c8 = n^2+m^2 t^2 = F8 :", C[8] == Nv^2 + Mv^2*Tv^2)
print("  => СВОБОДНЫЕ (не автоматические) клетки: c0,c4,c8 (кривая C рода 5) и c2,c6.")
print("  c2 = s(1+t^2) - 2mnt :", C[2] == (Mv^2+Nv^2)*(1+Tv^2)/2 - 2*Mv*Nv*Tv)
print("  c6 = s(1+t^2) + 2mnt :", C[6] == (Mv^2+Nv^2)*(1+Tv^2)/2 + 2*Mv*Nv*Tv)
print("  тождества: c2+c6 = 2c4 :", C[2]+C[6] == 2*C[4], "; c6-c2 = 4mnt :", C[6]-C[2] == 4*Mv*Nv*Tv)
print("  c2 = (c3+c7)/2 :", C[2] == (C[3]+C[7])/2, " ; c6 = (c1+c5)/2 :", C[6] == (C[1]+C[5])/2)
print("  при t=1: c2 = (m-n)^2, c6 = (m+n)^2, c0=c4=c8=m^2+n^2 :",
      C[2].subs({Tv:1}) == (Mv-Nv)^2, C[6].subs({Tv:1}) == (Mv+Nv)^2,
      C[0].subs({Tv:1}) == Mv^2+Nv^2, C[4].subs({Tv:1}) == Mv^2+Nv^2, C[8].subs({Tv:1}) == Mv^2+Nv^2)

# ---------------------------------------------------------------- B ---------
print()
print("="*78)
print("B. ЦЕЛОЧИСЛЕННАЯ МОДЕЛЬ  t = u/v,  умножение всех клеток на 4v^2 (квадрат)")
print("="*78)
print("""  H0 = 4(m^2 v^2 + n^2 u^2)                         [=4v^2 c0]
  H4 = 2(m^2+n^2)(u^2+v^2)                          [=4v^2 c4]
  H8 = 4(n^2 v^2 + m^2 u^2)                         [=4v^2 c8]
  H2 = 2[(mu-nv)^2 + (mv-nu)^2]                     [=4v^2 c2]   (красная)
  H6 = 2[(mu+nv)^2 + (mv+nu)^2]                     [=4v^2 c6]   (красная)
  c_i квадрат в Q  <=>  H_i квадрат в Q;  то же покомпонентно над Q_2.""")

def H_forms(m, n):
    m, n = ZZ(m), ZZ(n)
    return {
        'H0': lambda u, v: 4*(m*m*v*v + n*n*u*u),
        'H4': lambda u, v: 2*(m*m+n*n)*(u*u+v*v),
        'H8': lambda u, v: 4*(n*n*v*v + m*m*u*u),
        'H2': lambda u, v: 2*((m*u-n*v)**2 + (m*v-n*u)**2),
        'H6': lambda u, v: 2*((m*u+n*v)**2 + (m*v+n*u)**2),
    }

# проверка тождеств H2, H6 против прямой подстановки
_t = Rmn.gens()[2]
for (m0, n0) in [(13,8),(15,8),(19,16),(16,5),(3,2),(7,1),(5,3)]:
    Hs = H_forms(m0, n0); ok = True
    for u0 in range(-4, 5):
        for v0 in range(1, 6):
            if gcd(u0, v0) != 1: continue
            tt = QQ(u0)/QQ(v0)
            vals = [QQ(m0^2+n0^2)*(1+tt^2)/2 - 2*m0*n0*tt, QQ(m0^2+n0^2)*(1+tt^2)/2 + 2*m0*n0*tt,
                    m0^2+n0^2*tt^2, QQ(m0^2+n0^2)*(1+tt^2)/2, n0^2+m0^2*tt^2]
            hs = [Hs['H2'](u0,v0), Hs['H6'](u0,v0), Hs['H0'](u0,v0), Hs['H4'](u0,v0), Hs['H8'](u0,v0)]
            for a_, b_ in zip(vals, hs):
                if 4*v0^2*a_ != b_: ok = False
    print(f"  тождества H_i = 4v^2 c_i для (m,n)=({m0},{n0}): {ok}")

# ---------------------------------------------------------------- C ---------
print()
print("="*78)
print("C. ТОЧНОЕ РЕШЕНИЕ: есть ли точка над Q_2  (BFS по P^1(Z/2^k), доказательство)")
print("="*78)
print("""  Классы (u:v) mod 2^k, (u,v) примитивны. Для формы H и класса:
    r = H(u0,v0) mod 2^k; всякий лифт даёт H ≡ r (mod 2^k).
    Если r != 0 и e = v_2(r) <= k-3, то v_2(H)=e и единичная часть определена mod 8:
      квадрат <=> e чётно и (r >> e) ≡ 1 (mod 8)   — решается ТОЧНО.
    Иначе класс дробится на 4 подкласса mod 2^{k+1}.
  Ни одна из форм не обращается в 0 на P^1(Q_2) (нужно было бы -1 ∈ (Q_2*)^2),
  поэтому процедура конечна: каждый класс либо убивается, либо становится решением.""")

def two_adic_decide(m, n, keys, kmax=28, want_witness=True):
    """Возвращает ('SOLUTION', (u,v,k)) или ('NO_POINT', kfinal) или ('UNDECIDED', k)."""
    Hs = H_forms(m, n)
    forms = [Hs[k_] for k_ in keys]
    cur = [(1,0,1), (0,1,1), (1,1,1)]      # (u0,v0,k) mod 2^k = 2
    k = 1
    while cur and k <= kmax:
        nxt = []
        for (u0, v0, kk) in cur:
            mod = 1 << kk
            verdict = 'SQUARE_ALL'
            for f in forms:
                r = f(u0, v0) % mod
                if r == 0:
                    verdict = 'REFINE'; break
                e = ZZ(r).valuation(2)
                if e > kk - 3:
                    verdict = 'REFINE'; break
                if e % 2 == 1 or ((r >> e) % 8) != 1:
                    verdict = 'DEAD'; break
            if verdict == 'DEAD':
                continue
            if verdict == 'SQUARE_ALL':
                return ('SOLUTION', (u0, v0, kk))
            for a in range(2):
                for b in range(2):
                    nxt.append((u0 + a*mod, v0 + b*mod, kk+1))
        cur = nxt
        k += 1
    if not cur:
        return ('NO_POINT', k)
    return ('UNDECIDED', k)

PAIRS_EVEN = [(13,8),(15,8),(19,16),(16,5),(3,2),(5,2),(7,2),(9,2),(11,2),(7,4),(9,4),
              (5,4),(11,8),(17,8),(21,16),(9,8),(13,4),(15,4),(3,4),(1,4),(1,8),(5,8)]
PAIRS_ODD  = [(3,1),(5,1),(5,3),(7,1),(7,3),(7,5),(9,1),(11,1),(15,1),(19,5),(13,1),(9,5),(11,3)]

def cls(m, n):
    if (m % 2 == 1) and (n % 2 == 1): return "оба нечётны"
    e = n if n % 2 == 0 else m
    return f"чётный параметр ≡ {e % 8} mod 8 (≡ {e % 4} mod 4)"

print()
print("  СИСТЕМА S3 = {c0,c4,c8}  — это кривая C рода 5 (СЕМЬ квадратов):")
print(f"  {'(m,n)':>10} {'класс':>34}   вердикт над Q_2")
res3 = {}
for (m0,n0) in PAIRS_EVEN + PAIRS_ODD:
    if gcd(m0,n0) != 1: continue
    st = two_adic_decide(m0, n0, ['H0','H4','H8'])
    res3[(m0,n0)] = st
    tag = {"SOLUTION":"ТОЧКА ЕСТЬ (локально разрешима)","NO_POINT":"ТОЧЕК НЕТ — 2-адическое препятствие","UNDECIDED":"не решено"}[st[0]]
    extra = f"   свидетель (u,v)≡({st[1][0]},{st[1][1]}) mod 2^{st[1][2]}" if st[0]=="SOLUTION" else f"   (убито к k={st[1]})"
    print(f"  {str((m0,n0)):>10} {cls(m0,n0):>34}   {tag}{extra}")

print()
print("  СИСТЕМА S5 = {c0,c4,c8,c2,c6} — ВСЕ ДЕВЯТЬ клеток квадраты:")
print(f"  {'(m,n)':>10} {'класс':>34}   вердикт над Q_2")
res5 = {}
for (m0,n0) in PAIRS_EVEN + PAIRS_ODD:
    if gcd(m0,n0) != 1: continue
    st = two_adic_decide(m0, n0, ['H0','H4','H8','H2','H6'])
    res5[(m0,n0)] = st
    tag = {"SOLUTION":"ТОЧКА ЕСТЬ (локально разрешима)","NO_POINT":"ТОЧЕК НЕТ — 2-адическое препятствие","UNDECIDED":"не решено"}[st[0]]
    extra = f"   свидетель (u,v)≡({st[1][0]},{st[1][1]}) mod 2^{st[1][2]}" if st[0]=="SOLUTION" else f"   (убито к k={st[1]})"
    print(f"  {str((m0,n0)):>10} {cls(m0,n0):>34}   {tag}{extra}")

print()
print("  СИСТЕМА S2 = {c2,c6} — ТОЛЬКО красные клетки (гипотеза Grok в чистом виде):")
for (m0,n0) in PAIRS_EVEN + PAIRS_ODD:
    if gcd(m0,n0) != 1: continue
    st = two_adic_decide(m0, n0, ['H2','H6'])
    tag = {"SOLUTION":"ТОЧКА ЕСТЬ","NO_POINT":"ТОЧЕК НЕТ","UNDECIDED":"не решено"}[st[0]]
    print(f"  {str((m0,n0)):>10} {cls(m0,n0):>34}   {tag}")

# --- полная развёртка по классам (m mod 8, n mod 8) ------------------------
print()
print("  ПОЛНАЯ РАЗВЁРТКА: вердикт зависит только от (m mod 8, n mod 8)? — проверка")
print("  (2-адическая разрешимость S5 для всех допустимых пар m,n <= 40, gcd=1)")
byclass = {}
for m0 in range(1, 41):
    for n0 in range(1, 41):
        if gcd(m0,n0) != 1: continue
        st = two_adic_decide(m0, n0, ['H0','H4','H8','H2','H6'])
        key = (m0 % 8, n0 % 8)
        byclass.setdefault(key, set()).add(st[0])
bad = {k_: v_ for k_, v_ in byclass.items() if len(v_) > 1}
print("   классы (m mod 8, n mod 8) с неоднородным вердиктом:", bad if bad else "нет — вердикт зависит только от m,n mod 8")
print("   таблица вердиктов S5 по (m mod 8, n mod 8):")
print("      m\\n " + " ".join(f"{j:>5}" for j in range(8)))
for i in range(8):
    row = []
    for j in range(8):
        v_ = byclass.get((i,j))
        row.append("-" if v_ is None else ("ЕСТЬ" if 'SOLUTION' in v_ else ("НЕТ" if 'NO_POINT' in v_ else "?")))
    print(f"      {i:>3} " + " ".join(f"{x:>5}" for x in row))

# ---------------------------------------------------------------- D ---------
print()
print("="*78)
print("D. ТАБЛИЧНАЯ РАЗВЕДКА: вычеты красных клеток mod 8/16/32/64 и v_2")
print("="*78)

def scan_residues(m, n, B=160, cond_c048=False):
    """Собирает статистику по t=u/v, 1<=v<=B, |u|<=B, gcd(u,v)=1.
       cond_c048: учитывать только те (u,v), где c0,c4,c8 — КВАДРАТЫ В Q_2."""
    Hs = H_forms(m, n)
    stat = {'n':0, 'c2_sq2':0, 'c6_sq2':0, 'both_sq2':0,
            'r2_8':{}, 'r6_8':{}, 'r2_16':{}, 'v2_c2':{}, 'both_sq_Q':0, 'both_sq_Q_nondeg':[]}
    for v in range(1, B+1):
        for u in range(-B, B+1):
            if gcd(u, v) != 1: continue
            if u == 0: continue
            h0, h4, h8 = Hs['H0'](u,v), Hs['H4'](u,v), Hs['H8'](u,v)
            if cond_c048:
                ok2 = True
                for h in (h0,h4,h8):
                    e = ZZ(h).valuation(2)
                    if e % 2 or ((h >> e) % 8) != 1: ok2 = False; break
                if not ok2: continue
            h2, h6 = Hs['H2'](u,v), Hs['H6'](u,v)
            stat['n'] += 1
            def sq2(h):
                e = ZZ(h).valuation(2)
                return (e % 2 == 0) and ((h >> e) % 8 == 1)
            s2, s6 = sq2(h2), sq2(h6)
            stat['c2_sq2'] += int(s2); stat['c6_sq2'] += int(s6); stat['both_sq2'] += int(s2 and s6)
            e2 = ZZ(h2).valuation(2)
            stat['v2_c2'][e2] = stat['v2_c2'].get(e2,0)+1
            odd2 = h2 >> e2; odd6 = h6 >> ZZ(h6).valuation(2)
            stat['r2_8'][odd2 % 8] = stat['r2_8'].get(odd2 % 8,0)+1
            stat['r6_8'][odd6 % 8] = stat['r6_8'].get(odd6 % 8,0)+1
            stat['r2_16'][odd2 % 16] = stat['r2_16'].get(odd2 % 16,0)+1
            if ZZ(h2).is_square() and ZZ(h6).is_square():
                stat['both_sq_Q'] += 1
                if abs(u) != v and u != 0:
                    stat['both_sq_Q_nondeg'].append((u,v))
    return stat

print("  Обозначения: sq2 = квадрат в Q_2; odd(x) = x / 2^{v_2(x)}.")
print("  «c2,c6 обе квадраты в Q» — это уже семейство Весоловского (сужение G1).")
for tag, lst in [("ЧЁТНЫЙ параметр", [(13,8),(15,8),(19,16),(16,5),(7,4),(9,4),(3,2),(5,2),(7,2),(9,2),(11,2),(5,4)]),
                 ("ОБА НЕЧЁТНЫ (контроль)", [(3,1),(5,1),(5,3),(7,1),(15,1),(19,5)])]:
    print(f"\n  --- {tag} ---")
    print(f"  {'(m,n)':>9} {'всего':>7} {'c2 sq2':>7} {'c6 sq2':>7} {'обе sq2':>8}  HIST_HDR  {'v_2(H2)':>16}  обе кв. в Q")
    for (m0,n0) in lst:
        if gcd(m0,n0) != 1: continue
        st = scan_residues(m0, n0, B=120)
        hist = ",".join(f"{k_}:{v_}" for k_, v_ in sorted(st['r2_8'].items()))
        vh = ",".join(f"{k_}:{v_}" for k_, v_ in sorted(st['v2_c2'].items()))
        nd = st['both_sq_Q_nondeg'][:4]
        print(f"  {str((m0,n0)):>9} {st['n']:>7} {st['c2_sq2']:>7} {st['c6_sq2']:>7} {st['both_sq2']:>8}  {hist:>34}  {vh:>16}  {st['both_sq_Q']} {nd}")

print()
print("  То же, но ТОЛЬКО по (u,v), где c0,c4,c8 — квадраты в Q_2 (т.е. на 2-адической кривой C):")
print(f"  {'(m,n)':>9} {'подходящих':>11} {'c2 sq2':>7} {'c6 sq2':>7} {'обе sq2':>8}   odd(c2) mod 8")
for (m0,n0) in [(13,8),(15,8),(19,16),(16,5),(3,2),(5,2),(7,4),(9,4),(3,1),(5,1),(5,3),(7,1),(15,1),(19,5)]:
    if gcd(m0,n0) != 1: continue
    st = scan_residues(m0, n0, B=120, cond_c048=True)
    hist = ",".join(f"{k_}:{v_}" for k_, v_ in sorted(st['r2_8'].items()))
    print(f"  {str((m0,n0)):>9} {st['n']:>11} {st['c2_sq2']:>7} {st['c6_sq2']:>7} {st['both_sq2']:>8}   {hist}")

# ---------------------------------------------------------------- E ---------
print()
print("="*78)
print("E. РЕАЛЬНЫЕ рациональные t, где ОБЕ красные клетки — квадраты (семейство Весоловского)")
print("="*78)
print("  Условие c2=□, c6=□ эквивалентно (c2+c6=2c4, c6-c2=4mnt):")
print("  c2=(alpha)^2, c6=(beta)^2 => alpha^2+beta^2 = 2c4, beta^2-alpha^2 = 4mnt.")
print("  Вырожденные решения t=±1 есть ВСЕГДА: c2=(m∓n)^2, c6=(m±n)^2.")
print()
print(f"  {'(m,n)':>9}  невырожденные (u,v) с H2,H6 квадратами, |u|,v<=400")
for (m0,n0) in [(13,8),(15,8),(19,16),(16,5),(3,2),(5,2),(7,4),(9,4),(3,1),(5,1),(5,3),(7,1),(15,1),(19,5)]:
    if gcd(m0,n0) != 1: continue
    Hs = H_forms(m0,n0); found = []
    for v in range(1, 401):
        for u in range(-400, 401):
            if u == 0 or abs(u) == v: continue
            if gcd(u,v) != 1: continue
            if ZZ(Hs['H2'](u,v)).is_square() and ZZ(Hs['H6'](u,v)).is_square():
                found.append((u,v))
                if len(found) >= 6: break
        if len(found) >= 6: break
    print(f"  {str((m0,n0)):>9}  {found if found else 'не найдено в этом ящике (НЕ доказательство отсутствия)'}")

print()
print("="*78)
print("G. ПОЛНАЯ 2-АДИЧЕСКАЯ КЛАССИФИКАЦИЯ: добавляют ли красные клетки хоть что-то при p=2?")
print("="*78)
print("""  Процедура: BFS дробит P^1(Z_2) до классов, где ВСЕ ПЯТЬ форм решены точно.
  Классы, где c0,c4,c8 — квадраты (2-адическая часть кривой C), собираются;
  на каждом смотрим вердикт для c2 и c6. Если всюду «квадрат» — красные клетки
  при p=2 не дают НИЧЕГО нового. Это доказательство (конечный полный разбор).""")

def two_adic_classify(m, n, kmax=30):
    """Полный разбор. Возвращает (n_class_C, n_class_C_red_ok, примеры-исключения)."""
    Hs = H_forms(m, n)
    order = ['H0','H4','H8','H2','H6']
    cur = [(1,0,1), (0,1,1), (1,1,1)]
    good, bad = 0, []
    while cur:
        nxt = []
        for (u0, v0, kk) in cur:
            if kk > kmax:
                bad.append(('UNRESOLVED', u0, v0, kk)); continue
            mod = 1 << kk
            dec = {}
            refine = False
            for key in order:
                r = Hs[key](u0, v0) % mod
                if r == 0: refine = True; break
                e = ZZ(r).valuation(2)
                if e > kk - 3: refine = True; break
                dec[key] = (e % 2 == 0) and (((r >> e) % 8) == 1)
                if key in ('H0','H4','H8') and not dec[key]:
                    dec = None; break
            if dec is None:   # класс не лежит на 2-адической кривой C — не интересует
                continue
            if refine:
                for a in range(2):
                    for b in range(2):
                        nxt.append((u0 + a*mod, v0 + b*mod, kk+1))
                continue
            # класс целиком лежит на 2-адической C; смотрим красные
            if dec['H2'] and dec['H6']:
                good += 1
            else:
                bad.append((u0, v0, kk, dec['H2'], dec['H6']))
        cur = nxt
    return good, bad

print()
print(f"  {'(m,n)':>9}  классов на 2-адич. C  красные тоже квадраты?  исключения")
for (m0,n0) in [(13,8),(15,8),(19,16),(16,5),(7,4),(9,4),(5,4),(11,8),(3,4),
                (5,3),(7,1),(15,1),(19,5),(9,1),(11,3),(3,2),(5,2),(3,1),(5,1),(9,5)]:
    if gcd(m0,n0) != 1: continue
    g_, b_ = two_adic_classify(m0, n0)
    verdict = "ДА (всегда)" if (g_ > 0 and not b_) else ("C пуста над Q_2" if g_ == 0 and not b_ else "НЕТ")
    print(f"  {str((m0,n0)):>9}  {g_:>19}  {verdict:>22}  {b_[:2] if b_ else '—'}")

print()
print("  То же для ВСЕХ взаимно простых m,n <= 30 (сводка):")
tot_ok = tot_empty = tot_bad = 0
for m0 in range(1, 31):
    for n0 in range(1, 31):
        if gcd(m0,n0) != 1: continue
        g_, b_ = two_adic_classify(m0, n0)
        if b_: tot_bad += 1
        elif g_ == 0: tot_empty += 1
        else: tot_ok += 1
print(f"   пар, где C пуста над Q_2: {tot_empty};  пар, где красные автоматически квадраты в Q_2: {tot_ok};"
      f"  пар с исключением: {tot_bad}")

print()
print("="*78)
print("H. ЯВНЫЕ РАЦИОНАЛЬНЫЕ ТОЧКИ, ГДЕ ОБЕ КРАСНЫЕ КЛЕТКИ — ТОЧНЫЕ КВАДРАТЫ (чётный параметр!)")
print("="*78)
for (m0, n0, u0, v0) in [(13,8,27,11), (13,8,11,27), (15,8,289,41), (15,1,41,11)]:
    C9 = cells_G1(m0, n0); tt = QQ(u0)/QQ(v0)
    vals = [c(tt) for c in C9]
    sq = [(x*v0^2*1).is_square() for x in vals]
    print(f"  (m,n)=({m0},{n0}), t={u0}/{v0}:")
    print(f"     9 клеток (×{v0}^2): " + ", ".join(str(x*v0^2) for x in vals))
    print(f"     квадраты:           " + ", ".join(("ДА " if s else "нет") for s in sq))
    print(f"     красные c2,c6 = {vals[2]*v0^2} = {sqrt(vals[2]*v0^2) if (vals[2]*v0^2).is_square() else '?'}^2, "
          f"{vals[6]*v0^2} = {sqrt(vals[6]*v0^2) if (vals[6]*v0^2).is_square() else '?'}^2")

print()
print("="*78)
print("I. КРИВАЯ «ОБЕ КРАСНЫЕ КЛЕТКИ — КВАДРАТЫ»: род 1, ранг, бесконечность решений")
print("="*78)
print("""  c2 = s t^2 - 2mn t + s = □  — коника с рациональной точкой t=1 (c2=(m-n)^2);
  параметризуем её и подставляем в c6 = s t^2 + 2mn t + s = □ => квартика (род 1
  с рациональной точкой) => эллиптическая кривая. rank > 0 == решений бесконечно много.""")
Rk = PolynomialRing(QQ, 'k'); kk_ = Rk.gen()

def red_quartic(m0, n0):
    """y^2 = quart(k): кривая «c2 и c6 одновременно квадраты» после параметризации коники c2=□.
       Прямая через рациональную точку t=1 (c2=(m-n)^2) с наклоном k даёт
       t(k) = N(k)/D(k);  подстановка в c6 даёт quart(k) = c6(t(k))*D(k)^2."""
    s0 = QQ(m0**2+n0**2)/2
    N = -kk_**2 + 2*(m0-n0)*kk_ + (4*m0*n0 - m0**2 - n0**2)/2
    D = s0 - kk_**2
    return s0*N**2 + 2*m0*n0*N*D + s0*D**2, N, D

def jac_of_quartic(q):
    a,b,c,d,e = [q[i] for i in (4,3,2,1,0)]
    I = 12*a*e - 3*b*d + c**2
    J = 72*a*c*e + 9*b*c*d - 27*a*d**2 - 27*e*b**2 - 2*c**3
    return EllipticCurve(QQ, [-27*I, -27*J])

# контроль параметризации: t(k) действительно даёт c2 = квадрат
for (m0,n0) in [(13,8),(15,8),(16,5),(7,1)]:
    q, N, D = red_quartic(m0, n0); s0 = QQ(m0**2+n0**2)/2; ok = True
    for kv in [QQ(1)/3, QQ(2), QQ(-5)/7, QQ(11)/4]:
        if D(kv) == 0: continue
        tv = N(kv)/D(kv)
        if not (s0*tv**2 - 2*m0*n0*tv + s0).is_square(): ok = False
        if (s0*tv**2 + 2*m0*n0*tv + s0)*D(kv)**2 != q(kv): ok = False
    print(f"  контроль параметризации ({m0},{n0}): c2 всегда квадрат и quart = c6*D^2 — {ok}")

print()
print(f"  {'(m,n)':>9}  {'rank':>7}  кручение  quart((m-n)/2)=кв.?   E (min. model)")
for (m0,n0) in [(13,8),(15,8),(19,16),(16,5),(7,4),(9,4),(3,2),(5,2),(7,2),(11,2),(7,1),(15,1),(19,5),(5,3),(3,1),(5,1)]:
    if gcd(m0,n0) != 1: continue
    q, N, D = red_quartic(m0, n0)
    chk = q(QQ(m0-n0)/2).is_square()
    Em = jac_of_quartic(q).minimal_model()
    r = pari(Em.ainvs()).ellinit().ellrank()
    print(f"  {str((m0,n0)):>9}  [{r[0]},{r[1]}]{'':>3}  {Em.torsion_order():>7}  {str(chk):>16}   {Em.ainvs()}")

print()
print("  Прямой поиск k ∈ Q малой высоты с quart(k) = □  (=> ОБЕ красные клетки квадраты):")
for (m0,n0) in [(13,8),(15,8),(19,16),(16,5),(7,4),(9,4),(3,2),(5,2),(7,2),(11,2),(7,1),(15,1),(19,5),(5,3)]:
    if gcd(m0,n0) != 1: continue
    q, N, D = red_quartic(m0, n0); sols = []
    for den in range(1, 61):
        for nu in range(-900, 901):
            if gcd(nu, den) != 1: continue
            kv = QQ(nu)/QQ(den)
            if D(kv) == 0: continue
            val = q(kv)
            if val > 0 and val.is_square():
                tv = N(kv)/D(kv)
                if tv in (QQ(1), QQ(-1), QQ(0)): continue
                if tv not in [x[1] for x in sols]: sols.append((kv, tv))
                if len(sols) >= 3: break
        if len(sols) >= 3: break
    print(f"  {str((m0,n0)):>9}  невырожденные t: {[str(x[1]) for x in sols] if sols else 'нет в этом ящике (НЕ доказательство отсутствия)'}")


print()
print("="*78)
print("J. ПАРЫ С rank = 0: ПОЛНОЕ ПЕРЕЧИСЛЕНИЕ t, ГДЕ ОБЕ КРАСНЫЕ КЛЕТКИ — КВАДРАТЫ")
print("="*78)
print("""  Квартика Q: y^2 = quart(k) имеет рациональную точку => Q ≅ Jac(Q) = E над Q,
  значит #Q(Q) = #E(Q). При rank E = 0 это ровно #E(Q)_tors точек — конечное число.
  Старший коэффициент quart = (m+n)^2 — квадрат, значит 2 точки на бесконечности (они дают t=1).
  Ищем все аффинные k; если найдено ровно #E(Q)_tors - 2 аффинных точек, перечисление ПОЛНОЕ.""")

def count_mod_p(q, p):
    Fp = GF(p); cnt = 0
    for x in Fp:
        y2 = Fp(q(QQ(x.lift())))
        if y2 == 0: cnt += 1
        elif y2.is_square(): cnt += 2
    lead = Fp(q[4])
    if lead != 0 and lead.is_square(): cnt += 2
    elif lead == 0: cnt += 1
    return cnt

for (m0,n0) in [(19,16),(7,4),(3,2),(5,2),(19,5),(5,3),(7,1),(3,1),(5,1)]:
    if gcd(m0,n0) != 1: continue
    q, N, D = red_quartic(m0, n0)
    Em = jac_of_quartic(q).minimal_model()
    r = pari(Em.ainvs()).ellinit().ellrank()
    if r[1] != 0: continue
    tors = Em.torsion_order()
    # независимый контроль изоморфизма Q ≅ E: #Q(F_p) = #E(F_p) на нескольких p
    ctrl = []
    for pp in [11,13,17,19,23,29,31,37]:
        if ZZ(q.denominator()).gcd(pp) != 1: continue
        try:
            if Em.has_good_reduction(pp) and ZZ(q.discriminant()).gcd(pp) == 1:
                ctrl.append(count_mod_p(q, pp) == Em.change_ring(GF(pp)).cardinality())
        except Exception:
            pass
    pts = []
    for den in range(1, 400):
        for nu in range(-4000, 4001):
            if gcd(nu, den) != 1: continue
            kv = QQ(nu)/QQ(den)
            val = q(kv)
            if val == 0: pts.append((kv, 'y=0'))
            elif val > 0 and val.is_square(): pts.append((kv, '±y'))
    naff = sum(2 if x[1]=='±y' else 1 for x in pts)
    ts = sorted(set([N(kv)/D(kv) for kv,_ in pts if D(kv) != 0] + [QQ(1)]))
    inf_t = any(D(kv) == 0 for kv,_ in pts)
    complete = (naff + 2 == tors)
    print(f"  ({m0},{n0}): #E(Q)=tors={tors}; аффинных точек найдено {naff} (+2 на беск.) => "
          f"{'ПЕРЕЧИСЛЕНИЕ ПОЛНОЕ' if complete else 'НЕПОЛНО (найдено меньше)'}"
          f"; контроль #Q(F_p)=#E(F_p): {all(ctrl) if ctrl else 'н/д'} ({len(ctrl)} простых)")
    print(f"          все рациональные t с c2=□ и c6=□ : {[str(x) for x in ts]}"
          f"{', t=∞' if inf_t else ''}   (вырожденные: t ∈ {{0,±1,∞}})")


print()
print("="*78)
print("K. ВСЕ ВЗАИМНО ПРОСТЫЕ ПАРЫ m > n, m <= 20: ранг кривой «обе красные клетки квадраты»")
print("="*78)
print("""  rank = 0  =>  рациональных t с c2=□ и c6=□ только конечное число; ниже проверяется,
  что все они вырожденные (t ∈ {0,±1,∞}) => в этом G1-семействе НЕТ квадрата из ДЕВЯТИ квадратов.
  ВНИМАНИЕ: к задаче о СЕМИ квадратах (кривая C рода 5) это отношения не имеет.""")
rank0, rankpos, unknown = [], [], []
for m0 in range(2, 21):
    for n0 in range(1, m0):
        if gcd(m0,n0) != 1: continue
        q, N, D = red_quartic(m0, n0)
        Em = jac_of_quartic(q).minimal_model()
        r = pari(Em.ainvs()).ellinit().ellrank()
        if r[1] == 0: rank0.append((m0,n0))
        elif r[0] >= 1: rankpos.append((m0,n0,r[0],r[1]))
        else: unknown.append((m0,n0,r[0],r[1]))
print(f"  пар всего: {len(rank0)+len(rankpos)+len(unknown)}")
print(f"  rank = 0 (конечно решений)  : {len(rank0)}  -> {rank0}")
print(f"  rank >= 1 (бесконечно много): {len(rankpos)}  -> {[(a,b) for a,b,_,_ in rankpos]}")
print(f"  ранг не определён           : {len(unknown)}  -> {unknown}")
print()
print("  Разбивка rank=0 по чётности:")
ev = [(a,b) for a,b in rank0 if (a+b) % 2 == 1]
od = [(a,b) for a,b in rank0 if (a+b) % 2 == 0]
print(f"    ровно один из m,n чётен: {len(ev)} -> {ev}")
print(f"    оба нечётны            : {len(od)} -> {od}")
evp = [(a,b) for a,b,_,_ in rankpos if (a+b) % 2 == 1]
odp = [(a,b) for a,b,_,_ in rankpos if (a+b) % 2 == 0]
print("  Разбивка rank>=1 по чётности:")
print(f"    ровно один из m,n чётен: {len(evp)} -> {evp}")
print(f"    оба нечётны            : {len(odp)} -> {odp}")
print("  => если бы «чётность» давала запрет, в строке rank>=1 не было бы ни одной чётной пары.")


print()
print("="*78)
print("L. РЕКОНСТРУКЦИЯ СЕМЕЙСТВА ВЕСОЛОВСКОГО И ЕГО «ЧЁТНОГО n»")
print("="*78)
print("""  Условие магичности квадрата Весоловского (шесть клеток (xy∓z)^2,(xz±y)^2,(x+yz)^2,(yz-x)^2):
     c2 = (xy-z)^2 и c6 = (xy+z)^2  <=>  (x^2+y^2)(1+z^2) = 2(x^2 y^2 + z^2)
     <=>  z^2 = (2x^2y^2 - x^2 - y^2)/(x^2 + y^2 - 2).
  Тогда НЕквадратными остаются ровно центр c4 и две свободные клетки c0, c8 —
  то есть в точности наша кривая C рода 5 (задача о СЕМИ квадратах).
  Прямой поиск решений даёт цепочку x(n+1) = 4x(n) - x(n-1), y(n) = x(n-1), z — та же рекуррента.""")
X = [3, 11]; Z = [1, 4]
while X[-1] < 10**6:
    X.append(4*X[-1] - X[-2]); Z.append(4*Z[-1] - Z[-2])
print()
print(f"  {'n':>3} {'(x,y)':>18} {'z':>10} S_HDR   вердикт C над Q_2   (проверка условия)")
for i in range(1, len(X)-1):
    x_, y_, z_ = X[i], X[i-1], Z[i]
    ok_cond = ((x_^2+y_^2)*(1+z_^2) == 2*(x_^2*y_^2 + z_^2))
    s_ = QQ(x_^2+y_^2)/2
    st = two_adic_decide(x_, y_, ['H0','H4','H8'])
    par = "ЧЁТНОЕ n" if i % 2 == 0 else "нечётное n"
    print(f"  {i:>3} {str((x_,y_)):>18} {z_:>10} {ZZ(s_) % 8:>21}   "
          f"{'ТОЧКИ ЕСТЬ':>12} " if st[0]=='SOLUTION' else f"  {i:>3} {str((x_,y_)):>18} {z_:>10} {ZZ(s_) % 8:>21}   {'C(Q_2) = ПУСТО':>12} ", end="")
    print(f"  [{par}]  условие Весоловского: {ok_cond}")
print()
print("  ВЫВОД: в цепочке Весоловского ЧЁТНЫЕ индексы n — ровно те, где s ≡ 5 (mod 8),")
print("  и тогда центр c4 = s(1+t^2) не может быть квадратом ни при каком t ∈ P^1(Q_2).")
print("  Это и есть утверждение Бойера «при чётном n красные клетки не квадраты»:")
print("  n — ИНДЕКС рекурренты, красные клетки — c0,c4,c8 (центр и две свободные), а не c2,c6.")

print()
print("="*78)
print("F. КОНТРОЛЬНЫЕ ЧИСЛЕННЫЕ ФАКТЫ ДЛЯ ОТЧЁТА")
print("="*78)
for (m0,n0) in [(13,8),(15,8),(19,16),(16,5)]:
    s_ = QQ(m0^2+n0^2)/2
    print(f"  (m,n)=({m0},{n0}): при t=1  c0=c4=c8={m0^2+n0^2}  (mod 8 = {(m0^2+n0^2)%8}),"
          f"  c2=(m-n)^2={(m0-n0)^2}, c6=(m+n)^2={(m0+n0)^2}")
print("  m^2+n^2 ≡ 1 (mod 8) при m нечётном, n ≡ 0 (mod 4)  =>  m^2+n^2 — квадрат в Q_2")
print("  m^2+n^2 ≡ 5 (mod 8) при m нечётном, n ≡ 2 (mod 4)  =>  НЕ квадрат в Q_2")
print("  m^2+n^2 ≡ 2 (mod 8) при m,n нечётных              =>  v_2 = 1, НЕ квадрат в Q_2")
print()
print("Готово.")
