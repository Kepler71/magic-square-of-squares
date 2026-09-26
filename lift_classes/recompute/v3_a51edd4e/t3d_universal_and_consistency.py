# -*- coding: utf-8 -*-
"""t3d: (A) универсальный образ I_p = объединение I_p(r0,s0) по всем (r0,s0) в F_p^2 \ {0}, p = 3 mod 4;
(B) контроль согласия строгой верхней оценки (t3c) с ТОЧНЫМ выборочным образом на реальных рациональных z
    (t3-метод, классы в Q_p точны): выборка должна лежать в оценке; фиксируем, совпадают ли они;
(C) независимая сверка Sage для 126/451 при p=7: все z = 7^e a, a mod 7^5 (единицы), e in [-4,5],
    восемь произведений через QQ.is_padic_square; класс L -- тоже через Sage.
(D) F_2-структура: решения 8 уравнений классов (своя линейная алгебра).
Все циклы ограничены явно.
"""
import os, sys, json, time, random, itertools
os.environ.setdefault("DOT_SAGE", "/tmp/sage_a51edd4e")
from fractions import Fraction
from t3c_unit_filter import image_upper, chi, primes_upto
from t3_local import cell_classes, eight_ok, samples_a, T_CELLS, L_CELLS, CELLS, LINES

t0 = time.time()
out = {}
# (D) F_2
sol = []
for bits in range(1 << 8):
    cells = [c for c in CELLS if c != (0, 0)]
    v = {c: (bits >> k) & 1 for k, c in enumerate(cells)}; v[(0, 0)] = 0
    if all((v[l[0]] ^ v[l[1]] ^ v[l[2]]) == 0 for l in LINES):
        sol.append(v)
ok_pattern = True
for v in sol:
    S, U = v[(-1, -1)], v[(-1, 0)]
    pat = {(-1, -1): S, (-1, 0): U, (-1, 1): S ^ U, (0, -1): U, (0, 0): 0, (0, 1): U, (1, -1): S ^ U, (1, 0): U, (1, 1): S}
    ok_pattern &= (v == pat)
    ok_pattern &= ((v[T_CELLS[0]] ^ v[T_CELLS[1]] ^ v[T_CELLS[2]]) == S)
    ok_pattern &= ((v[L_CELLS[0]] ^ v[L_CELLS[1]] ^ v[L_CELLS[2]]) == (S ^ U))
out["F2"] = dict(n_solutions=len(sol), weights=sorted(sum(v.values()) for v in sol), pattern_ok=ok_pattern)
print("(D) F2:", out["F2"], flush=True)

# (A) универсальные образы
univ = {}
for p in [q for q in primes_upto(200) if q % 4 == 3 and q >= 7]:
    U = set()
    for r0 in range(p):
        for s0 in range(p):
            if r0 == 0 and s0 == 0: continue
            U |= image_upper(r0, s0, p)
            if len(U) == 4: break
        if len(U) == 4: break
    univ[p] = sorted(U)
    print("(A) p=%d: объединение I_p по всем (r,s) mod p = %s" % (p, sorted(U)), flush=True)
out["universal_I_p"] = {str(k): v for k, v in univ.items()}

# (B) согласие оценки и точной выборки
rng = random.Random(11)
slopes = [(126, 451), (-126, 451), (73, 362), (265, 298), (12, 85), (2, 3), (4, 9), (8, 9), (19, 30), (5, 8)]
cons = []
for (r, s) in slopes:
    for p in [q for q in primes_upto(80) if q % 4 == 3 and q >= 7]:
        emp = set()
        for e in range(-4, 5):
            for a in samples_a(p, r, s, rng, 6000):
                C = cell_classes(r, s, p, e, a)
                if C is None or not eight_ok(C): continue
                cT = C[T_CELLS[0]] ^ C[T_CELLS[1]] ^ C[T_CELLS[2]]
                cL = C[L_CELLS[0]] ^ C[L_CELLS[1]] ^ C[L_CELLS[2]]
                if (cT & 1) or (cL & 1):
                    emp.add(("ODD", cT, cL)); continue
                emp.add((1 if cT == 0 else -1, 1 if cL == 0 else -1))
        up = image_upper(r, s, p)
        sub = emp <= up
        cons.append(dict(slope="%d/%d" % (r, s), p=p, empirical=sorted(map(str, emp)), upper=sorted(up),
                         subset=sub, equal=(emp == up)))
        if not sub or emp != up:
            print("(B) %d/%d p=%d: выборка %s, оценка %s, subset=%s" % (r, s, p, sorted(map(str, emp)), sorted(up), sub), flush=True)
nsub = sum(c["subset"] for c in cons); neq = sum(c["equal"] for c in cons)
print("(B) пар (наклон,p): %d; выборка внутри оценки: %d; совпадение: %d" % (len(cons), nsub, neq), flush=True)
out["consistency"] = dict(n=len(cons), subset=nsub, equal=neq, rows=cons)

# (C) Sage-сверка 126/451 при p=7
from sage.all import QQ
r, s, p = 126, 451, 7
cnt = passed = Lnontriv = 0
for e in range(-4, 6):
    for a in range(1, 7 ** 5):
        if a % 7 == 0: continue
        z = Fraction(a) * Fraction(7) ** e
        f = {c: 1 + (c[0] * r + c[1] * s) * z for c in CELLS}
        if any(v == 0 for v in f.values()): continue
        cnt += 1
        if all(QQ(f[l[0]] * f[l[1]] * f[l[2]]).is_padic_square(7) for l in LINES):
            passed += 1
            L = f[L_CELLS[0]] * f[L_CELLS[1]] * f[L_CELLS[2]]
            T = f[T_CELLS[0]] * f[T_CELLS[1]] * f[T_CELLS[2]]
            if not (QQ(L).is_padic_square(7) and QQ(T).is_padic_square(7)):
                Lnontriv += 1
    print("(C) e=%d: проверено %d, прошло %d, нетривиальный [T] или [L]: %d  t=%.1fs" % (e, cnt, passed, Lnontriv, time.time() - t0), flush=True)
out["sage_126_451_p7"] = dict(tested=cnt, passed=passed, nontrivial_TL=Lnontriv,
                             chi7_577=int(QQ(577).is_padic_square(7)))
print("(C) 577 -- квадрат в Q_7?", QQ(577).is_padic_square(7))
json.dump(out, open("t3d_universal_and_consistency.json", "w"), indent=1, default=str)
print("время %.1f" % (time.time() - t0))
