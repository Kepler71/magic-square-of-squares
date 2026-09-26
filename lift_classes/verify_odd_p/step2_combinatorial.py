#!/usr/bin/env python3
"""Скептическая проверка Шага 2 (нечётное p) теоремы Codex об остаточных классах.
Чисто комбинаторные части (без рациональных точек):
  Part 1. Лемма о прямой: при z-единице множество клеток 1+(ir+js)z = 0 mod p
          содержит <= 3 клетки сетки {-1,0,1}^2 (p нечётно, (r,s) != (0,0) mod p).
  Part 2. Пространство чётностей 8 линий арифметической сетки (центр = 0):
          ровно матрица S,U,S+U; веса 0,4,6,6; противоположные клетки равны.
  Part 3. Абстрактная модель полюса: НАДМНОЖЕСТВО всех реальных конфигураций
          оценок (v(r),v(s),v(r+s),v(r-s), n, выбор (k,0)/(0,k) при v(lambda)=n).
          Проверяем: если вектор чётностей в пространстве, то
          S=1 => p | (r-s),  S+U=1 => p | (r+s);  p|r или p|s => S=U=0.
Все циклы конечны и с явными границами.
"""
import itertools, sys, time

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
IDX = {c: k for k, c in enumerate(CELLS)}
LINES = ([[(i, j) for j in (-1, 0, 1)] for i in (-1, 0, 1)] +
         [[(i, j) for i in (-1, 0, 1)] for j in (-1, 0, 1)] +
         [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]])
assert len(LINES) == 8

def primes_upto(N):
    return [q for q in range(3, N + 1) if all(q % d for d in range(2, int(q**0.5) + 1))]

# ---------------- Part 1 ----------------
t0 = time.time()
maxcnt = {}
for p in primes_upto(211):
    m = 0
    for a in range(p):          # a = r z mod p
        for b in range(p):      # b = s z mod p
            if a == 0 and b == 0:
                continue
            cnt = sum(1 for (i, j) in CELLS if (1 + i * a + j * b) % p == 0)
            m = max(m, cnt)
    maxcnt[p] = m
    assert m <= 3, (p, m)
print("Part1 line lemma: max zero-cells per prime (p<=211):",
      {p: maxcnt[p] for p in list(maxcnt)[:6]}, "... all <=3; time %.1fs" % (time.time() - t0), flush=True)

# ---------------- Part 2 ----------------
def in_space(vec):
    return all(sum(vec[IDX[c]] for c in L) % 2 == 0 for L in LINES)

sols = [v for v in itertools.product((0, 1), repeat=9) if in_space(v)]
sols_c0 = [v for v in sols if v[IDX[(0, 0)]] == 0]
print("Part2: solutions of 8 line eqs over F2: %d total, %d with centre 0" % (len(sols), len(sols_c0)))

def matrix(S, U):
    m = {(-1, -1): S, (-1, 0): U, (-1, 1): (S + U) % 2,
         (0, -1): U, (0, 0): 0, (0, 1): U,
         (1, -1): (S + U) % 2, (1, 0): U, (1, 1): S}
    return tuple(m[c] for c in CELLS)

expected = {matrix(S, U) for S in (0, 1) for U in (0, 1)}
assert set(sols_c0) == expected, "matrix claim FAILS"
weights = sorted(sum(v) for v in sols_c0)
assert weights == [0, 4, 6, 6], weights
for v in sols_c0:
    for (i, j) in CELLS:
        assert v[IDX[(i, j)]] == v[IDX[(-i, -j)]]
# no nonzero solution supported on <=3 cells, and in particular on any line mod p
assert all(sum(v) == 0 or sum(v) >= 4 for v in sols_c0)
print("Part2: matrix S,U,S+U confirmed; weights", weights, "; opposite cells equal", flush=True)
# note: without the centre constraint
extra = [v for v in sols if v[IDX[(0, 0)]] == 1]
print("Part2 (info): solutions with centre parity 1:", len(extra), "(irrelevant: centre = 1)", flush=True)

# T and L classes
T_cells = [(-1, 0), (1, 1), (0, -1)]
L_cells = [(-1, 0), (1, -1), (0, 1)]
for S in (0, 1):
    for U in (0, 1):
        v = matrix(S, U)
        assert sum(v[IDX[c]] for c in T_cells) % 2 == S
        assert sum(v[IDX[c]] for c in L_cells) % 2 == (S + U) % 2
print("Part2: [T]=S, [L]=S+U confirmed", flush=True)

# ---------------- Part 3 ----------------
# lambda for cell (i,j) is i r + j s; pairs: r: (+-1,0); s: (0,+-1); r+s: (1,1),(-1,-1); r-s: (1,-1),(-1,1)
PAIRS = {'r': [(1, 0), (-1, 0)], 's': [(0, 1), (0, -1)],
         'p': [(1, 1), (-1, -1)], 'm': [(1, -1), (-1, 1)]}

def realizable(a, b, c, d, constrained=True):
    if not constrained:
        return True
    if min(a, b) != 0:
        return False
    if a > 0 and (b, c, d) != (0, 0, 0):
        return False
    if b > 0 and (a, c, d) != (0, 0, 0):
        return False
    if a == 0 and b == 0 and min(c, d) != 0:
        return False
    return True

def run_model(constrained, VMAX=7, NMAX=6, KMAX=3):
    checked = inspace = viol = 0
    nontriv = {'S': 0, 'SU': 0}
    for a, b, c, d in itertools.product(range(VMAX + 1), repeat=4):
        if not realizable(a, b, c, d, constrained):
            continue
        val = {'r': a, 's': b, 'p': c, 'm': d}
        for n in range(1, NMAX + 1):
            eq = [lam for lam in 'rspm' if val[lam] == n]
            # for each equal-valuation pair: choose which cell carries k, and k in 0..KMAX
            choices = list(itertools.product(*[[(0, k) for k in range(KMAX + 1)] +
                                               [(k, 0) for k in range(1, KMAX + 1)] for _ in eq]))
            for ch in choices:
                vv = {}
                for lam in 'rspm':
                    e = val[lam]
                    c1, c2 = PAIRS[lam]
                    if e < n:
                        vv[c1] = vv[c2] = e - n
                    elif e > n:
                        vv[c1] = vv[c2] = 0
                    else:
                        k1, k2 = ch[eq.index(lam)]
                        vv[c1], vv[c2] = k1, k2
                vv[(0, 0)] = 0
                P = tuple(vv[cc] % 2 for cc in CELLS)
                checked += 1
                if not in_space(P):
                    continue
                inspace += 1
                S = P[IDX[(1, 1)]]
                U = P[IDX[(1, 0)]]
                SU = P[IDX[(1, -1)]]
                assert SU == (S + U) % 2
                bad = False
                if S == 1:
                    nontriv['S'] += 1
                    if d == 0:
                        bad = True
                if SU == 1:
                    nontriv['SU'] += 1
                    if c == 0:
                        bad = True
                if (a > 0 or b > 0) and (S or U):
                    bad = True
                if bad:
                    viol += 1
                    if constrained:
                        raise AssertionError(("violation", a, b, c, d, n, ch, P))
    return checked, inspace, viol, nontriv

t0 = time.time()
res = run_model(True)
print("Part3 (constrained = real odd-p valuations): configs %d, in-space %d, violations %d, nontrivial %s; %.1fs"
      % (res[0], res[1], res[2], res[3], time.time() - t0), flush=True)
t0 = time.time()
res_u = run_model(False, VMAX=4, NMAX=4, KMAX=2)
print("Part3 negative control (valuation constraints dropped): configs %d, in-space %d, 'violations' %d; %.1fs"
      % (res_u[0], res_u[1], res_u[2], time.time() - t0), flush=True)
print("Part3 conclusion: the model check is non-vacuous iff the unconstrained run shows violations:",
      res_u[2] > 0, flush=True)
print("ALL COMBINATORIAL CHECKS PASSED", flush=True)
