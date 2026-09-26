# (а) Независимая символьная проверка тождеств распределения и F2-структуры
# записки Codex «восемь линий» (RESULT_EIGHT_LINES_BRIDGE_2026-09-25_FROM_CODEX.md).
# Запуск: env DOT_SAGE=/tmp/claude_verify_bridge_sage python3 a_identities.py
from sage.all import *
import itertools, json
from pathlib import Path
OUT = Path(__file__).resolve().parent
res = {}

# --- 1. Общий магический квадрат: решаем систему сами, а не подставляем форму Codex ---
R = PolynomialRing(QQ, ['x%d' % i for i in range(9)] + ['S'])
X = R.gens()[:9]; S = R.gens()[9]
M = [[X[0], X[1], X[2]], [X[3], X[4], X[5]], [X[6], X[7], X[8]]]
lines_pos = [[(0,0),(0,1),(0,2)], [(1,0),(1,1),(1,2)], [(2,0),(2,1),(2,2)],
             [(0,0),(1,0),(2,0)], [(0,1),(1,1),(2,1)], [(0,2),(1,2),(2,2)],
             [(0,0),(1,1),(2,2)], [(0,2),(1,1),(2,0)]]
Amat = matrix(QQ, [[1 if (i//3, i%3) in L else 0 for i in range(9)] + [-1] for L in lines_pos])
ker = Amat.right_kernel()
res['magic_solution_dim'] = int(ker.dimension())
res['incidence_rank_Q'] = int(matrix(QQ, [r[:9] for r in Amat.rows()]).rank())
res['incidence_rank_F2'] = int(matrix(GF(2), [[int(c) for c in r[:9]] for r in Amat.rows()]).rank())
res['membership'] = [int(sum(1 for L in lines_pos if (i//3, i%3) in L)) for i in range(9)]
# форма Codex
P = PolynomialRing(QQ, ['a', 'b', 'c']); a, b, c = P.gens()
mag = [a-b, a+b+c, a-c, a+b-c, a, a-b+c, a+c, a-b-c, a+b]
vecs = [vector(QQ, [mag[i].monomial_coefficient(g) for i in range(9)] + [3*(g == a)]) for g in (a, b, c)]
assert all(v in ker for v in vecs)
assert matrix(QQ, vecs).rank() == 3 == ker.dimension()   # форма Codex = всё пространство решений
res['codex_form_is_general'] = True
corners = [mag[i] for i in (0, 2, 6, 8)]; edges = [mag[i] for i in (1, 3, 5, 7)]
assert sum(corners) == 4*a and sum(edges) == 4*a
assert all(mag[i] + mag[8-i] == 2*a for i in range(4))
C2 = sum((v-a)**2 for v in corners); E2 = sum((v-a)**2 for v in edges)
C4 = sum((v-a)**4 for v in corners); E4 = sum((v-a)**4 for v in edges)
assert C2 == 2*(b**2+c**2) and E2 == 2*C2
assert E4 - 2*C4 == 24*b**2*c**2
res['C2'] = str(C2); res['E2'] = str(E2); res['E4_minus_2C4'] = str(E4-2*C4)
# корни: значение клетки = r^2, r^4 = значение^2; центр a = r_center^2
roots_id = sum(v**2 for v in edges) - 2*sum(v**2 for v in corners) + 4*a**2
assert roots_id == 0
res['root4_identity_zero'] = True
# Тождество как функция от корней с a = r5^2: проверка через кольцо корней
Rr = PolynomialRing(QQ, ['r%d' % i for i in range(9)]); r = Rr.gens()
# не независимое условие: выводится из линейных. Проверяем, что оно в идеале, порождённом 8 линейными равенствами в r_i^2
Sm = Rr.gen(4)**2
lin = [sum(r[3*i+j]**2 for (i, j) in L) - 3*Sm for L in lines_pos]
I = Rr.ideal(lin)
f = sum(r[i]**4 for i in (1, 3, 5, 7)) - 2*sum(r[i]**4 for i in (0, 2, 6, 8)) + 4*r[4]**4
res['root4_identity_in_linear_ideal'] = bool(f in I)
# положительность: |b|+|c|<a  (min клетки = a-|b|-|c|)
res['positivity_min_cell'] = 'a-|b|-|c| (клетки a±b±c присутствуют)'

# --- 2. F2: ядро восьми AP-линий, классы T, L ---
GRID = [[0,1,2],[3,4,5],[6,7,8],[0,3,6],[1,4,7],[2,5,8],[0,4,8],[2,4,6]]
A2 = matrix(GF(2), [[int(i in L) for i in range(9)] for L in GRID])
K2 = A2.right_kernel()
kvecs = sorted(tuple(int(x) for x in v) for v in K2)
expected = sorted(tuple([s, t, s^t, t, 0, t, s^t, t, s]) for s in (0, 1) for t in (0, 1))
assert kvecs == expected
res['grid_kernel_dim'] = int(K2.dimension())
# AP-сетка: клетка (i,j) = a + i*b + j*c, i,j = -1,0,1, построчно
apgrid = [a + i*b + j*c for i in (-1, 0, 1) for j in (-1, 0, 1)]
idx = {v: k for k, v in enumerate(apgrid)}
assert set(apgrid) == set(mag)
T = [a-b, a+b+c, a-c]; L = [a-b, a+b-c, a+c]
Tv = vector(GF(2), [int(apgrid[k] in T) for k in range(9)]); Lv = vector(GF(2), [int(apgrid[k] in L) for k in range(9)])
for s_ in (0, 1):
    for t_ in (0, 1):
        v = vector(GF(2), [s_, t_, s_^t_, t_, 0, t_, s_^t_, t_, s_])
        assert Tv.dot_product(v) == s_ and Lv.dot_product(v) == (s_ ^ t_)
res['class_T_eq_s_class_L_eq_s_plus_t'] = True
assert matrix(GF(2), list(A2.rows()) + [Tv, Lv]).rank() == 9
# центральные линии магической расстановки = центральные линии AP-сетки (как множества)
magic_lines = [[mag[3*i+j] for (i, j) in Lp] for Lp in lines_pos]
grid_lines = [[apgrid[k] for k in G] for G in GRID]
common = [sorted(map(str, m)) for m in magic_lines if any(set(m) == set(g) for g in grid_lines)]
res['common_lines_count'] = len(common)
assert len(common) == 4 and all(a in m for m in [[mag[3*i+j] for (i, j) in Lp] for Lp in lines_pos[4:5]])
# формула sqrt(A22)
Q9 = PolynomialRing(QQ, ['A%d' % i for i in range(9)]); Av = Q9.gens()
R1 = Av[0]*Av[1]*Av[2]; R3 = Av[6]*Av[7]*Av[8]; C2_ = Av[1]*Av[4]*Av[7]; D1 = Av[0]*Av[4]*Av[8]; D2 = Av[2]*Av[4]*Av[6]
num2 = R1*R3*C2_*D1*D2
den = Av[4]*R1*R3
assert num2 == den**2 * Av[4]
res['sqrtA22_formula_squared_ok'] = True

# --- 3. Пример по модулю 37 ---
p = 37; av, bv, cv = 1, 14, 16
g = [(av + i*bv + j*cv) % p for i in (-1, 0, 1) for j in (-1, 0, 1)]
sq = set(x*x % p for x in range(1, p))
prods = [prod(g[k] for k in G) % p for G in GRID]
Tm = (av-bv)*(av+bv+cv)*(av-cv) % p; Lm = (av-bv)*(av+bv-cv)*(av+cv) % p
res['mod37'] = dict(grid=g, distinct=len(set(g)) == 9, nonzero=0 not in g,
                    eight_products_residues=all(x in sq for x in prods), T=int(Tm), T_residue=Tm in sq, L=int(Lm), L_residue=Lm in sq)
assert g == [8, 24, 3, 22, 1, 17, 36, 15, 31] and all(x in sq for x in prods) and Tm == 14 and Tm not in sq and Lm == 36

# --- 4. Род: независимо через Риман-Гурвиц для (Z/2)^r-накрытий с явными подсчётами ---
# при a=1, b=z, c=kz: 8 нецентральных линейных форм, различные нули; бесконечность
An = [[int(k in G) for k in range(9) if k != 4] for G in GRID]
rk = matrix(GF(2), An).rank()
def genus_of_subcover(rows):
    Mx = matrix(GF(2), rows); r_ = Mx.rank(); d = 2**r_
    br = 0
    for j in range(8):
        if any(row[j] for row in rows): br += 1
    if any(sum(row) % 2 for row in rows): br += 1
    return r_, d, 1 + d*(br/QQ(4) - 1) if True else None
r8, d8, g8 = genus_of_subcover(An)
rX, dX, gX = genus_of_subcover([[int(i == j) for j in range(8)] for i in range(8)])
res['genus'] = dict(Y_rank=int(r8), Y_degree=int(d8), Y_genus=str(g8), X_degree=int(dX), X_genus=str(gX))
assert (d8, g8, dX, gX) == (64, 81, 256, 321)
assert 2*gX - 2 == 4*(2*g8 - 2)
# носители автоморфизмов X/Y
supp = sorted(sum(int(x) for x in (list(v[:4]) + list(v[5:]))) for v in K2 if v != 0)
res['X_over_Y_supports'] = supp
assert supp == [4, 6, 6]

(OUT / 'a_identities.json').write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str))
print('ALL (a) CHECKS PASSED')
