# Независимая проверка шага 1: матрица классов над F2 и отождествление T,L.
# Все циклы конечны: 2^9 векторов, 8 линий.
from itertools import product
import sympy as sp

cells = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]   # сетка f_ij = 1 + i*b + j*c
idx = {c: k for k, c in enumerate(cells)}
grid_lines = [[(i, j) for j in (-1, 0, 1)] for i in (-1, 0, 1)] \
           + [[(i, j) for i in (-1, 0, 1)] for j in (-1, 0, 1)] \
           + [[(t, t) for t in (-1, 0, 1)], [(t, -t) for t in (-1, 0, 1)]]
assert len(grid_lines) == 8

def rank_f2(rows):
    rows = [r[:] for r in rows]; rk = 0; n = len(rows[0])
    for col in range(n):
        piv = next((k for k in range(rk, len(rows)) if rows[k][col]), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for k in range(len(rows)):
            if k != rk and rows[k][col]:
                rows[k] = [(x + y) % 2 for x, y in zip(rows[k], rows[rk])]
        rk += 1
    return rk

M = [[1 if c in L else 0 for c in cells] for L in grid_lines]
print('rank over F2 of 8 grid-line equations on 9 cells:', rank_f2(M))
sols = [v for v in product((0, 1), repeat=9)
        if all(sum(v[idx[c]] for c in L) % 2 == 0 for L in grid_lines)]
print('number of solutions:', len(sols), '(expect 4 = 2^2)')
assert len(sols) == 4
assert all(v[idx[(0, 0)]] == 0 for v in sols), 'center forced trivial'
# Параметризация S = класс угла (-1,-1), U = класс (-1,0)
for v in sols:
    S, U = v[idx[(-1, -1)]], v[idx[(-1, 0)]]
    pred = {(-1,-1): S, (-1,0): U, (-1,1): (S+U)%2, (0,-1): U, (0,0): 0, (0,1): U,
            (1,-1): (S+U)%2, (1,0): U, (1,1): S}
    assert all(v[idx[c]] == pred[c] for c in cells), v
print('matrix [[S,U,S+U],[U,0,U],[S+U,U,S]] confirmed for all 4 solutions')
# над произвольной F2-группой (Q*/Q*^2 бесконечномерна) утверждение покоординатное: ок.

# Магическая расстановка из моста (a=1): проверяем, что она магическая и что T,L её линии.
a, b, c, z, r, s = sp.symbols('a b c z r s')
magic = [[a-b, a+b+c, a-c], [a+b-c, a, a-b+c], [a+c, a-b-c, a+b]]
mlines = [row for row in magic] + [[magic[i][j] for i in range(3)] for j in range(3)] \
       + [[magic[t][t] for t in range(3)], [magic[t][2-t] for t in range(3)]]
assert all(sp.expand(sum(L) - 3*a) == 0 for L in mlines)
gridvals = {(i, j): a + i*b + j*c for (i, j) in cells}
assert sorted(map(str, [sp.expand(x) for row in magic for x in row])) == \
       sorted(map(str, [sp.expand(x) for x in gridvals.values()]))
def cls(expr_list):
    # класс произведения = сумма классов клеток (как векторы в (S,U))
    tot = [0, 0]
    inv = {sp.expand(v): k for k, v in gridvals.items()}
    for e in expr_list:
        k = inv[sp.expand(e)]
        S_, U_ = {(-1,-1):(1,0),(1,1):(1,0),(-1,1):(1,1),(1,-1):(1,1),
                  (-1,0):(0,1),(1,0):(0,1),(0,-1):(0,1),(0,1):(0,1),(0,0):(0,0)}[k]
        tot = [(tot[0]+S_)%2, (tot[1]+U_)%2]
    return tuple(tot)
names = ['row1','row2','row3','col1','col2','col3','diag','anti']
for nm, L in zip(names, mlines):
    print('magic', nm, [str(e) for e in L], 'class (S,U)-coeffs =', cls(L))
T = (1 - r*z)*(1 + (r+s)*z)*(1 - s*z)
Lp = (1 - r*z)*(1 + (r-s)*z)*(1 + s*z)
sub = {a: 1, b: r*z, c: s*z}
top = sp.prod([e.subs(sub) for e in magic[0]])
left = sp.prod([magic[i][0].subs(sub) for i in range(3)])
assert sp.expand(T - top) == 0 and sp.expand(Lp - left) == 0
print('T == magic top row, L == magic left column (a=1,b=rz,c=sz): OK')
print('[T] =', cls(magic[0]), ' [L] =', cls([magic[i][0] for i in range(3)]))
assert cls(magic[0]) == (1, 0) and cls([magic[i][0] for i in range(3)]) == (1, 1)
# Все 4 внешние магические линии: классы S,S,S+U,S+U; центральные: 0.
# Восстановление: T,L квадраты => S=0,S+U=0 => U=0.
print('ALL F2 CHECKS PASSED')
