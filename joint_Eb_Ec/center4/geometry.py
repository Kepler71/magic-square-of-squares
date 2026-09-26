# -*- coding: utf-8 -*-
"""Геометрия «общего центра»: роды кривых (Риман–Гурвиц + контроль рядом Гильберта в Singular),
кривая D_4 (общий x на четырёх E_u, только произведения) и её якобиан (разложение по характерам).
Запуск: env DOT_SAGE=<scratch> python3 geometry.py"""
from sage.all import *
import json, time
t0 = time.time()
out = {}

def rh_genus(nforms, nbranch):
    # (Z/2)^{n}/diag-накрытие P^1 (n квадратных корней из линейных форм), ветвление в nbranch точках, индекс 2
    G = 2 ** (nforms - 1)
    return (G * (-2) + nbranch * (G // 2)) // 2 + 1

out['RH'] = {
  'slope_fixed_projective_9forms': rh_genus(9, 9),
  # аффинная карта a=1: 8 форм 1+s_t b, ветвление в 8 конечных точках И в b=oo (все 8 корней), группа (Z/2)^8
  'slope_fixed_affine_a=1_check': (2**8*(-2) + 9*2**7)//2 + 1,
  'b_over_a_fixed_projective_7forms': rh_genus(7, 7),
  # аффинная карта a=1, b фиксировано: 6 форм 1±c, 1±b±c, ветвление в 6 точках и в c=oo, группа (Z/2)^6
  'b_fixed_affine_a=1_check': (2**6*(-2) + 7*2**5)//2 + 1,
}
print(out['RH'])

# контроль: ряд Гильберта кривой фиксированного наклона k=c/b=3/7 в P^8 (9 квадратичных соотношений на 2 параметра -> 7 квадрик)
k = QQ(3)/7
R = PolynomialRing(QQ, ['z%d' % t for t in range(9)], order='degrevlex')
zs = R.gens()
forms = [(1, i + j * k) for i in (-1, 0, 1) for j in (-1, 0, 1)]   # клетка = a + (i + j k) b  -> (коэф a, коэф b)
# соотношения: z_t^2 = a + s_t b; исключаем a, b: z_t^2 - z_0'^2 ... берём базис ядра
M = matrix(QQ, [[f[0], f[1]] for f in forms])      # 9x2
K = M.left_kernel().basis()                          # 7 векторов: sum c_t (a + s_t b) = 0
I = R.ideal([sum(c[t] * zs[t] ** 2 for t in range(9)) for c in K])
hp = I.hilbert_polynomial()
out['slope_curve_k=3/7'] = dict(n_quadrics=len(K), dim=I.dimension(), hilbert_polynomial=str(hp))
co = hp.list()
deg = co[1] if len(co) == 2 else None
g = 1 - co[0]
out['slope_curve_k=3/7'].update(degree=int(deg), genus=int(g))
print('кривая наклона 3/7:', out['slope_curve_k=3/7'], f'{time.time()-t0:.1f}s', flush=True)

# D_4: y_u^2 = x(x^2-u^2), u in {b,c,b+c,b-c}; разложение якобиана по характерам (Z/2)^4
b, c = 7, 3   # «типичная» пара (без совпадений)
us = [b, c, b + c, b - c]
Px = PolynomialRing(QQ, 'x'); x = Px.gen()
comp = []
from itertools import combinations
tot = 0
for r in range(1, 5):
    for S in combinations(range(4), r):
        f = prod(x**2 - us[t]**2 for t in S) * (x if r % 2 else 1)
        C = HyperellipticCurve(f)
        comp.append(dict(S=[['b','c','b+c','b-c'][t] for t in S], deg=int(f.degree()), genus=int(C.genus())))
        tot += C.genus()
out['D4'] = dict(b=b, c=c, genus_sum=int(tot), RH_genus=int((16*(-2) + 10*8)//2 + 1), components=comp)
print('D4: сумма родов компонент', tot, ' RH:', out['D4']['RH_genus'])
# j-инвариант F_{u,v}: w^2=(x^2-u^2)(x^2-v^2), якобиан y^2 = X(X+(u+v)^2)(X+(u-v)^2) — зависит только от наклона
kk = polygen(QQ, 'k')
def jF(u, v):
    A = (u + v) ** 2; B = (u - v) ** 2
    # y^2 = X^3 + (A+B) X^2 + A B X
    a2, a4 = A + B, A * B
    c4 = 16 * (a2 ** 2 - 3 * a4)
    disc = 16 * a4 ** 2 * (a2 ** 2 - 4 * a4)
    return c4 ** 3 / disc
out['j_F_bc_as_function_of_k=c/b'] = str(factor(jF(1, kk)))
E1 = EllipticCurve([0, (1+kk.parent()(0)) , 0, 0, 0]) if False else None
# контроль: сравнение с Sage для k=3/7
Ft = EllipticCurve(QQ, [0, (1 + k) ** 2 + (1 - k) ** 2, 0, (1 + k) ** 2 * (1 - k) ** 2, 0])
out['j_check_k=3/7'] = dict(sage=str(Ft.j_invariant()), formula=str(jF(1, k)), equal=bool(Ft.j_invariant() == jF(1, k)))
print(out['j_check_k=3/7'])
out['time_s'] = round(time.time() - t0, 1)
json.dump(out, open('geometry.json', 'w'), indent=1, default=str)
print(json.dumps(out, indent=1, default=str)[:3000])
