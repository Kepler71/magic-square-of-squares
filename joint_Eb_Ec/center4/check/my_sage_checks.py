# Проверяющий: Sage-проверки к утверждениям 2, 8, 9, 10.
from sage.all import *
import json
out = {}
# --- утверждение 8: Саллоус, 2P, 3P
F = json.load(open('/home/kep/magicKube/fable_symmetry/s1_sallows.json'))
a, b, c, lam = F['a'], F['b'], F['c'], F['lam']
sq = [[a-b+lam, a+b+c, a-c], [a+b-c, a+lam, a-b+c], [a+c, a-b-c, a+b+lam]]
out['sallows_all_squares'] = all(ZZ(x).is_square() for r in sq for x in r)
rows = [sum(r) for r in sq]; cols = [sum(sq[i][j] for i in range(3)) for j in range(3)]
d1 = sq[0][0]+sq[1][1]+sq[2][2]; d2 = sq[0][2]+sq[1][1]+sq[2][0]
out['sallows_lines'] = dict(rows=rows, cols=cols, diag=d1, anti=d2)
six = lambda S: [S[0][1], S[0][2], S[1][0], S[1][2], S[2][0], S[2][1]]
for k, S in enumerate(F['new_seven_line_squares']):
    S = [[QQ(x) for x in r] for r in S]
    rs = [sum(r) for r in S]; cs = [sum(S[i][j] for i in range(3)) for j in range(3)]
    dd = S[0][0]+S[1][1]+S[2][2]; aa = S[0][2]+S[1][1]+S[2][0]
    allsq = all(x.is_square() for r in S for x in r)
    out[f'sq{k+2}P'] = dict(six_same=six(S) == six(sq), all_squares=allsq, magic_lines=sum(1 for v in rs+cs+[dd, aa] if v == rs[0]), center_minus_a=str(S[1][1]-a))
# три прямые через (неквадратный) центр 4420 и /4
A4 = 4420
us = [u for u in range(1, A4) if ZZ(A4-u).is_square() and ZZ(A4+u).is_square()]
out['lines_through_4420'] = us
# ранги
def rk(N):
    E = EllipticCurve([-N**2, 0]); N0 = ZZ(N).squarefree_part()
    E0 = EllipticCurve([-N0**2, 0])
    r = E0.rank(only_use_mwrank=True)
    return dict(N=int(N0), N_mod8=int(N0 % 8), root_number=int(E0.root_number()), rank=int(r), sel2=int(E0.selmer_rank()))
bb, cc = b, c
out['ranks'] = {k: rk(abs(v)) for k, v in dict(b=bb, c=cc, bpc=bb+cc, bmc=bb-cc, ac=A4*cc, abpc=A4*(bb+cc), abmc=A4*(bb-cc)).items()}
# AB1: центр 425^2, B=41496, C=-138600 (данные автора; сверяем, какие прямые полные)
A1, B1, C1 = 425**2, 41496, -138600
cells = {(i, j): A1 + i*B1 + j*C1 for i in (-1, 0, 1) for j in (-1, 0, 1)}
out['AB1_square_cells'] = sorted([str(k) for k, v in cells.items() if v > 0 and ZZ(v).is_square()])
out['AB1_full_lines'] = [u for u in (B1, C1, B1+C1, B1-C1) if ZZ(A1-u).is_square() and ZZ(A1+u).is_square()]
# --- утверждение 9/10: формула корневого числа и 2-чётность на малых N
bad_w = []; bad_par = []
for N in range(1, 400):
    if not ZZ(N).is_squarefree(): continue
    E = EllipticCurve([-N**2, 0]); w = E.root_number()
    pred = 1 if N % 8 in (1, 2, 3) else -1
    if w != pred: bad_w.append(N)
    s2 = E.selmer_rank()  # mwrank: ранг 2-Сельмера без 2-кручения
    if (-1)**s2 != w: bad_par.append((N, s2, w))
out['rootnumber_formula_fail_N_lt_400'] = bad_w
out['two_parity_fail_N_lt_400'] = bad_par
E1 = EllipticCurve([-1, 0]); out['selmer_rank_E1'] = int(E1.selmer_rank())
# --- утверждение 2: якобиан w^2=(x^2-u^2)(x^2-v^2) и j-инвариант
R = PolynomialRing(QQ, 'k'); k = R.gen()
jac_ok = []
for (u, v) in [(3, 7), (1, 5), (2, 9), (5, 11)]:
    # инварианты квартики y^2 = x^4 + 0 x^3 + c2 x^2 + 0 x + e
    a4, b3, c2, d1, e0 = 1, 0, -(u*u+v*v), 0, u*u*v*v
    I = 12*a4*e0 - 3*b3*d1 + c2**2
    J = 72*a4*c2*e0 + 9*b3*c2*d1 - 27*a4*d1**2 - 27*e0*b3**2 - 2*c2**3
    Ej = EllipticCurve([-27*I, -27*J])
    Ec = EllipticCurve([0, (u+v)**2 + (u-v)**2, 0, (u+v)**2*(u-v)**2, 0])
    kk = QQ(v)/QQ(u)
    jf = 16*(kk**4+14*kk**2+1)**3/(kk**2*(kk**2-1)**4)
    jac_ok.append(dict(u=u, v=v, iso_QQ=Ej.is_isomorphic(Ec), j_quartic=str(Ej.j_invariant()), j_formula=str(jf), j_eq=Ej.j_invariant() == jf))
out['F_uv_jacobian'] = jac_ok
# Гильбертов многочлен кривой наклона c/b = 3/7
P = PolynomialRing(QQ, ['y%d' % i for i in range(9)] + ['s', 't'], order='degrevlex')
ys = P.gens()[:9]; s_, t_ = P.gens()[9:]
kq = QQ(3)/7
forms = [s_ + (i + j*kq)*t_ for i in (-1, 0, 1) for j in (-1, 0, 1)]
# исключаем s,t: y_k^2 = L_k(s,t); 7 линейных соотношений между y^2
M = matrix(QQ, [[f.coefficient({s_: 1}), f.coefficient({t_: 1})] for f in forms])
Kr = M.left_kernel().basis()
Q = PolynomialRing(QQ, ['y%d' % i for i in range(9)])
yq = Q.gens()
I7 = Q.ideal([sum(QQ(vv[i])*yq[i]**2 for i in range(9)) for vv in Kr])
out['slope_curve'] = dict(n_quadrics=len(Kr), dim=int(I7.dimension()), hilbert_poly=str(I7.hilbert_polynomial()))
print(json.dumps(out, indent=1, default=str))
json.dump(out, open('my_sage_checks.json', 'w'), indent=1, default=str)
