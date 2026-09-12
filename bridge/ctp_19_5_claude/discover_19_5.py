"""Независимое построение сертификата Кассельса-Тейта для G1 (m,n)=(19,5).

Пишется с нуля: НЕ загружает descent/ctp_quartic.sage и не использует чужие числа.
Sage нужен только на этапе ОТКРЫТИЯ (решение коники, ell2cover как источник кандидатов).
Проверка сертификата делается отдельным скриптом на стандартном Python.

Ранг эллиптической кривой НЕ вычисляется и НЕ используется.
"""
from sage.all import *
import json, time, itertools, random
from pathlib import Path

OUT = Path('/home/kep/magicKube/bridge/ctp_19_5_claude')
OUT.mkdir(parents=True, exist_ok=True)
random.seed(1905)
set_random_seed(1905)

# ---------------------------------------------------------------- инварианты квартик
def quartI(g):
    a, b, c, d, e = g
    return 12*a*e - 3*b*d + c*c

def quartJ(g):
    a, b, c, d, e = g
    return 72*a*c*e - 27*a*d*d - 27*b*b*e + 9*b*c*d - 2*c**3

def hessian(g):
    a, b, c, d, e = g
    return [3*b*b - 8*a*c,
            4*(b*c - 6*a*d),
            2*(2*c*c - 24*a*e - 3*b*d),
            4*(c*d - 6*b*e),
            3*d*d - 8*c*e]

def ev(coeffs, x, z):
    """coeffs = [c_0..c_k] при x^k, x^(k-1) z, ..., z^k"""
    k = len(coeffs) - 1
    return sum(coeffs[i]*x**(k-i)*z**i for i in range(k+1))

# ---------------------------------------------------------------- исходная кривая
m, n = 19, 5
s = QQ(m**2 + n**2)/2
b_par = s*m**2*n**2
E_roots = [-b_par, -s*m**4, -s*n**4]          # e1, e2, e3 в порядке из задачи
print('m,n =', m, n, ' s =', s, ' b =', b_par)
print('E_roots =', E_roots)
assert s == 193 and b_par == 1741825

# требуемый квадратный класс: проверяем тождества X-e_i = c_i * F_j напрямую
Rt = PolynomialRing(QQ, 't'); t = Rt.gen()
F0 = m**2 + n**2*t**2
F4 = s*(1 + t**2)
F8 = n**2 + m**2*t**2
X = b_par*t**2
assert X - E_roots[0] == QQ(m*n)**2 * F4 / s * s / s * s   # см. ниже, аккуратно:
# X - e1 = b(1+t^2) = s m^2 n^2 (1+t^2) = (mn)^2 * [s(1+t^2)] = (mn)^2 F4
assert X - E_roots[0] == QQ(m*n)**2 * F4
assert X - E_roots[1] == s*m**2 * F0
assert X - E_roots[2] == s*n**2 * F8
print('OK: X-e1=(mn)^2 F4, X-e2=s m^2 F0, X-e3=s n^2 F8')
assert not QQ(s).is_square()
DELTA = [QQ(1), QQ(s), QQ(s)]   # squarefree представители

# ---------------------------------------------------------------- модель y^2=x^3+Ax+B
sh = -sum(E_roots)/3
assert sh in ZZ
r_shift = [ZZ(e + sh) for e in E_roots]
gg = gcd(r_shift)
u2 = ZZ(1)
for p, k in factor(gg):
    u2 *= p**(2*(k//2))
u = ZZ(u2).sqrt()
rr = [ZZ(v/u2) for v in r_shift]
assert sum(rr) == 0
A = rr[0]*rr[1] + rr[0]*rr[2] + rr[1]*rr[2]
B = -rr[0]*rr[1]*rr[2]
print('shift =', sh, ' u^2 =', u2, ' roots =', rr)
print('A,B =', A, B)
M = EllipticCurve([0, 0, 0, A, B])
print('M =', M, ' minimal?', M.is_minimal() if hasattr(M, 'is_minimal') else '?')
print('M.minimal_model ainvs =', M.minimal_model().ainvs())

I = -48*A
J = -1728*B
assert I == M.c4() and J == 2*M.c6()
phis = [-12*ri for ri in rr]
Rx = PolynomialRing(QQ, 'X'); Xv = Rx.gen()
for ph in phis:
    assert ph**3 - 3*I*ph + J == 0
DISC = 16*(4*I**3 - J**2)/27
assert DISC != 0
print('I =', I); print('J =', J); print('Disc =', DISC)

# ---------------------------------------------------------------- L = Q^3 (полное 2-кручение)
Dl = [prod(rr[i] - rr[j] for j in range(3) if j != i) for i in range(3)]
# коэффициенты лагранжевых базисных многочленов L_i(X) = (X^2 + p1_i X + p0_i)/D_i
Lc2 = [QQ(1)/Dl[i] for i in range(3)]
Lc1 = [QQ(-sum(rr[j] for j in range(3) if j != i))/Dl[i] for i in range(3)]
Lc0 = [QQ(prod(rr[j] for j in range(3) if j != i))/Dl[i] for i in range(3)]
# контроль интерполяции
for k in range(3):
    w = [QQ(rr[i]**k) for i in range(3)]
    c2 = sum(w[i]*Lc2[i] for i in range(3))
    c1 = sum(w[i]*Lc1[i] for i in range(3))
    c0 = sum(w[i]*Lc0[i] for i in range(3))
    poly = c2*Xv**2 + c1*Xv + c0
    assert poly == Xv**k % prod(Xv - ri for ri in rr)
print('OK: интерполяция в L корректна')

def sqfree(q):
    q = QQ(q)
    return QQ(q.numerator()*q.denominator()).squarefree_part()

def sameclass(a, c):
    a, c = QQ(a), QQ(c)
    return a != 0 and c != 0 and (a*c).is_square()

# ---------------------------------------------------------------- квартика из класса delta
def quartic_from_delta(dd, tag=''):
    """dd = [d1,d2,d3] рациональные, класс в (Q*/Q*^2)^3 = L*/L*^2 с квадратной нормой.
    Возвращает бинарную квартику g с I(g)=I, J(g)=J и z(g) в классе dd."""
    dd = [QQ(x) for x in dd]
    assert (dd[0]*dd[1]*dd[2]).is_square(), 'норма не квадрат'
    cc = [dd[i]*Lc2[i] for i in range(3)]      # коника sum c_i T_i^2 = 0
    den = lcm([QQ(c).denominator() for c in cc])
    cc_int = [ZZ(c*den) for c in cc]
    gI = gcd(cc_int)
    cc_int = [ZZ(c/gI) for c in cc_int]
    cc_red = [ZZ(QQ(c).squarefree_part()) for c in cc_int]   # квадратные множители убираем
    # T_i -> T_i / sqrt(c_i/c_red_i):  масштаб координат, запомним
    scal = [QQ(cc_int[i]/cc_red[i]).sqrt() for i in range(3)]
    for i in range(3):
        assert cc_red[i]*scal[i]**2 == cc_int[i]
    con = Conic(QQ, cc_red)
    ok, pt = con.has_rational_point(point=True)
    assert ok, 'коника без точки: класс не в Селмере?'
    par = con.parametrization(pt)[0]
    defs = par.defining_polynomials()
    Rxz = PolynomialRing(QQ, ['xx', 'zz'])
    xx, zz = Rxz.gens()
    Tred = [Rxz(p) for p in defs]
    # проверяем, что параметризация лежит на конике
    assert sum(cc_red[i]*Tred[i]**2 for i in range(3)) == 0
    T = [Tred[i]/scal[i] for i in range(3)]
    assert sum(cc_int[i]*T[i]**2 for i in range(3)) == 0
    assert sum(cc[i]*T[i]**2 for i in range(3)) == 0        # [delta t^2]_2 = 0
    graw = -sum(dd[i]*T[i]**2*Lc1[i] for i in range(3))     # g = -[delta t^2]_1
    # чистим общий квадратный множитель содержания (это собственная эквивалентность при
    # одновременном учёте нормировки инвариантов ниже -- нормировка всё равно фиксирует масштаб)
    gc = [QQ(graw.coefficient({xx: 4-i, zz: i})) for i in range(5)]
    assert graw == sum(gc[i]*xx**(4-i)*zz**i for i in range(5))
    Ig, Jg = quartI(gc), quartJ(gc)
    assert Ig != 0 or Jg != 0
    lam2 = QQ(Jg*I)/QQ(J*Ig)
    gc = [QQ(x)/lam2 for x in gc]
    assert quartI(gc) == I and quartJ(gc) == J, (quartI(gc), I, quartJ(gc), J)
    gc = make_z_unit(gc)
    zv = z_inv(gc)
    for i in range(3):
        assert sameclass(zv[i], dd[i]), ('класс z(g) не совпал', tag, i, zv[i], dd[i])
    return gc

def z_inv(g):
    a, b, c, d, e = g
    return [(4*a*ph + 3*b*b - 8*a*c)/3 for ph in phis]

def act(g, mat):
    """g -> g(alpha x + beta z, gamma x + delta z), det=1 => I,J сохраняются"""
    Rxz = PolynomialRing(QQ, ['xx', 'zz'])
    xx, zz = Rxz.gens()
    al, be, ga, de = mat
    P = sum(g[i]*(al*xx+be*zz)**(4-i)*(ga*xx+de*zz)**i for i in range(5))
    return [QQ(P.coefficient({xx: 4-i, zz: i})) for i in range(5)]

def make_z_unit(g):
    if all(v != 0 for v in z_inv(g)):
        return g
    for trial in range(200):
        k = random.randint(-4, 4); l = random.randint(-4, 4)
        mat = (1, k, l, 1 + k*l)
        assert mat[0]*mat[3] - mat[1]*mat[2] == 1
        g2 = act(g, mat)
        if all(v != 0 for v in z_inv(g2)):
            assert quartI(g2) == quartI(g) and quartJ(g2) == quartJ(g)
            return g2
    raise RuntimeError('не удалось сделать z(g) единицей')

t0 = time.time()
g1 = quartic_from_delta(DELTA, 'target')
print('g1 =', g1, ' за', time.time()-t0, 'c')

record = {
    'm': m, 'n': n, 's': str(s), 'b': str(b_par),
    'E_roots': [str(x) for x in E_roots],
    'shift': str(sh), 'u2': str(u2),
    'model_roots': [str(x) for x in rr], 'A': str(A), 'B': str(B),
    'I': str(I), 'J': str(J), 'phi_roots': [str(x) for x in phis],
    'discriminant': str(DISC),
    'delta': [str(x) for x in DELTA],
    'g1': [str(x) for x in g1],
}
(OUT/'discover_19_5.json').write_text(json.dumps(record, indent=2)+'\n')
print('--- этап 1 записан ---', flush=True)
