# Fable, 25.09.2026. Собственная реализация (без импорта кода Codex и старого кода Fable).
# Кривые y_T^2 = ∏_{c∈T}(1+cz) (|T|∈{3,4}), общая модель E (a1=a3=0), мёбиусовы матрицы x(Q_T)=M(z),
# функции u± = x(Q_a ± Q_b) на вспомогательной кривой D = T_aΔT_b (коника при |D|=2, эллиптическая при |D|∈{3,4}).
from sage.all import *

def cells(p, q):
    p, q = ZZ(p), ZZ(q); assert 0 < p < q
    return sorted([-(q+p), -q, -(q-p), -p, p, q-p, q, q+p])

def target_model(E):
    """модель с a1=a3=0: минимальная, если она такова, иначе [0,b2,0,8b4,16b6] минимальной."""
    Em = E.minimal_model()
    if Em.a1() == 0 and Em.a3() == 0: return Em
    b2, b4, b6, b8 = Em.b_invariants()
    return EllipticCurve(QQ, [0, b2, 0, 8*b4, 16*b6])

def _Fz():
    return PolynomialRing(QQ, 'z').fraction_field()

class Fac:
    """E_T: Y^2 = h(X), X = L w, Y = L c0 w^2 y_T, w = 1/(z + 1/c0), c0 = min T. Тождество проверяется в Q(z)."""
    def __init__(self, T):
        self.T = sorted(QQ(c) for c in T); assert len(self.T) in (3, 4) and len(set(self.T)) == len(self.T)
        self.c0 = self.T[0]
        R = PolynomialRing(QQ, 'w'); w = R.gen()
        g = prod((self.c0 - c)*w + c*self.c0 for c in self.T[1:])
        if len(self.T) == 3: g *= self.c0*w
        self.L = g.leading_coefficient()
        h = R(self.L**2 * g(w/self.L)); assert h.degree() == 3 and h.is_monic()
        self.h = h
        self.E = EllipticCurve(QQ, [0, h[2], 0, h[1], h[0]])
        z = _Fz().gen(); X, Yc = self.XY(z)
        assert Yc**2 * prod(1 + c*z for c in self.T) == h(X)
    def XY(self, z):
        w = 1/(z + 1/self.c0)
        return self.L*w, self.L*self.c0*w**2

class Image:
    """Q_T(z) ∈ Et: x = (X−r)/u^2 = M(z) (M целая примитивная), y = ycoef(z)·y_T."""
    def __init__(self, fac, Et):
        self.fac, self.Et = fac, Et
        u, r, s, t = fac.E.isomorphism_to(Et).tuple(); assert s == 0 and t == 0
        self.u, self.r = u, r
        c0, L = fac.c0, fac.L
        M = matrix(QQ, [[-r*c0, L*c0 - r], [u**2*c0, u**2]])
        M *= lcm(c.denominator() for c in M.list()); M /= gcd(ZZ(c) for c in M.list())
        self.M = M.change_ring(ZZ); assert self.M.det() != 0
        z = _Fz().gen(); x, yc = self.xy(z)
        assert x == (M[0,0]*z + M[0,1])/(M[1,0]*z + M[1,1])
        assert yc**2 * prod(1 + c*z for c in fac.T) == x**3 + Et.a2()*x**2 + Et.a4()*x + Et.a6()
    def xy(self, z):
        X, Yc = self.fac.XY(z); return (X - self.r)/self.u**2, Yc/self.u**3
    def z_of_x(self, x):
        """обратный мёбиус; x=None ↔ ∞; ответ None ↔ z=∞."""
        M = self.M
        if x is None: return None if M[1,0] == 0 else -M[1,1]/M[1,0]
        den = -M[1,0]*x + M[0,0]
        return None if den == 0 else (M[1,1]*x - M[0,1])/den
    def mobius_norm(self):
        return ZZ(2*max(abs(c) for c in self.M.list()))

def pair_AB(ima, imb):
    """x(Q_a ± Q_b) = A ± B·y_D в Q(z); y_{T_a} y_{T_b} = ∏_{I}(1+cz)·y_D, D = T_aΔT_b, I = T_a∩T_b."""
    assert ima.Et == imb.Et
    z = _Fz().gen(); Et = ima.Et
    xa, ya = ima.xy(z); xb, yb = imb.xy(z)
    Ta, Tb = set(ima.fac.T), set(imb.fac.T); I = sorted(Ta & Tb); D = sorted(Ta ^ Tb)
    ya2 = ya**2 * prod(1 + c*z for c in Ta); yb2 = yb**2 * prod(1 + c*z for c in Tb)
    A = (ya2 + yb2)/(xa - xb)**2 - Et.a2() - xa - xb
    B = -2*ya*yb*prod(1 + c*z for c in I)/(xa - xb)**2
    return A, B, D, I

def exact_pair_check(ima, imb, A, B, D, zval):
    """точный контроль формулы сложения: точки Q_a, Q_b над числовым полем Q(√a,√b), сложение средствами Sage."""
    z = QQ(zval)
    a = prod(1 + c*z for c in ima.fac.T); b = prod(1 + c*z for c in imb.fac.T)
    R = PolynomialRing(QQ, 't'); t = R.gen()
    Ka = QQ if a.is_square() else NumberField(t**2 - a, 'sa')
    Rb = PolynomialRing(Ka, 's'); s = Rb.gen()
    bb = Ka(b)
    if bb.is_square(): Kb = Ka
    else: Kb = Ka.extension(s**2 - bb, 'sb')
    sa = Kb(a).sqrt(); sb = Kb(b).sqrt()
    EK = ima.Et.change_ring(Kb)
    xa, ya = ima.xy(z); xb, yb = imb.xy(z)
    Pa = EK([Kb(xa), Kb(ya)*sa]); Pb = EK([Kb(xb), Kb(yb)*sb])
    d = prod(1 + c*z for c in D); I = sorted(set(ima.fac.T) & set(imb.fac.T))
    yD = sa*sb/prod(1 + c*z for c in I); assert yD**2 == d
    Az = A.numerator()(z)/A.denominator()(z); Bz = B.numerator()(z)/B.denominator()(z)
    got = {(Pa + Pb)[0], (Pa - Pb)[0]}; want = {Az + Bz*yD, Az - Bz*yD}
    return got == want

def rational_function_eval(f, z):
    """f ∈ Q(z) → значение в элементе z произвольного поля (через списки коэффициентов)."""
    num, den = f.numerator(), f.denominator()
    return sum(QQ(c)*z**i for i, c in enumerate(num.list())) / sum(QQ(c)*z**i for i, c in enumerate(den.list()))

class Aux:
    """вспомогательная кривая y_D^2 = ∏_{c∈D}(1+cz), |D|∈{3,4}: модель E3 (a1=a3=0) и u± ∈ Q(E3) = Q(x)[y]/(y^2−f3)."""
    def __init__(self, A, B, D):
        self.D = sorted(QQ(c) for c in D)
        self.fac = Fac(self.D); self.E3 = target_model(self.fac.E)
        u3, r3, s, t = self.fac.E.isomorphism_to(self.E3).tuple(); assert s == 0 and t == 0
        K = FunctionField(QQ, 'x'); x = K.gen(); Ry = PolynomialRing(K, 'y'); y = Ry.gen()
        f3 = x**3 + self.E3.a2()*x**2 + self.E3.a4()*x + self.E3.a6()
        Lf = K.extension(y**2 - f3, 'y'); y = Lf.gen()
        self.K, self.Lf, self.x, self.y, self.f3 = K, Lf, x, y, f3
        X3 = u3**2*x + r3; L3, d0 = self.fac.L, self.fac.c0
        self.z = L3/X3 - 1/d0
        self.yD = u3**3*L3*y/(d0*X3**2)
        assert self.yD**2 == prod(1 + c*self.z for c in self.D)
        Az = rational_function_eval(A, self.z); Bz = rational_function_eval(B, self.z)
        self.up = Az + Bz*self.yD; self.um = Az - Bz*self.yD
    def two_torsion(self):
        return [P for P in self.E3.torsion_points() if not P.is_zero() and (2*P).is_zero()]
    def basis_3O_plus_S(self, S):
        """базис L(3O+S): S=None ↔ O: (1,x,y,x^2); S=(xs,ys): (1,x,y,(y+ys)/(x−xs))."""
        x, y = self.x, self.y
        if S is None: return [self.Lf(1), self.Lf(x), y, self.Lf(x**2)]
        xs, ys = QQ(S[0]), QQ(S[1])
        return [self.Lf(1), self.Lf(x), y, (y + ys)/(x - xs)]

def coeffs_xy(elt):
    """элемент a(x)+b(x)y поля Q(x)[y]/(y^2−f3) → (a, b) ∈ Q(x)."""
    l = list(elt.list()); assert len(l) <= 2
    l += [0]*(2 - len(l)); return l[0], l[1]

def solve_linear_forms(u, basis):
    """ищет (n_i, h_i): Σ n_i φ_i = u·Σ h_i φ_i в поле функций; возвращает список решений (базис ядра) как пары векторов."""
    cols = []
    for phi in basis: cols.append(coeffs_xy(phi))
    for phi in basis: a, b = coeffs_xy(u*phi); cols.append((-a, -b))
    den = lcm([c.denominator() for col in cols for c in col])
    polys = [[(c*den).numerator() for c in col] for col in cols]
    deg = max(c.degree() for col in polys for c in col)
    rows = []
    for part in (0, 1):
        for k in range(deg + 1):
            rows.append([QQ(col[part][k]) for col in polys])
    Mt = matrix(QQ, rows)
    ker = Mt.right_kernel().basis()
    return [(vector(QQ, v[:len(basis)]), vector(QQ, v[len(basis):])) for v in ker]

def vanishing_forms(basis, degree):
    """все однородные формы данной степени от координат v_i = φ_i, обращающиеся в нуль на кривой."""
    nvar = len(basis)
    V = PolynomialRing(QQ, ['v%d' % i for i in range(nvar)]); vs = V.gens()
    mons = [prod(v**int(e) for v, e in zip(vs, ee)) for ee in IntegerVectors(degree, nvar)]
    vals = [prod(phi**int(e) for phi, e in zip(basis, ee)) for ee in IntegerVectors(degree, nvar)]
    cols = [coeffs_xy(val) for val in vals]
    den = lcm([c.denominator() for col in cols for c in col])
    polys = [[(c*den).numerator() for c in col] for col in cols]
    deg = max(c.degree() for col in polys for c in col)
    rows = [[QQ(col[part][k]) for col in polys] for part in (0, 1) for k in range(deg + 1)]
    ker = matrix(QQ, rows).right_kernel().basis()
    out = []
    for v in ker:
        f = sum(QQ(c)*m for c, m in zip(v, mons))
        f *= lcm(c.denominator() for c in f.coefficients()); f /= gcd(ZZ(c) for c in f.coefficients())
        out.append(V(f))
    return V, out
