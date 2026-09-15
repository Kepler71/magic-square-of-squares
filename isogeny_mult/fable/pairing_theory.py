# Fable 15.09. Спаривание ⟨Q_a,Q_b⟩ на общей кривой E для двух множителей T_a, T_b одного класса.
# Теория: x(Q_a ± Q_b) = u_±(R), u_± ∈ ℚ(E₃), E₃ = E_{T_aΔT_b} (третья кривая, |Δ| ∈ {3,4} ⇒ род 1), R = P_{T_aΔT_b}(z);
# u_−(R) = u_+(−R); h(u_+(R)) = 2ĥ(R) − ⟨R,S⟩ + O(1), S = сумма полюсов u_+ (точка E₃(ℚ)).
# ⇒ 4⟨Q_a,Q_b⟩ = h(u_+(R)) − h(u_+(−R)) + O(B_E) = −2⟨R,S⟩ + O(1).
# Здесь: строим u_+ символьно, находим S, проверяем ограниченность Λ(R)+2⟨R,S⟩ на многих R ∈ E₃(ℚ).
from sage.all import *
import sys, json
sys.path.insert(0, '/home/kep/magicKube/joint_sieve'); sys.path.insert(0, '/home/kep/magicKube/isogeny_mult/fable')
from classes import CLASSES, BIG, cells_of
from g1census import weier


def naive_h(x):
    x = QQ(x); return float(RR(max(abs(x.numerator()), abs(x.denominator()))).log())


class Fac:
    """модель weier для T; short Weierstrass-образ E' и изоморфизм; X,Y как функции z, y_T."""
    def __init__(self, T):
        self.T = sorted(ZZ(t) for t in T)
        E, L, c0 = weier([QQ(t) for t in self.T]); self.E, self.L, self.c0 = E, QQ(L), QQ(c0)
        a = E.a_invariants(); self.a2, self.a4, self.a6 = QQ(a[1]), QQ(a[3]), QQ(a[4])

    def XY(self, z, yT):
        """координаты P_T(z) на модели weier (yT = ∏√(1+λz), любая алгебра, где z, yT живут)."""
        w = 1 / (z + 1 / self.c0)
        return self.L * w, self.L * self.c0 * w**2 * yT

    def z_of_X(self, X): return -1 / self.c0 + self.L / X


def setup_pair(Ta, Tb):
    fa, fb = Fac(Ta), Fac(Tb)
    D = sorted(set(fa.T) ^ set(fb.T)); I = sorted(set(fa.T) & set(fb.T))
    assert len(D) in (1, 2, 3, 4)
    f3 = Fac(D) if len(D) >= 3 else None
    # общая кривая: короткая модель E' и изоморфизмы (x,y) ↦ (u²x + r, u³y) из E_a, E_b (a1=a3=0 ⇒ s=t=0)
    Es = fa.E.short_weierstrass_model()
    ia = fa.E.isomorphism_to(Es); ib = fb.E.isomorphism_to(Es)
    for iso in (ia, ib):
        u, r, s_, t_ = iso.tuple(); assert s_ == 0 and t_ == 0
    # функция u_+ на E₃: в ℚ(z)[y3], y3² = ∏_{λ∈D}(1+λz); y_a y_b = ∏_{λ∈I}(1+λz)·y3
    R = PolynomialRing(QQ, 'z'); z = R.gen(); Fz = R.fraction_field()
    S3 = PolynomialRing(Fz, 'y3'); y3 = S3.gen()
    f3pol = prod(1 + lam * z for lam in D); q3 = S3.quotient(y3**2 - f3pol, 'Y3'); Y3 = q3.gen()
    def img(fac, iso, ysq):
        u, r, s_, t_ = iso.tuple()
        Xw, Yw = fac.XY(Fz(z), 1)  # Yw = коэффициент при y_T
        # x' = (X − r)/u², y' = Y/u³   (Sage: iso.tuple()=(u,r,s,t) с x = u²x'+r)
        return (Xw - r) / u**2, Yw / u**3
    xa, ya = img(fa, ia, None); xb, yb = img(fb, ib, None)   # ya, yb — множители при y_Ta, y_Tb
    ya2 = ya**2 * prod(1 + lam * z for lam in fa.T); yb2 = yb**2 * prod(1 + lam * z for lam in fb.T)
    yayb = ya * yb * prod(1 + lam * z for lam in I)    # умножить на y3
    lam2 = (ya2 + yb2 - 2 * yayb * Y3) / (xa - xb)**2   # λ² для Q_a + Q_b  (y_b − y_a)²/(x_b − x_a)²
    uplus = lam2 - xa - xb            # short Weierstrass: x(P+Q) = λ² − x1 − x2
    uminus = (ya2 + yb2 + 2 * yayb * Y3) / (xa - xb)**2 - xa - xb
    A = (ya2 + yb2) / (xa - xb)**2 - xa - xb; B = -2 * yayb / (xa - xb)**2   # u_± = A ± B·y3
    return dict(fa=fa, fb=fb, f3=f3, D=D, I=I, Es=Es, ia=ia, ib=ib, A=A, B=B, z=z)


def eval_u(P, R3, sign=+1):
    """u_±(R) для R ∈ E₃(ℚ) (модель weier f3): z = z(R), y3 из Y-координаты."""
    f3 = P['f3']
    if R3.is_zero(): return None
    X, Y = R3[0], R3[1]
    if X == 0: return None
    zv = f3.z_of_X(X); w = X / f3.L; y3 = Y / (f3.L * f3.c0 * w**2)
    try:
        Av = P['A'].numerator()(zv) / P['A'].denominator()(zv); Bv = P['B'].numerator()(zv) / P['B'].denominator()(zv)
    except ZeroDivisionError:
        return None
    return Av + sign * Bv * y3


def pole_sum(P):
    """сумма полюсов u_+ как точка E₃(ℚ): через поле функций Sage."""
    f3 = P['f3']; z = P['z']
    # поле функций E₃ в координатах (z, y3): y3² = f3pol(z); полюсы u_+ = A + B y3
    K = FunctionField(QQ, 'z'); zz = K.gen(); Ry = PolynomialRing(K, 'yy'); yy = Ry.gen()
    f3pol = prod(1 + lam * zz for lam in P['D'])
    L = K.extension(yy**2 - f3pol, 'y3'); y3 = L.gen()
    An, Ad = P['A'].numerator(), P['A'].denominator(); Bn, Bd = P['B'].numerator(), P['B'].denominator()
    A = K(An)(zz) / K(Ad)(zz) if False else L(K(str(An)) / K(str(Ad)))
    Bf = L(K(str(Bn)) / K(str(Bd)))
    u = A + Bf * y3
    Dpol = u.divisor_of_poles(); Dzer = u.divisor_of_zeros()
    return u, Dpol, Dzer, L


if __name__ == '__main__':
    r, s = int(sys.argv[1]), int(sys.argv[2])
    ci = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    cl = CLASSES[BIG[ci]]
    Ts = [cells_of(S, r, s) for S in cl]
    print('класс', ci, 'кратность', len(cl), Ts)
    pairs = [(i, j) for i in range(len(Ts)) for j in range(i + 1, len(Ts)) if len(set(Ts[i]) ^ set(Ts[j])) in (3, 4)]
    print('пары с третьей кривой рода 1:', [(Ts[i], Ts[j], sorted(set(Ts[i]) ^ set(Ts[j]))) for i, j in pairs])
    for i, j in pairs[:2]:
        P = setup_pair(Ts[i], Ts[j])
        u, Dp, Dz, L = pole_sum(P)
        print('  пара', Ts[i], Ts[j], 'третья', P['D'], 'deg полюсов', Dp.degree(), 'полюсы:', Dp, '\n   нули:', Dz)
