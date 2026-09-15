# Fable 15.09. Символьная проверка над ℚ(k): для пары (T_a,T_b) одного класса с |T_aΔT_b| ∈ {3,4}
# полюсы u_+ = x(Q_a + Q_b) — две сопряжённые точки S₁, S₂ третьей кривой E₃ (над z₁,z₂ — корнями x_a(z) = x_b(z)),
# и S₁ + S₂ ∈ E₃[2] (⇒ сумма полюсов 2(S₁+S₂) = O, спаривание ⟨Q_a,Q_b⟩ = O(1)).
# Считаем над K = ℚ(k)(√disc): точки S₁, S₂ и 2(S₁+S₂).
from sage.all import *
import sys, itertools
sys.path.insert(0, '/home/kep/magicKube/isogeny_mult/fable')
from classes import CLASSES, BIG

Kk = FunctionField(QQ, 'k'); k = Kk.gen()
def lam(g): return g[0] * 1 + g[1] * k            # s = 1, r = k
def weierK(T):
    c0 = T[0]; others = T[1:]; R = PolynomialRing(Kk, 'w'); w = R.gen()
    g = R.prod((c0 - c) * w + c * c0 for c in others)
    if len(T) == 3: g = c0 * w * g
    L = g.leading_coefficient(); h = (g(w / L) * L**2).monic(); co = h.list()
    return EllipticCurve(Kk, [0, co[2], 0, co[1], co[0]]), L, c0
def cellsK(S): return [lam(g) for g in S if g != (0, 0)]

def check_pair(Sa, Sb, verbose=True):
    Ta, Tb = cellsK(Sa), cellsK(Sb); D = [c for c in Ta if c not in Tb] + [c for c in Tb if c not in Ta]
    Ea, La, ca = weierK(Ta); Eb, Lb, cb = weierK(Tb); E3, L3, c3 = weierK(D)
    assert Ea.is_isomorphic(Eb)
    Es = Ea.short_weierstrass_model(); ia = Ea.isomorphism_to(Es); ib = Eb.isomorphism_to(Es)
    R = PolynomialRing(Kk, 'z'); z = R.gen(); Fz = R.fraction_field()
    def xs(iso, L, c0):
        u, r, s_, t_ = iso.tuple(); w = 1 / (z + 1 / c0); X = L * w
        return (X - r) / u**2
    xa, xb = xs(ia, La, ca), xs(ib, Lb, cb)
    q = (xa - xb).numerator(); q = R(q)
    f3z = prod(1 + c * z for c in D)
    if q.degree() == 1:
        # второй корень — z = ∞ (обе T трёхклеточные, Q_a(∞), Q_b(∞) — 2-кручение): полюсы u_+ над ∞ — пара ±, сумма O;
        # конечный корень z₁ обязан быть точкой ветвления E₃ (S₁ ∈ E₃[2]), иначе сумма полюсов 2S₁ ≠ O
        # (численно, 19/60: полюсы u_+ = пара точек над z₁ = −1/λ, λ ∈ T_a∩T_b, + пара над ∞; обе пары вида ±R ⇒ сумма O)
        # (численно, 19/60 и 204/247, пары (4,3) и (4,4) с c_a = c_b ∈ I, полюсы — пары ±R над z₁ = −1/λ (λ ∈ I) и над общим началом z = −1/c_a)
        z1 = -q[0] / q[1]; I = [c for c in Ta if c in Tb]
        ok = any(z1 == -1 / c for c in I) and f3z(z1) != 0 and ((len(Ta) == 3 and len(Tb) == 3 and len(D) == 4) or ca == cb)
        if verbose: print('  ', Ta, Tb, 'Δ =', D, ' deg q = 1: z₁ = −1/λ, λ ∈ T_a∩T_b, второй полюс над ∞ или над общим O (пары ±R):', ok, flush=True)
        return ok
    assert q.degree() == 2, q.degree()
    disc = q.discriminant()
    # поле ℚ(k)(√disc)
    Rs = PolynomialRing(Kk, 's'); s = Rs.gen()
    sq = disc.sqrt() if disc.is_square() else None
    if sq is not None:
        Kbig = Kk; roots = [(-q[1] + sq) / (2 * q[2]), (-q[1] - sq) / (2 * q[2])]
    else:
        Kbig = Kk.extension(s**2 - disc, 'sd'); sd = Kbig.gen()
        roots = [(-q[1] + sd) / (2 * q[2]), (-q[1] - sd) / (2 * q[2])]
    E3K = E3.change_ring(Kbig)
    pts = []
    for zr in roots:
        f3 = prod(1 + c * zr for c in D)      # y3² = f3
        # y3 из условия Q_a = −Q_b: y_a y_b = −(...)… проще: точка над zr с y3 — оба знака; выбираем по u_+ полюсу: проверяем оба
        w = 1 / (zr + 1 / c3); X = L3 * w
        Y2 = X**3 + E3.a2() * X**2 + E3.a4() * X + E3.a6()
        pts.append((X, Y2, zr))
    # Y = L3 c3 w² y3; знак y3 для полюса u_+: Q_a + Q_b = O ⇔ y_a'/... Условие: (y_b − y_a)/(x_b − x_a) → ∞ при y_a ≠ y_b, т.е. y_a = −y_b (в E_s-координатах).
    # y_a' y_b' = (u_a u_b)^{-3} Y_a Y_b = (u_a u_b)^{-3} (La ca wa²)(Lb cb wb²) ∏_{T_a∩T_b}(1+λz) · y3.  Нужен знак y3 с y_a' = −y_b'.
    ua = ia.tuple()[0]; ub = ib.tuple()[0]
    I = [c for c in Ta if c in Tb]
    res = []; nbranch = 0
    for X, Y2, zr in pts:
        if any(zr == -1 / c for c in I):
            nbranch += 1; continue            # корень — общая точка ветвления C_a, C_b: полюсы u_+ над ним — пара ±R, сумма O
        wa = 1 / (zr + 1 / ca); wb = 1 / (zr + 1 / cb)
        ya2 = (La * ca * wa**2)**2 * prod(1 + c * zr for c in Ta) / ua**6
        yb2 = (Lb * cb * wb**2)**2 * prod(1 + c * zr for c in Tb) / ub**6
        assert ya2 == yb2, 'x_a = x_b, но y_a² ≠ y_b²'
        # y_a' y_b' = coef · y3, y_a' = −y_b' ⇔ coef·y3 = −ya2 ⇔ y3 = −ya2/coef
        coef = (La * ca * wa**2) * (Lb * cb * wb**2) * prod(1 + c * zr for c in I) / (ua * ub)**3
        y3 = -ya2 / coef
        w3 = 1 / (zr + 1 / c3); Y = L3 * c3 * w3**2 * y3
        assert Y**2 == Y2, 'точка не на E₃'
        res.append(E3K(X, Y))
    Ssum = sum(res, E3K(0))
    two = 2 * Ssum
    if verbose: print('  ', Ta, Tb, 'Δ =', D, f' корней-ветвлений {nbranch};', 'сумма непарных полюсов S₁+S₂ ∈ E₃[2]:', two.is_zero(), '(= O:', Ssum.is_zero(), ')', flush=True)
    return two.is_zero()

if __name__ == '__main__':
    allok = True; n = 0
    for ci, c in enumerate(BIG):
        cl = CLASSES[c]
        cells = lambda S: frozenset(g for g in S if g != (0, 0))
        for i, j in itertools.combinations(range(len(cl)), 2):
            if len(cells(cl[i]) ^ cells(cl[j])) in (3, 4):
                ok = check_pair(cl[i], cl[j]); allok &= ok; n += 1
    print('пар рода 1 проверено:', n, 'все S₁+S₂ ∈ E₃[2]:', allok)
