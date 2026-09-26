# Claude, 26.09.2026. Проверка результатов dem_rigorous.py по сохранённым JSON (не импортирует dem_rigorous.py).
# 1) Модель: для сохранённых Emin, M, T — тождество в Q(z): (4x³+b2x²+2b4x+b6)|_{x=M(z)} · ∏_T(1+λz) — квадрат в Q(z)
#    (значит, каждое z с ∏_T(1+λz) = □ даёт рациональную точку Emin с x = M(z)).
# 2) Сертификат Безу: формы = формула удвоения Emin, тождества A_j F + B_j G = R t^j точные, K и L пересчитаны.
# 3) Точное целочисленное неравенство для M0 (без логарифмов и плавающей точки):
#    H_k^{3(2M0+1)} > K^{2M0+1} · Kmax^{2·4^k} · c^{6·4^k},  H_k = H(x(2^k P0)), c = max строковых сумм M и adj M.
#    Оно равносильно (2M0+1)·(log H_k − log K/3)/4^k > (2/3)log Kmax + 2 log c, а левая часть ≤ (2M0+1)·ĥ(P0).
# 4) Случай B другим способом: закон группы Sage над функциональным полем Q(E), норма z(P)+z(P+t) в Q(x), корни.
# 5) Перебор случая A заново по сохранённым M, P0 (кручение пересчитывается), проверка всех восьми клеток.
from sage.all import *
import json, sys, time


def say(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def cells(r, s):
    base = [s, r, s - r, s + r]
    return sorted(set(base + [-c for c in base]))


def all_sq(z, S):
    for l in S:
        v = 1 + l * z
        if v < 0 or not v.is_square():
            return False
    return True


def Hproj(x):
    if x is None:
        return ZZ(1)
    x = QQ(x); return max(abs(x.numerator()), abs(x.denominator()))


def is_square_ratfun(f):
    """f ∈ Q(z): квадрат ли (все кратности чётны, старший коэффициент — квадрат)."""
    num, den = f.numerator(), f.denominator()
    for p in (num, den):
        fa = p.factor()
        if any(e % 2 for _, e in fa):
            return False
    u = num.leading_coefficient() / den.leading_coefficient()
    return u > 0 and QQ(u).is_square()


def verify_one(o, S):
    T = [ZZ(c) for c in o['T']]
    Emin = EllipticCurve(QQ, [QQ(v) for v in o['Emin']])
    assert Emin.is_minimal() and all(v in ZZ for v in Emin.a_invariants())
    M = matrix(ZZ, 2, 2, [ZZ(v) for v in o['M']]); a, b, c, d = M.list(); assert M.det() != 0
    b2, b4, b6, b8 = Emin.b_invariants()
    # 1) модель
    Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); Fz = Rz.fraction_field()
    X = Fz((a * z + b) / (c * z + d))
    f = 4 * X**3 + b2 * X**2 + 2 * b4 * X + b6
    assert is_square_ratfun(f * Rz.prod(1 + l * z for l in T)), 'модель z ↦ x_min не подтверждена'
    # 2) сертификат
    cert = o['cert']
    Rt = PolynomialRing(ZZ, 't'); t = Rt.gen()
    F = Rt([ZZ(v) for v in cert['F']]); G = Rt([ZZ(v) for v in cert['G']]); R = ZZ(cert['R']); dd = cert['d']
    assert dd == 4 and R != 0
    assert F == t**4 - b4 * t**2 - 2 * b6 * t - b8 and G == 4 * t**3 + b2 * t**2 + 2 * b4 * t + b6
    norms = []
    for w in cert['witnesses']:
        A = Rt([ZZ(v) for v in w['A']]); B = Rt([ZZ(v) for v in w['B']])
        assert A.degree() <= dd - 1 and B.degree() <= dd - 1
        assert A * F + B * G == R * t**w['j']
        norms.append(sum(abs(v) for v in A.list() + B.list()))
    assert sorted(w['j'] for w in cert['witnesses']) == [0, 2 * dd - 1]
    K = max(norms); L = max(sum(abs(v) for v in F.list()), sum(abs(v) for v in G.list()))
    assert K == ZZ(cert['K']) and L == ZZ(cert['L'])
    Kmax = max(K, L)
    # формула удвоения на точках
    x0, y0 = [QQ(v) for v in o['P0'].strip('()').split(',')]
    P0 = Emin(x0, y0); assert P0.order() == oo
    Q = P0
    for _ in range(3):
        assert (2 * Q)[0] == F(Q[0]) / G(Q[0]); Q = 2 * Q
    # 3) точное неравенство для M0
    k = 7; Q = P0
    for _ in range(k):
        Q = 2 * Q
    Hk = Hproj(None if Q.is_zero() else Q[0])
    cint = max(abs(a) + abs(b), abs(c) + abs(d), abs(d) + abs(b), abs(c) + abs(a))
    M0 = int(o['M0'])
    lhs = Hk**(3 * (2 * M0 + 1))
    rhs = K**(2 * M0 + 1) * Kmax**(2 * 4**k) * cint**(6 * 4**k)
    exact_ok = lhs > rhs
    # 4) случай B через функциональное поле
    tors = Emin.torsion_points()
    assert len(tors) == o['torsion']['n']
    Kx = PolynomialRing(QQ, 'x').fraction_field(); xg = Kx.gen()
    a1, a2, a3, a4, a6 = Emin.a_invariants()
    RY = PolynomialRing(Kx, 'Y'); Y = RY.gen()
    K2 = Kx.extension(Y**2 + a1 * xg * Y + a3 * Y - (xg**3 + a2 * xg**2 + a4 * xg + a6), 'y')
    EK = Emin.change_ring(K2); Pg = EK(K2(xg), K2.gen())
    def zf(xx):
        return (d * xx - b) / (-c * xx + a)
    zB = set(); zero_norm = False
    for tt in tors:
        if tt.is_zero():
            continue
        Pt = Pg + EK(tt)
        s_ = zf(K2(xg)) + zf(Pt[0])
        nrm = s_.norm()
        if nrm == 0:
            zero_norm = True; continue
        for rt, _ in nrm.numerator().roots(QQ):
            if -c * rt + a != 0:
                zB.add(QQ(zf(rt)))
    zB_main = set(QQ(v) for cb in o['caseB'] for v in cb.get('z', []))
    zB_nontriv = set(v for v in zB if v != 0)
    ztors = set()
    for tt in tors:
        if tt.is_zero():
            ztors.add(QQ(d) / QQ(-c)); continue
        if -c * tt[0] + a != 0:
            ztors.add(QQ((d * tt[0] - b) / (-c * tt[0] + a)))
    zBm_nontriv = set(v for v in zB_main if v != 0)
    # 5) перебор A заново
    sols = []; ncand = 0
    Qp = Emin(0); Qm = Emin(0)
    for n in range(0, M0 + 1):
        for Qn in ([Qp] if n == 0 else [Qp, Qm]):
            for tt in tors:
                R_ = Qn + tt
                xx = None if R_.is_zero() else R_[0]
                if xx is None:
                    num, den = QQ(d), QQ(-c)
                else:
                    num, den = d * xx - b, -c * xx + a
                if den == 0:
                    continue
                zz = num / den; ncand += 1
                if zz != 0 and (all_sq(zz, S) or all_sq(-zz, S)):
                    sols.append(str(zz))
        Qp = Qp + P0; Qm = Qm - P0
    for zz in zB | zB_main:
        if zz != 0 and (all_sq(zz, S) or all_sq(-zz, S)):
            sols.append('B:' + str(zz))
    return dict(T=[int(v) for v in T], model_identity=True, cert_ok=True, exact_M0_inequality=bool(exact_ok),
                digits_lhs=int(lhs.ndigits()), caseB_norm_zero=zero_norm,
                caseB_z_funcfield=sorted(str(v) for v in zB_nontriv), caseB_z_main=sorted(str(v) for v in zBm_nontriv),
                caseB_funcfield_subset_main=zB_nontriv <= zBm_nontriv, caseB_equal=(zB_nontriv == zBm_nontriv),
                caseB_main_extras_are_torsion_z=(zBm_nontriv - zB_nontriv) <= ztors,
                recount_A=ncand, nondeg=sols, closes=bool(exact_ok) and not zero_norm and not sols)


if __name__ == '__main__':
    out = {}
    for sl in sys.argv[1:]:
        r, s = map(int, sl.split('/'))
        D = json.load(open(f'/home/kep/magicKube/rigorous_census/chunk1/rig_{r}_{s}.json'))
        S = cells(r, s); assert D['cells'] == [int(v) for v in S]
        res = []
        for o in D['results']:
            if not o.get('ok'):
                continue
            v = verify_one(o, S)
            say(sl, v['T'], 'модель+сертификат ок; точное нер-во M0:', v['exact_M0_inequality'], f"({v['digits_lhs']} цифр)",
                'B: функц.поле ⊆ основной:', v['caseB_funcfield_subset_main'], 'лишние = z(кручения):', v['caseB_main_extras_are_torsion_z'], 'нетрив. B (функц.поле):', v['caseB_z_funcfield'],
                'перебор A:', v['recount_A'], 'невырожд.:', v['nondeg'], 'закрывает:', v['closes'])
            res.append(v)
        out[sl] = dict(results=res, closed=any(v['closes'] for v in res))
        say('ИТОГ', sl, 'закрыт (перепроверка):', out[sl]['closed'])
    json.dump(out, open('/home/kep/magicKube/rigorous_census/chunk1/verify_results.json', 'w'), ensure_ascii=False, indent=1)
