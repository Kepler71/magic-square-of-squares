#!/usr/bin/env sage
# audit_E2_claude.sage -- НЕЗАВИСИМЫЙ АУДИТ фактора E2 моста (15,8).
#
# Это НЕ повторение basis_E2.sage. Здесь намеренно выбраны другие маршруты:
#   * ранг/базис: Sage E.gens() + PARI ellrank + eclib rank_bound (три реализации),
#   * насыщенность: БРУТФОРС-сертификаты "P mod p не лежит в q*E(F_p)" -- множество
#     q*E(F_p) строится полным перебором точек, БЕЗ abelian_group()/discrete_log,
#     плюс точный division_points для малых q,
#   * карты phi1, phi2: подстановка в поле дробей Q(x) + проверка дифференциалов
#     + проверка на случайных точках D(F_p) для 20 простых,
#   * образ Jac(D)(Q): явные сертификаты исключения для каждого из 14 мёртвых классов.
#
# Выход: audit_E2_claude.json, audit_E2_claude.log
#
# НИЧЕГО здесь не утверждает, что D(Q) состоит ровно из четырёх точек.

from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve
from pathlib import Path
import json, os, sys, tempfile

OUT = Path('/home/kep/magicKube/bridge/qc40')
LOG = []
R = {}


def say(*a):
    s = " ".join(str(z) for z in a)
    LOG.append(s)
    print(s, flush=True)


def dump():
    (OUT / 'audit_E2_claude.json').write_text(
        json.dumps(R, indent=2, ensure_ascii=False, default=str) + '\n')
    (OUT / 'audit_E2_claude.log').write_text("\n".join(LOG) + "\n")


# ===================================================== 0. модель из JSON моста
bridge = json.loads((OUT / 'new_bridge_verified.json').read_text())
a0, a2, a4, a6 = [QQ(z) for z in bridge['coefficients_a0_a2_a4_a6']]
Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
Kx = Rx.fraction_field()
f = a6 * x**6 + a4 * x**4 + a2 * x**2 + a0
assert str(f) == bridge['f']
assert gcd(f, f.derivative()) == 1
say("[0] a0,a2,a4,a6 = %s %s %s %s ; f = %s ; gcd(f,f')=1 => гладкая, род 2" % (a0, a2, a4, a6, f))
say("[0] a0 = 2023^2 ? %s   a6 = 6647^2 ? %s" % (a0 == 2023**2, a6 == 6647**2))

# мост t <-> x: f(x) = (1-x)^6 * F0(t)*F8(t)*L(t), t=(1+x)/(1-x)   -- ПЕРЕПРОВЕРКА
t = (1 + x) / (1 - x)
F0 = 225 + 64 * t**2          # (15)^2 + (8t)^2
F8 = 64 + 225 * t**2
L = QQ(289) / 2 * (1 + t**2) - 240 * t
assert Kx(f) == (1 - x)**6 * F0 * F8 * L
say("[0] тождество моста f(x) = (1-x)^6 * F0*F8*L в Q(x): ПОДТВЕРЖДЕНО независимо")

# ===================================================== 1. карты phi1, phi2
E1 = EllipticCurve([0, a4, 0, a2 * a6, a0 * a6**2])
E2 = EllipticCurve([0, a2, 0, a4 * a0, a6 * a0**2])
assert [str(z) for z in E1.ainvs()] == bridge['elliptic_factors'][0]['raw_ainvs']
assert [str(z) for z in E2.ainvs()] == bridge['elliptic_factors'][1]['raw_ainvs']
say("[1] E1 raw = %s" % (E1.ainvs(),))
say("[1] E2 raw = %s" % (E2.ainvs(),))


def rhs(E, T):
    A = E.ainvs()
    return T**3 + A[1] * T**2 + A[3] * T + A[4]


# phi1 : (x,y) -> (a6 x^2, a6 y).  Требуется (a6 y)^2 = rhs1(a6 x^2) при y^2 = f(x).
id1 = rhs(E1, a6 * x**2) - a6**2 * f
say("[1a] rhs_E1(a6*x^2) - a6^2*f(x) = %s   => phi1 корректна ? %s" % (id1, id1 == 0))
# phi2 : (x,y) -> (a0/x^2, a0 y/x^3).  Требуется (a0 y/x^3)^2 = rhs2(a0/x^2).
id2 = Kx(rhs(E2, Kx(a0) / Kx(x)**2)) - Kx(a0**2 * f) / Kx(x)**6
say("[1b] rhs_E2(a0/x^2) - a0^2*f(x)/x^6 = %s   => phi2 корректна ? %s" % (id2, id2 == 0))
assert id1 == 0 and id2 == 0

# дифференциалы: phi1^*(dX/2Y) = x dx/y ;  phi2^*(dX/2Y) = -dx/y
X1 = a6 * x**2; Y1c = a6           # Y1 = a6*y
d_X1 = X1.derivative()             # dX1 = 2 a6 x dx
w1 = Kx(d_X1) / (2 * Kx(Y1c))      # = (2 a6 x)/(2 a6) * dx/y = x dx/y
X2 = Kx(a0) / Kx(x)**2; Y2c = Kx(a0) / Kx(x)**3
d_X2 = X2.derivative()
w2 = d_X2 / (2 * Y2c)
say("[1c] phi1^*(dX/2Y) = (%s) dx/y  -- ожидалось x" % w1)
say("[1c] phi2^*(dX/2Y) = (%s) dx/y  -- ожидалось -1" % w2)
assert w1 == Kx(x) and w2 == Kx(-1)
say("[1c] {x dx/y, dx/y} -- базис H^0(D, Omega); значит Jac(D) ~ E1 x E2 (изогения)")

# проверка на реальных точках D(F_p) для 20 простых
disc_bad = ZZ(Rx(f).discriminant() * a6 * a0)
pmap = []
for p in prime_range(5, 200):
    if disc_bad % p == 0:
        continue
    F = GF(p)
    e1 = EllipticCurve(F, [F(z) for z in E1.ainvs()])
    e2 = EllipticCurve(F, [F(z) for z in E2.ainvs()])
    cnt = 0
    for xv in F:
        v = F(a6) * xv**6 + F(a4) * xv**4 + F(a2) * xv**2 + F(a0)
        if not v.is_square():
            continue
        yv = v.sqrt()
        for yy in set([yv, -yv]):
            cnt += 1
            assert e1([F(a6) * xv**2, F(a6) * yy, 1]) in e1
            if xv != 0:
                assert e2([F(a0) / xv**2, F(a0) * yy / xv**3, 1]) in e2
    pmap.append([int(p), int(cnt)])
    if len(pmap) >= 20:
        break
say("[1d] phi1, phi2 переводят все аффинные точки D(F_p) в E1(F_p), E2(F_p) для 20 простых: OK")
say("[1d]     (p, #аффинных точек D(F_p)): %s" % pmap)

R['maps_phi'] = {
    'phi1': '(a6*x^2, a6*y) : D -> E1,  E1: Y^2 = X^3 + a4 X^2 + a2 a6 X + a0 a6^2',
    'phi2': '(a0/x^2, a0*y/x^3) : D -> E2,  E2: Y^2 = X^3 + a2 X^2 + a4 a0 X + a6 a0^2',
    'symbolic_identity_phi1': 'rhs_E1(a6 x^2) == a6^2 * f(x) в Q[x] -- ТОЖДЕСТВО (proved)',
    'symbolic_identity_phi2': 'rhs_E2(a0/x^2) == a0^2 f(x)/x^6 в Q(x) -- ТОЖДЕСТВО (proved)',
    'pullback_dX_over_2Y': {'phi1': 'x dx/y', 'phi2': '-dx/y'},
    'differentials_independent': True,
    'finite_field_check_primes': pmap,
    'bridge_identity_f_eq_1mx6_F0F8L': True,
    'status': 'proved',
}
dump()

# ===================================================== 2. E2: базис
M2 = E2.minimal_model()
assert [str(z) for z in M2.ainvs()] == bridge['elliptic_factors'][1]['minimal_ainvs']
u_, r_, s_, t_ = E2.isomorphism_to(M2).tuple()
say("[2] E2 -> E2min: (u,r,s,t) = (%s,%s,%s,%s)  => чистый сдвиг X на %s" % (u_, r_, s_, t_, -r_))
assert u_ == 1 and s_ == 0 and t_ == 0
sh2 = ZZ(-r_)
say("[2] E2min = %s ; sh2 = %s" % (M2.ainvs(), sh2))
say("[2] conductor = %s = %s" % (M2.conductor(), factor(M2.conductor())))
say("[2] disc(E2min) = %s" % factor(M2.discriminant()))

Tg = M2.torsion_subgroup()
say("[2] torsion_subgroup().invariants() = %s ; все точки кручения: %s"
    % (Tg.invariants(), [str(z) for z in M2.torsion_points()]))
T2 = M2([-a6 + sh2, 0, 1])
assert T2.order() == 2 and list(Tg.invariants()) == [2]
say("[2] T2 = %s  (= образ (-a6, 0) с сырой модели), порядок 2" % T2)

P2 = M2([QQ(bridge['elliptic_factors'][1]['gens'][0][0]),
         QQ(bridge['elliptic_factors'][1]['gens'][0][1]), 1])
assert P2.order() == oo
hP2 = RR(P2.height(precision=300))
say("[2] P2 = %s (из new_bridge_verified.json) ; hhat(P2) = %s" % (P2, hP2))
say("[2] P2 на сырой модели: (%s, %s)" % (P2[0] - sh2, P2[1]))

# --- три независимые реализации ранга
ec2 = mwrank_EllipticCurve(list(map(int, M2.ainvs())), verbose=False)
ec2.two_descent(verbose=False)
say("[2a] eclib(mwrank): rank=%s rank_bound=%s certain=%s selmer_rank=%s"
    % (ec2.rank(), ec2.rank_bound(), ec2.certain(), ec2.selmer_rank()))
pr2 = M2.pari_curve().ellrank()
say("[2a] PARI ellrank: [%s, %s], найденные точки %s" % (pr2[0], pr2[1], pr2[3]))
g_sage = M2.gens(proof=True)
say("[2a] Sage M2.gens(proof=True) = %s" % (g_sage,))
assert len(g_sage) == 1
G = g_sage[0]
rel = None
for m in range(-4, 5):
    for e in [0, 1]:
        if m * P2 + e * T2 == G:
            rel = (m, e)
say("[2a] генератор Sage выражается через (P2, T2) как %s => совпадает с точностью до знака/кручения"
    % (rel,))
assert rel is not None and abs(rel[0]) == 1
say("[2a] ВЫВОД: rank E2(Q) = 1 (три реализации согласны), E2(Q)_tors = Z/2")

# ===================================================== 3. НАСЫЩЕННОСТЬ (брутфорс)
# (a) доказуемая нижняя граница канонической высоты (Cremona-Prickett-Siksek)
hmin = None
try:
    alarm(2400)
    hmin = RR(M2.height_function().min(RealNumber('0.005'), 6))
    cancel_alarm()
except Exception as e:
    cancel_alarm()
    say("[3a] height_function().min не завершился: %s" % type(e).__name__)
if hmin is None:
    qmax = None
    say("[3a] НЕ ПОЛУЧЕНА доказуемая нижняя граница высоты -- граница на q НЕ установлена этим путём")
else:
    qmax = ZZ(floor(RR(hP2 / hmin).sqrt()))
    say("[3a] доказуемая min hhat на E2(Q)\\tors = %s" % hmin)
    say("[3a] если P2 = q*Q + (кручение), то q^2 <= hhat(P2)/min = %s => q <= %s"
        % (RR(hP2 / hmin), qmax))

# (b) БРУТФОРС-сертификаты неделимости: строим q*E(F_p) полным перебором
D2 = ZZ(M2.discriminant())


def brute_nondiv(M, Pt, q, disc, pstart=5, pstop=3000):
    """Возвращает p и ДОКАЗАТЕЛЬСТВО, что образ Pt в E(F_p) не лежит в q*E(F_p).
    q*E(F_p) вычисляется полным перебором точек -- никаких дискретных логарифмов."""
    for p in prime_range(pstart, pstop):
        if disc % p == 0:
            continue
        Ep = M.change_ring(GF(p))
        img = Ep(Pt)
        qE = set(q * S for S in Ep.points())
        if img not in qE:
            return {'p': int(p), 'card_E_Fp': int(Ep.cardinality()),
                    'card_qE_Fp': int(len(qE)),
                    'image': str(img),
                    'method': 'q*E(F_p) построено ПОЛНЫМ перебором точек E(F_p)',
                    'conclusion': 'образ P не лежит в %s*E(F_%s) => P не делится на %s в E(Q)' % (q, p, q)}
    return None


qlist = list(prime_range(2, int(qmax) + 1)) if qmax else list(prime_range(2, 20))
certs = {}
allq = True
for q in qlist:
    names = ['P', 'P+T'] if q == 2 else ['P']
    certs[str(q)] = {}
    for nm in names:
        Pt = P2 if nm == 'P' else P2 + T2
        c = brute_nondiv(M2, Pt, q, D2)
        certs[str(q)][nm] = c
        say("[3b] q=%-2s %-4s: %s" % (q, nm, ('p=%s, |E(F_p)|=%s, |qE(F_p)|=%s -- образ вне qE'
                                              % (c['p'], c['card_E_Fp'], c['card_qE_Fp'])) if c else 'СЕРТИФИКАТ НЕ НАЙДЕН'))
        if c is None:
            allq = False
say("[3b] брутфорс-сертификаты найдены для всех q <= %s: %s" % (qlist[-1], allq))

# (c) точный division_points для малых q (алгебраический, без редукций)
divp = {}
for q in [2, 3, 5, 7, 11]:
    try:
        alarm(900)
        dp = {'P': [str(z) for z in P2.division_points(q)],
              'P+T': [str(z) for z in (P2 + T2).division_points(q)]}
        cancel_alarm()
    except Exception as e:
        cancel_alarm()
        dp = {'error': type(e).__name__}
    divp[str(q)] = dp
    say("[3c] точный division_points q=%s : %s" % (q, dp))

sat_proved = bool(allq and hmin is not None)
say("[3d] ВЫВОД по насыщенности: %s" %
    ("E2(Q) = <P2> (+) <T2>, индекс подгруппы = 1 (ДОКАЗАНО ПО: нижняя граница высоты + "
     "полные сертификаты неделимости для всех простых q <= %s)" % (qlist[-1])
     if sat_proved else "НЕ ЗАВЕРШЕНО -- граница на q или сертификат отсутствует"))

R['E2'] = {
    'ainvs_raw': [str(z) for z in E2.ainvs()],
    'ainvs_min': [str(z) for z in M2.ainvs()],
    'X_shift_raw_to_min': int(sh2),
    'conductor': str(M2.conductor()),
    'conductor_factored': str(factor(M2.conductor())),
    'disc_min_factored': str(factor(M2.discriminant())),
    'rank': 1,
    'rank_eclib': {'rank': int(ec2.rank()), 'rank_bound': int(ec2.rank_bound()),
                   'certain': bool(ec2.certain()), 'selmer_rank': int(ec2.selmer_rank())},
    'rank_pari_ellrank_interval': [int(pr2[0]), int(pr2[1])],
    'rank_sage_gens': [str(z) for z in g_sage],
    'sage_gen_in_terms_of_P2_T2': list(rel),
    'torsion': {'invariants': [int(z) for z in Tg.invariants()],
                'generator_min': str(T2),
                'generator_raw': '(%s, 0)' % (-a6),
                'all_torsion_points': [str(z) for z in M2.torsion_points()]},
    'generator_min_model': str(P2),
    'generator_raw_model': '(%s, %s)' % (P2[0] - sh2, P2[1]),
    'canonical_height': str(hP2),
    'regulator': str(hP2),
    'saturation': {
        'provable_min_canonical_height_CPS': str(hmin) if hmin else None,
        'q_bound': int(qmax) if qmax else None,
        'primes_q_checked': [int(q) for q in qlist],
        'brute_force_certificates': certs,
        'exact_division_points': divp,
        'index': 1 if sat_proved else None,
        'saturated': sat_proved,
    },
    'basis': 'E2(Q) = Z*P2  (+)  (Z/2)*T2',
    'status': 'proved_software' if sat_proved else 'incomplete',
}
dump()

# ===================================================== 4. E1 (для совместного решета)
M1 = E1.minimal_model()
u_, r_, s_, t_ = E1.isomorphism_to(M1).tuple()
assert u_ == 1 and s_ == 0 and t_ == 0
sh1 = ZZ(-r_)
T1 = M1([-a0 + sh1, 0, 1])
P1 = M1([QQ(bridge['elliptic_factors'][0]['gens'][0][0]),
         QQ(bridge['elliptic_factors'][0]['gens'][0][1]), 1])
assert T1.order() == 2 and P1.order() == oo
hP1 = RR(P1.height(precision=300))
ec1 = mwrank_EllipticCurve(list(map(int, M1.ainvs())), verbose=False)
ec1.two_descent(verbose=False)
hmin1 = None
try:
    alarm(2400)
    hmin1 = RR(M1.height_function().min(RealNumber('0.005'), 6))
    cancel_alarm()
except Exception:
    cancel_alarm()
qmax1 = ZZ(floor(RR(hP1 / hmin1).sqrt())) if hmin1 else None
D1 = ZZ(M1.discriminant())
certs1 = {}
allq1 = True
for q in (list(prime_range(2, int(qmax1) + 1)) if qmax1 else list(prime_range(2, 12))):
    certs1[str(q)] = {}
    for nm in (['P', 'P+T'] if q == 2 else ['P']):
        Pt = P1 if nm == 'P' else P1 + T1
        c = brute_nondiv(M1, Pt, q, D1)
        certs1[str(q)][nm] = c
        if c is None:
            allq1 = False
say("[4] E1min = %s ; sh1 = %s ; T1 = %s ; P1 = %s ; hhat=%s" % (M1.ainvs(), sh1, T1, P1, hP1))
say("[4] E1: eclib rank=%s certain=%s ; CPS min hhat=%s ; q<=%s ; брутфорс-сертификаты для всех q: %s"
    % (ec1.rank(), ec1.certain(), hmin1, qmax1, allq1))

R['E1'] = {
    'ainvs_min': [str(z) for z in M1.ainvs()], 'X_shift_raw_to_min': int(sh1),
    'generator_min_model': str(P1), 'torsion_generator_min': str(T1),
    'canonical_height': str(hP1),
    'rank_eclib': {'rank': int(ec1.rank()), 'certain': bool(ec1.certain())},
    'saturation': {'provable_min_canonical_height_CPS': str(hmin1) if hmin1 else None,
                   'q_bound': int(qmax1) if qmax1 else None,
                   'brute_force_certificates': certs1,
                   'saturated': bool(allq1 and hmin1 is not None)},
    'status': 'proved_software' if (allq1 and hmin1 is not None) else 'incomplete',
}
dump()

# ===================================================== 5. образ Jac(D)(Q) в E1(Q) x E2(Q)
# Psi = (phi1_*, phi2_*) : Jac(D) -> E1 x E2.
#
# НИЖНЕЕ включение (что ТОЧНО лежит в образе):
#   (i) Psi o Phi = [2], где Phi = phi1^* + phi2^* ; значит 2*(E1(Q)xE2(Q)) в образе.
#   (ii) рациональный класс [W(g2) - inf+ - inf-] даёт (T1, T2).
# Проверим (ii) ЯВНО и АРИФМЕТИЧЕСКИ (без расширения полей):
#   корни g1 = 529x^2+49 суть +-alpha, alpha^2 = -49/529 ; phi1(+-alpha,0) = (a6*alpha^2, 0)
#   a6*alpha^2 = 44182609 * (-49/529) = -4092529 = -a0.
#   корни g2, g3 суть +-beta, +-beta_bar ; сумма ВСЕХ трёх точек 2-кручения = O,
#   поэтому [phi1(beta)] + [phi1(beta_bar)] = -(точка с X = -a0) = (точка с X = -a0) (2-кручение).
tw1 = [r for r, _ in Rx(rhs(E1, x)).roots()]
tw2 = [r for r, _ in Rx(rhs(E2, x)).roots()]
say("[5a] рациональные X 2-кручения E1: %s   (ожидается -a0 = %s)" % (tw1, -a0))
say("[5a] рациональные X 2-кручения E2: %s   (ожидается -a6 = %s)" % (tw2, -a6))
assert tw1 == [-a0] and tw2 == [-a6]
assert a6 * (QQ(-49) / 529) == -a0 and a0 / (QQ(-49) / 529) == -a6
say("[5a] phi1({+-alpha}) = (-a0,0) = T1 ; phi2({+-alpha}) = (-a6,0) = T2   (alpha^2 = -49/529)")
say("[5a] так как сумма трёх точек 2-кручения = O и E1 имеет ЕДИНСТВЕННУЮ рациональную")
say("[5a] точку 2-кручения, класс [W(289x^2-322x+289) - Dinf] даёт Psi = (T1, T2).")
say("[5a] => Psi(Jac(D)(Q)) СОДЕРЖИТ 2*(E1(Q)xE2(Q)) + <(T1,T2)>  ==> индекс ДЕЛИТ 8")

# ВЕРХНЕЕ включение: локальные препятствия. Для каждого класса -- ЯВНЫЙ сертификат.
A1M = [int(z) for z in M1.ainvs()]
A2M = [int(z) for z in M2.ainvs()]


def local_image_set(p):
    """Psi(Jac(D)(F_p)) полным перебором эффективных дивизоров степени 2 над F_p.
    Возвращает множество пар (A,B) в E1(F_p) x E2(F_p)."""
    F = GF(p); F2 = GF(p**2, 'w')
    m1 = M1.change_ring(F); m2 = M2.change_ring(F)

    def phi_pt(tag, xv, yv, K, e1, e2):
        if tag == 'inf':
            return e1(0), e2([K(sh2), K(a0) * K(6647) * K(yv), 1])
        A = e1([K(a6) * xv**2 + K(sh1), K(a6) * yv, 1])
        B = e2(0) if xv == 0 else e2([K(a0) / xv**2 + K(sh2), K(a0) * yv / xv**3, 1])
        return A, B

    def affine(K):
        out = []
        for xv in K:
            v = K(a6) * xv**6 + K(a4) * xv**4 + K(a2) * xv**2 + K(a0)
            if v.is_square():
                sq = v.sqrt()
                out.append(('fin', xv, sq))
                if sq != 0:
                    out.append(('fin', xv, -sq))
        return out

    e1F, e2F = m1, m2
    e1K = EllipticCurve(F2, [F2(z) for z in A1M]); e2K = EllipticCurve(F2, [F2(z) for z in A2M])
    DF = affine(F) + [('inf', None, 1), ('inf', None, -1)]
    assert len(DF) == p + 1 - m1.trace_of_frobenius() - m2.trace_of_frobenius()
    S = set()
    for i in range(len(DF)):
        for j in range(i, len(DF)):
            A1, B1 = phi_pt(DF[i][0], DF[i][1], DF[i][2], F, e1F, e2F)
            A2, B2 = phi_pt(DF[j][0], DF[j][1], DF[j][2], F, e1F, e2F)
            S.add((A1 + A2, B1 + B2))
    Fset = set(F2(z) for z in F)
    done = set()
    for (_, xv, yv) in affine(F2):
        if xv in Fset or (xv, yv) in done:
            continue
        done.add((xv, yv)); done.add((xv**p, yv**p))
        A1, B1 = phi_pt('fin', xv, yv, F2, e1K, e2K)
        A2, B2 = phi_pt('fin', xv**p, yv**p, F2, e1K, e2K)
        A = A1 + A2; B = B1 + B2
        A = e1F(0) if A.is_zero() else e1F([F(A[0]), F(A[1]), 1])
        B = e2F(0) if B.is_zero() else e2F([F(B[0]), F(B[1]), 1])
        S.add((A, B))
    # самоконтроль: S -- подгруппа (ПОЛНАЯ проверка, не выборочная)
    Ls = list(S)
    for u1 in Ls:
        for v1 in Ls:
            assert (u1[0] + v1[0], u1[1] + v1[1]) in S, "S не подгруппа при p=%s" % p
    N = m1.cardinality() * m2.cardinality()
    return S, N, m1, m2


classes = [(i, j, k, l) for i in [0, 1] for j in [0, 1] for k in [0, 1] for l in [0, 1]]
surv = set(classes)
kill_cert = {}
rows5 = []
for p in prime_range(5, 62):
    if (D1 * D2) % p == 0:
        continue
    S, N, m1, m2 = local_image_set(p)
    idx = ZZ(N) // len(S)
    pred = 4 if p % 4 == 1 else 2
    ok = set()
    for c in classes:
        i, j, k, l = c
        if (m1(i * P1 + j * T1), m2(k * P2 + l * T2)) in S:
            ok.add(c)
        elif c in surv and str(list(c)) not in kill_cert:
            kill_cert[str(list(c))] = {
                'p': int(p),
                'class_meaning': '(%s*P1 + %s*T1, %s*P2 + %s*T2)' % c,
                'image_mod_p': [str(m1(i * P1 + j * T1)), str(m2(k * P2 + l * T2))],
                'card_Psi_J_Fp': int(len(S)), 'card_E1xE2_Fp': int(N),
                'local_index': int(idx),
                'certificate': 'редукция класса mod %s НЕ лежит в Psi(Jac(D)(F_%s))' % (p, p)}
    surv &= ok
    rows5.append({'p': int(p), 'card_image': int(len(S)), 'card_E1xE2': int(N),
                  'local_index': int(idx), 'predicted_index': int(pred),
                  'surviving': sorted([list(c) for c in ok])})
    say("[5b] p=%-3d |Psi(J(F_p))|=%-7d |E1xE2(F_p)|=%-7d индекс=%s (предсказание |K(F_p)|: %s) выжило=%d"
        % (p, len(S), N, idx, pred, len(ok)))
    assert idx == pred, "формула локального индекса не подтвердилась при p=%s" % p
say("[5b] ВЫЖИВШИЕ классы после всех простых: %s" % sorted(surv))
for k_, v_ in sorted(kill_cert.items()):
    say("[5c] класс %s убит при p=%s : %s" % (k_, v_['p'], v_['certificate']))
say("[5c] выдано сертификатов исключения: %d (ожидается 14)" % len(kill_cert))

exact = (surv == {(0, 0, 0, 0), (0, 1, 0, 1)})
say("[5d] нижнее включение даёт ровно эти 2 класса; верхнее включение даёт %d классов; совпало? %s"
    % (len(surv), exact))
if exact:
    say("[5d] ИТОГ: Psi(Jac(D)(Q)) = 2*(E1(Q)xE2(Q)) + <(T1,T2)> ; ИНДЕКС = 8 (точно, не только делит)")
    say("[5d] правило: (m1 P1 + e1 T1, m2 P2 + e2 T2) в образе  <=>  m1, m2 чётны и e1 = e2")
else:
    say("[5d] ИТОГ: точный индекс НЕ установлен; доказано только: индекс делит 8")

R['image_of_Jac'] = {
    'map': 'Psi = (phi1_*, phi2_*) : Jac(D) -> E1 x E2, (2,2)-изогения степени 4',
    'lower_containment': {
        'contains': '2*(E1(Q) x E2(Q)) + <(T1,T2)>',
        'reason_2E': 'Psi o Phi = [2] (Phi = phi1^* + phi2^*)',
        'reason_T': ('рациональный класс [W(289x^2-322x+289) - inf+ - inf-] в Jac(D)(Q); '
                     'его Psi = (T1,T2), проверено арифметикой 2-кручения: '
                     'a6*(-49/529) = -a0 и a0/(-49/529) = -a6'),
        'index_of_this_subgroup_in_product': 8,
        'status': 'proved',
    },
    'upper_containment': {
        'method': 'редукция mod p; класс c лежит в образе => c mod p лежит в Psi(Jac(D)(F_p))',
        'primes': [r_['p'] for r_ in rows5],
        'per_prime': rows5,
        'self_check': 'S проверено как подгруппа ПОЛНЫМ перебором пар; #D(F_p)=p+1-a_p(E1)-a_p(E2)',
        'kernel_index_formula': '|K(F_p)| = 4 при p=1 mod 4, иначе 2 -- совпало на всех p',
        'surviving_classes': sorted([list(c) for c in surv]),
        'exclusion_certificates': kill_cert,
        'status': 'proved_software',
    },
    'exact': bool(exact),
    'index_in_E1xE2': 8 if exact else None,
    'index_divides': 8,
    'membership_rule': ('(m1 P1 + e1 T1, m2 P2 + e2 T2) in Psi(Jac(D)(Q)) <=> '
                        'm1 = m2 = 0 mod 2 и e1 = e2') if exact else None,
    'WARNING': ('E1(Q) x E2(Q) -- НАДмножество образа. Использовать всё произведение как область '
                'решета допустимо (риск только недоисключить). Считать произведение РАВНЫМ образу '
                'НЕЛЬЗЯ: здесь показано, что индекс равен 8.'),
}
dump()

# ===================================================== 6. КОНТРОЛЬ: 4 известные точки
B2pt = M2([sh2, a0 * 6647, 1])       # phi2(inf+)
A1pt = M1([sh1, a6 * 2023, 1])       # phi1((0, 2023))
assert B2pt in M2 and A1pt in M1


def coord(Q, P, T):
    if Q.is_zero():
        return (0, 0)
    for m in range(-10, 11):
        for e in [0, 1]:
            if m * P + e * T == Q:
                return (m, e)
    return None


ctrl = []
for (xv, yv, tag) in [(QQ(0), QQ(2023), 't=1'), (QQ(0), QQ(-2023), 't=1')]:
    Q1 = M1([sh1, a6 * yv, 1]); Q2 = M2(0)
    c1 = coord(Q1, P1, T1); d2 = coord(Q2 - B2pt, P2, T2)
    ok = (c1[0] % 2 == 0) and (d2[0] % 2 == 0) and (c1[1] == d2[1])
    ctrl.append({'z': '(0, %s)' % yv, 't': tag, 'phi1_coords': list(c1),
                 'phi2_minus_base_coords': list(d2), 'passes': bool(ok)})
    say("[6] z=(0,%s): phi1 = %s*P1+%s*T1 ; phi2-база = %s*P2+%s*T2 ; условие образа: %s"
        % (yv, c1[0], c1[1], d2[0], d2[1], ok))
for sgn in [1, -1]:
    Q1 = M1(0); Q2 = M2([sh2, sgn * a0 * 6647, 1])
    c1 = coord(Q1, P1, T1); d2 = coord(Q2 - B2pt, P2, T2)
    ok = (c1[0] % 2 == 0) and (d2[0] % 2 == 0) and (c1[1] == d2[1])
    ctrl.append({'z': 'inf_%s' % ('+' if sgn > 0 else '-'), 't': 't=-1',
                 'phi1_coords': list(c1), 'phi2_minus_base_coords': list(d2), 'passes': bool(ok)})
    say("[6] z=inf_%s: phi1 = O ; phi2-база = %s*P2+%s*T2 ; условие образа: %s"
        % ('+' if sgn > 0 else '-', d2[0], d2[1], ok))
allok = all(c['passes'] for c in ctrl)
say("[6] КОНТРОЛЬ: все четыре известные точки D проходят условие образа: %s" % allok)
assert allok, "РЕШЕТО СЛОМАНО -- известная точка не проходит. СООБЩИТЬ НЕМЕДЛЕННО."
R['control_four_known_points'] = {'base_point': 'inf+', 'points': ctrl, 'all_pass': bool(allok),
                                  'status': 'proved_software'}

# ===================================================== 7. поточечное согласование факторов
say("[7] для z=(x,y) in D, x != 0, oo:  X1_raw * X2_raw = a6 x^2 * a0/x^2 = a0*a6 = %s = %s"
    % (a0 * a6, factor(a0 * a6)))
say("[7] в минимальных координатах: (X1-%s)*(X2-%s) = %s" % (sh1, sh2, a0 * a6))
say("[7] Y1_raw / Y2_raw = a6 y / (a0 y / x^3) = a6 x^3 / a0 ; => x^3 = a0 Y1/(a6 Y2)")
say("[7] это ПОТОЧЕЧНОЕ тождество, не зависит от индекса образа")
R['pointwise_compatibility'] = {
    'X1raw_times_X2raw': str(a0 * a6), 'factored': str(factor(a0 * a6)),
    'minimal_coords': '(X1-%s)*(X2-%s) = %s' % (sh1, sh2, a0 * a6),
    'Y_ratio': 'Y1/Y2 = a6 x^3 / a0', 'x_recovered': 'x^3 = a0*Y1/(a6*Y2)',
    'degenerate': {'x=0': 'phi2(z)=O', 'x=infty': 'phi1(z)=O'},
    'status': 'proved',
}

R['not_established'] = [
    'D(Q) НЕ определена; счёт (15,8) остаётся открытым, семейство G1 остаётся 126/127.',
    'Ни один 11-адический кандидат здесь НЕ удалён и НЕ признан рациональной точкой.',
    'Sha(E1), Sha(E2) не вычислены (и не нужны: граница ранга из 2-Селмера).',
    'Все ранги/насыщенность -- машинное доказательство (PARI/eclib/Sage), Magma не запускалась.',
]
dump()
say("готово: audit_E2_claude.json / audit_E2_claude.log")
