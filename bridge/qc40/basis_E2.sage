#!/usr/bin/env sage
# basis_E2.sage -- полный базис и НАСЫЩЕННОСТЬ фактора E2 для моста (15,8);
#                  символьная проверка карт phi1, phi2 : D -> E1, E2;
#                  ТОЧНЫЙ образ Jac(D)(Q) в E1(Q) x E2(Q) и граница индекса.
#
# Запуск:  sage /home/kep/magicKube/bridge/qc40/basis_E2.sage
# Выход:   basis_E2.json, basis_E2.log  (в каталоге скрипта)
#
# Метки статуса внутри JSON: proved / proved_software / numerical / observation.
# НИЧЕГО в этом файле не утверждает, что D(Q) состоит ровно из четырёх точек.

from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve
from pathlib import Path
import json, io, contextlib, random, os, sys, tempfile

OUT = Path('/home/kep/magicKube/bridge/qc40')   # sage run_file не выставляет __file__
R = {}
LOG = []
random.seed(int(20260912))


def say(*a):
    s = " ".join(str(z) for z in a)
    LOG.append(s)
    print(s, flush=True)


def jdefault(o):
    return str(o)


def capture_fd(func):
    """eclib печатает на уровне C, поэтому redirect_stdout не ловит его.
    Перехватываем сам дескриптор 1."""
    sys.stdout.flush()
    saved = os.dup(1)
    tf = tempfile.TemporaryFile(mode='w+b')
    try:
        os.dup2(tf.fileno(), 1)
        res = func()
    finally:
        sys.stdout.flush()
        os.dup2(saved, 1)
        os.close(saved)
        tf.seek(0)
        txt = tf.read().decode('utf-8', 'replace')
        tf.close()
    return res, txt


# ================================================================ 0. модель D
a0, a2, a4, a6 = QQ(4092529), QQ(47287151), QQ(37608911), QQ(44182609)
Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()
Kx = Rx.fraction_field()
f = a6 * x**6 + a4 * x**4 + a2 * x**2 + a0

t = (1 + x) / (1 - x)
s = QQ(289) / 2
h = (225 + 64 * t * t) * (64 + 225 * t * t) * (s * (1 + t * t) - 240 * t)
assert Kx(f) == (1 - x)**6 * h
assert gcd(f, f.derivative()) == 1
assert a0 == 2023**2 and a6 == 6647**2
say("[0] D: y^2 = %s ; gcd(f,f')=1 (род 2), тождество моста F0*F8*L подтверждено" % f)

# f распадается на три РАЦИОНАЛЬНЫХ квадратичных множителя -> Jac(D)[2](Q) = (Z/2)^2
g1 = 529 * x**2 + 49
g2 = 289 * x**2 - 322 * x + 289
g3 = 289 * x**2 + 322 * x + 289
assert g1 * g2 * g3 == 529 * 289 * 289 * f / a6 * 1
assert f == a6 / QQ(529 * 289 * 289) * g1 * g2 * g3
assert all(gg.is_irreducible() for gg in [g1, g2, g3])
disc_f = Rx(f).discriminant()
bad_f = sorted([int(p) for p, _ in factor(ZZ(disc_f * a6))])
say("[0] f = (a6/(529*289^2)) * (529x^2+49)(289x^2-322x+289)(289x^2+322x+289), все три неприводимы над Q")
say("[0] => Jac(D)[2](Q) = (Z/2)^2 (три рациональных квадратичных множителя)")
say("[0] простые деления disc(f)*a6: %s" % bad_f)

R['curve_D'] = {
    'f': str(f),
    'a0a2a4a6': [str(a0), str(a2), str(a4), str(a6)],
    'genus': 2,
    'bridge_identity': 'proved',
    'quadratic_factors': [str(g1), str(g2), str(g3)],
    'jac_two_torsion_rational': '(Z/2)^2',
    'bad_primes_of_model': bad_f,
    'status': 'proved',
}

# ================================ 1. символьная проверка карт phi1, phi2
# E1: Y^2 = X^3 + a4 X^2 + (a2 a6) X + a0 a6^2 ,  phi1(x,y) = (a6 x^2, a6 y)
# E2: Y^2 = X^3 + a2 X^2 + (a4 a0) X + a6 a0^2 ,  phi2(x,y) = (a0/x^2, a0 y/x^3)
E1 = EllipticCurve([0, a4, 0, a2 * a6, a0 * a6**2])
E2 = EllipticCurve([0, a2, 0, a4 * a0, a6 * a0**2])


def rhs(E, T):
    A = E.ainvs()
    return T**3 + A[1] * T**2 + A[3] * T + A[4]


P2r = PolynomialRing(QQ, ['xx', 'yy'])
xx, yy = P2r.gens()
fxx = a6 * xx**6 + a4 * xx**4 + a2 * xx**2 + a0
I = P2r.ideal(yy**2 - fxx)

r1sym = (a6 * yy)**2 - rhs(E1, a6 * xx**2)
ok1 = (r1sym in I)
q1 = P2r(r1sym).quo_rem(yy**2 - fxx)[0] if ok1 else None
say("[1a] phi1: (a6 y)^2 - RHS_E1(a6 x^2) = (%s)*(y^2 - f)  ->  %s ; множитель = a6^2 ? %s"
    % (q1, ok1, q1 == a6**2))

lhs2 = (a0 * yy)**2
rhs2 = (a0**3 + a2 * a0**2 * xx**2 + a4 * a0 * a0 * xx**4 + a6 * a0**2 * xx**6)
r2sym = lhs2 - rhs2
ok2 = (r2sym in I)
q2 = P2r(r2sym).quo_rem(yy**2 - fxx)[0] if ok2 else None
say("[1b] phi2: x^6*[(a0 y/x^3)^2 - RHS_E2(a0/x^2)] = (%s)*(y^2 - f)  ->  %s ; множитель = a0^2 ? %s"
    % (q2, ok2, q2 == a0**2))
assert ok1 and ok2

assert a6**2 * f == rhs(E1, a6 * x**2)
assert a0**2 * Kx(f) / x**6 == rhs(E2, Kx(a0 / x**2))
say("[1c] те же тождества в поле дробей Q(x): подтверждены (независимая форма)")

F2f = P2r.fraction_field()
xg, yg = F2f(xx), F2f(yy)
phi1sym = (a6 * xg**2, a6 * yg)
phi2sym = (a0 / xg**2, a0 * yg / xg**3)
sub_sigma = lambda z: z.subs({xx: -xx, yy: yy})
sub_tau = lambda z: z.subs({xx: -xx, yy: -yy})
inv1 = all(sub_sigma(c) == c for c in phi1sym)
inv2 = all(sub_tau(c) == c for c in phi2sym)
say("[1d] phi1 инвариантна под sigma:(x,y)->(-x,y):  %s   => E1 = D/sigma, deg phi1 = 2" % inv1)
say("[1d] phi2 инвариантна под tau:(x,y)->(-x,-y):   %s   => E2 = D/tau,   deg phi2 = 2" % inv2)
assert inv1 and inv2
say("[1e] phi1^*(dX/2Y) = x dx/y ,  phi2^*(dX/2Y) = -dx/y  -- базис H^0(D,Omega); независимы")

num_checks = []
for xv in [QQ(2), QQ(3), QQ(1) / 2, QQ(-5) / 3, QQ(7), QQ(-11) / 7]:
    yv2 = f(xv)
    if yv2.is_square():
        yv = yv2.sqrt()
        ok = ((a6 * yv)**2 == rhs(E1, a6 * xv**2)) and \
             ((a0 * yv / xv**3)**2 == rhs(E2, a0 / xv**2))
        num_checks.append([str(xv), 'rational y', bool(ok)])
    else:
        ok1n = (a6**2 * yv2 == rhs(E1, a6 * xv**2))
        ok2n = (a0**2 * yv2 / xv**6 == rhs(E2, a0 / xv**2))
        num_checks.append([str(xv), 'y^2 only', bool(ok1n and ok2n)])
say("[1f] числовые подстановки: %s" % num_checks)
assert all(c[2] for c in num_checks)

R['maps'] = {
    'phi1': '(a6*x^2, a6*y)', 'phi2': '(a0/x^2, a0*y/x^3)',
    'E1_ainvs_raw': [str(z) for z in E1.ainvs()],
    'E2_ainvs_raw': [str(z) for z in E2.ainvs()],
    'symbolic_identity_phi1': 'proved, cofactor a6^2',
    'symbolic_identity_phi2': 'proved, cofactor a0^2',
    'phi1_quotient_by': 'sigma(x,y)=(-x,y)',
    'phi2_quotient_by': 'tau(x,y)=(-x,-y)',
    'degrees': [2, 2],
    'pullback_differentials': ['x dx / y', '- dx / y'],
    'numeric_substitutions': num_checks,
    'status': 'proved',
}

# ================================================ 2. E2: полный базис и ранг
M2 = E2.minimal_model()
iso2 = E2.isomorphism_to(M2)
r2shift = ZZ(15762384)
assert M2([r2shift, a0 * 6647, 1]) == M2(iso2(E2([0, a0 * 6647, 1])))
say("[2] E2 raw  = %s" % (E2.ainvs(),))
say("[2] E2 min  = %s   (X_min = X_raw + %s)" % (M2.ainvs(), r2shift))
say("[2] disc(E2min) = %s" % factor(M2.discriminant()))
say("[2] conductor   = %s = %s" % (M2.conductor(), factor(M2.conductor())))

T2grp = M2.torsion_subgroup()
T2 = M2([-28420225, 0, 1])
assert T2.order() == 2 and list(T2grp.invariants()) == [2]
say("[2] torsion E2(Q)_tors = %s, генератор T2 = %s" % (T2grp.invariants(), T2))

P2 = M2([-2479007, 87162492480, 1])
assert P2.order() == oo
hP2 = RR(P2.height(precision=200))
say("[2] P2 = %s ; hhat(P2) = %s" % (P2, hP2))

ec2 = mwrank_EllipticCurve(list(map(int, M2.ainvs())), verbose=False)
ec2.two_descent(verbose=False)
say("[2] eclib: rank=%s rank_bound=%s certain=%s selmer_rank=%s gens=%s"
    % (ec2.rank(), ec2.rank_bound(), ec2.certain(), ec2.selmer_rank(), ec2.gens()))
assert ec2.rank() == 1 and ec2.certain()

pr2 = M2.pari_curve().ellrank()
say("[2] PARI ellrank: интервал [%s,%s], точки %s" % (pr2[0], pr2[1], pr2[3]))
say("[2] root_number(E2) = %s (согласуется с нечётным рангом; НЕ используется как довод)"
    % M2.root_number())

# насыщение eclib с ЗАХВАТОМ его собственной границы индекса (не хардкод!)
sat2, vtxt2 = capture_fd(lambda: M2.saturation([P2], verbose=True))
vlines2 = [z for z in vtxt2.strip().splitlines() if z.strip()]
say("[2] eclib saturation: pts=%s index=%s reg=%s" % (sat2[0], sat2[1], sat2[2]))
for z in vlines2:
    say("[2]   eclib| %s" % z)
sat2_bound = None
for z in vlines2:
    if 'Saturation index bound' in z:
        sat2_bound = int(z.split('=')[-1].strip())
sat2_primes = None
for z in vlines2:
    if 'Saturation primes are now' in z:
        sat2_primes = [int(w) for w in z.split('[')[-1].split(']')[0].split()]
say("[2] ЗАХВАЧЕНО из eclib: граница индекса = %s, проверенные простые = %s"
    % (sat2_bound, sat2_primes))
assert sat2[1] == 1

gens2 = M2.gens(proof=True)
say("[2] M2.gens(proof=True) = %s" % gens2)
assert len(gens2) == 1

# ================ 3. НЕЗАВИСИМАЯ проверка насыщенности: нижняя граница высоты + mod p
# (a) доказуемая нижняя граница канонической высоты неторсионных точек
hmin2 = None
try:
    alarm(1800)
    hmin2 = RR(M2.height_function().min(RealNumber('0.01'), 5))
    cancel_alarm()
    say("[3a] доказуемая нижняя граница min hhat на E2(Q)\\tors = %s" % hmin2)
except Exception as e:
    cancel_alarm()
    say("[3a] height_function().min не завершился: %s" % type(e).__name__)
qmax2 = int(RR(hP2 / hmin2).sqrt()) if hmin2 else 20
say("[3a] => если P2 = q*Q, то q^2 <= hhat(P2)/min = %s, т.е. q <= %s"
    % (RR(hP2 / hmin2) if hmin2 else 'н/д', qmax2))

# (b) сертификаты неделимости через редукции mod p


def nondiv_cert(M, Pt, q, disc):
    """Ищет p хорошей редукции, при котором образ Pt не лежит в q*E(F_p).
    Найденное p -- КОНЕЧНЫЙ сертификат того, что Pt не делится на q в E(Q)."""
    for p in prime_range(5, 8000):
        if disc % p == 0:
            continue
        Ep = M.change_ring(GF(p))
        A = Ep.abelian_group()
        inv = list(A.invariants())
        v = [ZZ(z) for z in A.discrete_log(Ep(Pt))]
        if any(vi % gcd(q, ZZ(ni)) != 0 for ni, vi in zip(inv, v)):
            return {'p': int(p), 'group': [int(z) for z in inv],
                    'discrete_log': [int(z) for z in v],
                    'reason': 'image not in %s*E(F_%s)' % (q, p)}
    return None


D2 = ZZ(M2.discriminant())
qs2 = [q for q in prime_range(2, max(qmax2 + 1, 20))]
sat_certs2 = {}
for q in qs2:
    names = ['P', 'P+T'] if q == 2 else ['P']
    sat_certs2[str(q)] = {}
    for nm in names:
        Pt = P2 if nm == 'P' else P2 + T2
        c = nondiv_cert(M2, Pt, q, D2)
        sat_certs2[str(q)][nm] = c
        say("[3b] E2, q=%s, %s: сертификат неделимости %s" % (q, nm, c))
        assert c is not None, "не найден сертификат неделимости E2 q=%s %s" % (q, nm)

divpts2 = {}
for q in [2, 3, 5, 7]:
    try:
        alarm(300)
        divpts2[str(q)] = {'P': [str(z) for z in P2.division_points(q)],
                           'P+T': [str(z) for z in (P2 + T2).division_points(q)]}
    except Exception as e:
        divpts2[str(q)] = {'error': type(e).__name__}
    finally:
        cancel_alarm()
    say("[3c] E2 division_points q=%s : %s" % (q, divpts2[str(q)]))

R['E2'] = {
    'ainvs_raw': [str(z) for z in E2.ainvs()],
    'ainvs_min': [str(z) for z in M2.ainvs()],
    'X_min_minus_X_raw': int(r2shift),
    'conductor': str(M2.conductor()),
    'conductor_factored': str(factor(M2.conductor())),
    'disc_min_factored': str(factor(M2.discriminant())),
    'torsion_invariants': [int(z) for z in T2grp.invariants()],
    'torsion_generator': str(T2),
    'generator_min_model': str(P2),
    'generator_raw_model': str([P2[0] - r2shift, P2[1]]),
    'canonical_height': str(hP2),
    'eclib': {'rank': int(ec2.rank()), 'rank_bound': int(ec2.rank_bound()),
              'certain': bool(ec2.certain()), 'selmer_rank': int(ec2.selmer_rank()),
              'gens': [str(g) for g in ec2.gens()]},
    'pari_ellrank_interval': [int(pr2[0]), int(pr2[1])],
    'eclib_saturation_index': int(sat2[1]),
    'eclib_saturation_index_bound_CAPTURED': sat2_bound,
    'eclib_saturation_primes_checked_CAPTURED': sat2_primes,
    'gens_proof_true': [str(g) for g in gens2],
    'provable_min_canonical_height': str(hmin2) if hmin2 else None,
    'q_bound_from_height': int(qmax2),
    'independent_nondivisibility_certificates': sat_certs2,
    'division_points': divpts2,
    'basis_statement': 'E2(Q) = <P2> (+) <T2>, P2 бесконечного порядка, T2 порядка 2',
    'saturated': True,
    'status': 'proved_software',
    'status_note': ('ранг 1 подтверждён двумя независимыми реализациями (eclib certain=True, PARI ellrank '
                    '[1,1]); насыщенность подтверждена (i) eclib saturation с его собственной строгой границей '
                    'индекса, (ii) независимо: доказуемой нижней границей высоты плюс сертификаты '
                    'неделимости mod p для всех q до этой границы'),
}

# ========================================== 4. E1: то же самое, для совместного решета
M1 = E1.minimal_model()
r1shift = ZZ(12536304)
assert M1([r1shift, a6 * 2023, 1]) == M1(E1.isomorphism_to(M1)(E1([0, a6 * 2023, 1])))
T1grp = M1.torsion_subgroup()
T1 = M1([8443775, 0, 1])
P1 = M1([34384993, 286391046720, 1])
hP1 = RR(P1.height(precision=200))
assert T1.order() == 2 and list(T1grp.invariants()) == [2]
ec1 = mwrank_EllipticCurve(list(map(int, M1.ainvs())), verbose=False)
ec1.two_descent(verbose=False)
sat1, vtxt1 = capture_fd(lambda: M1.saturation([P1], verbose=True))
vlines1 = [z for z in vtxt1.strip().splitlines() if z.strip()]
sat1_bound = None
sat1_primes = None
for z in vlines1:
    if 'Saturation index bound' in z:
        sat1_bound = int(z.split('=')[-1].strip())
    if 'Saturation primes are now' in z:
        sat1_primes = [int(w) for w in z.split('[')[-1].split(']')[0].split()]
hmin1 = None
try:
    alarm(1800)
    hmin1 = RR(M1.height_function().min(RealNumber('0.01'), 5))
    cancel_alarm()
except Exception:
    cancel_alarm()
qmax1 = int(RR(hP1 / hmin1).sqrt()) if hmin1 else 10
D1 = ZZ(M1.discriminant())
sat_certs1 = {}
for q in prime_range(2, max(qmax1 + 1, 14)):
    names = ['P', 'P+T'] if q == 2 else ['P']
    sat_certs1[str(q)] = {}
    for nm in names:
        Pt = P1 if nm == 'P' else P1 + T1
        c = nondiv_cert(M1, Pt, q, D1)
        sat_certs1[str(q)][nm] = c
        assert c is not None
say("[4] E1 min = %s ; tors %s ; T1 = %s ; P1 = %s ; hhat = %s" % (M1.ainvs(), T1grp.invariants(), T1, P1, hP1))
say("[4] E1 eclib rank=%s certain=%s ; saturation index=%s, граница eclib=%s (простые %s) ; min hhat=%s, q<=%s"
    % (ec1.rank(), ec1.certain(), sat1[1], sat1_bound, sat1_primes, hmin1, qmax1))
say("[4] E1 сертификаты неделимости: %s" % {k: {n: (v['p'] if v else None) for n, v in d.items()}
                                            for k, d in sat_certs1.items()})
assert ec1.rank() == 1 and ec1.certain() and sat1[1] == 1

R['E1'] = {
    'ainvs_raw': [str(z) for z in E1.ainvs()],
    'ainvs_min': [str(z) for z in M1.ainvs()],
    'X_min_minus_X_raw': int(r1shift),
    'torsion_invariants': [int(z) for z in T1grp.invariants()],
    'torsion_generator': str(T1),
    'generator_min_model': str(P1),
    'generator_raw_model': str([P1[0] - r1shift, P1[1]]),
    'canonical_height': str(hP1),
    'eclib': {'rank': int(ec1.rank()), 'rank_bound': int(ec1.rank_bound()),
              'certain': bool(ec1.certain())},
    'eclib_saturation_index': int(sat1[1]),
    'eclib_saturation_index_bound_CAPTURED': sat1_bound,
    'eclib_saturation_primes_checked_CAPTURED': sat1_primes,
    'provable_min_canonical_height': str(hmin1) if hmin1 else None,
    'q_bound_from_height': int(qmax1),
    'independent_nondivisibility_certificates': sat_certs1,
    'basis_statement': 'E1(Q) = <P1> (+) <T1>, T1 порядка 2',
    'saturated': True,
    'status': 'proved_software',
}

# ============================ 5. Psi o Phi = [2]: прямое доказательство на дивизорах
# Psi = (phi1_*, phi2_*) : Jac(D) -> E1 x E2 ;  Phi = phi1^* + phi2^* : E1 x E2 -> Jac(D).
# phi1^*([Q]-[O]) = [(x0,y0) + (-x0,y0) - inf+ - inf-],  x0^2 = X(Q)/a6, y0 = Y(Q)/a6.
#   phi1_* этого  = Q + Q = 2Q                        (deg phi1 = 2)
#   phi2_* этого  = (a0/x0^2, a0 y0/x0^3) + (a0/x0^2, -a0 y0/x0^3) - B2 - (-B2) = O
# phi2^*([Q']-[O]) = [(x0,y0) + (-x0,-y0) - P0+ - P0-],  a0/x0^2 = X(Q'), y0 = Y(Q') x0^3/a0.
#   phi2_* этого  = 2Q' ;  phi1_* этого = (a6x0^2, a6y0) + (a6x0^2, -a6y0) - A1 - (-A1) = O
# Оба крест-члена равны нулю ТОЖДЕСТВЕННО (сокращение противоположных точек), без ссылки
# на Hom(E2,E1)=0. Отсюда Psi o Phi = [2] на E1 x E2.
say("[5] Psi o Phi = [2] на E1 x E2: крест-члены phi2_* phi1^* и phi1_* phi2^* обнуляются")
say("[5]     тождественно (прообраз состоит из точки и её отражения, их образы противоположны).")

# второй, независимый довод: Hom_Q(E2,E1) = 0, так как кривые не изогенны над Q
noniso = None
for p in prime_range(5, 500):
    if D1 % p == 0 or D2 % p == 0:
        continue
    t1 = M1.change_ring(GF(p)).trace_of_frobenius()
    t2 = M2.change_ring(GF(p)).trace_of_frobenius()
    if t1 != t2:
        noniso = {'p': int(p), 'a_p(E1)': int(t1), 'a_p(E2)': int(t2)}
        break
say("[5] второй довод: E1 не изогенна E2 над Q, сертификат %s => Hom_Q(E2,E1)=0" % noniso)

# численный контроль Psi o Phi = [2] на конкретных точках (над расширением)
pp_checks = []
for Q in [P1, P1 + T1, 3 * P1, 3 * P1 + T1]:
    Xr = Q[0] - r1shift
    if Xr == 0:                      # x0 = 0: вырожденный прообраз, пропускаем
        continue
    Yr = Q[1]
    K = QQ.extension(Rx([-(Xr / a6), 0, 1]), 'w')
    w = K.gen()                                   # w^2 = X_raw/a6 = x0^2
    y0 = K(Yr / a6)
    E2K = EllipticCurve(K, [0, a2, 0, a4 * a0, a6 * a0**2])
    Bp = E2K([0, a0 * 6647, 1])
    S = E2K([a0 / w**2, a0 * y0 / w**3, 1]) + E2K([a0 / w**2, -a0 * y0 / w**3, 1]) - Bp - (-Bp)
    pp_checks.append([str(Q), bool(S.is_zero())])
    say("[5] численно: phi2_*(phi1^*(%s)) = O ?  %s" % (Q, S.is_zero()))
assert all(c[1] for c in pp_checks)

R['Psi_Phi'] = {
    'Psi': '(phi1_*, phi2_*) : Jac(D) -> E1 x E2  ((2,2)-изогения, deg 4)',
    'Phi': 'phi1^* + phi2^* : E1 x E2 -> Jac(D)  (deg 4)',
    'Psi_after_Phi': '[2]',
    'cross_terms_vanish': 'тождественно, прямой расчёт на дивизорах (не нужна гипотеза Hom=0)',
    'second_argument_non_isogenous': noniso,
    'numeric_checks': pp_checks,
    'status': 'proved',
}

# ====================== 6. ТОЧНЫЙ образ Jac(D)(Q) в E1(Q) x E2(Q)
# Верхнее включение (что содержится в образе):
#   Psi(Phi(E1(Q)xE2(Q))) = 2*(E1(Q) x E2(Q))  \subset  Psi(Jac(D)(Q)).
#   Плюс рациональный 2-кручёный класс [W(g2) - inf+ - inf-] in Jac(D)(Q).
# Нижнее включение (чего в образе НЕТ): локальные препятствия mod p.

# 6a. явный класс 2-кручения: Psi([W(g2) - D_inf]) = (T1, T2)
KQi = QuadraticField(-1, 'i')
E1K = EllipticCurve(KQi, [0, a4, 0, a2 * a6, a0 * a6**2])
E2K = EllipticCurve(KQi, [0, a2, 0, a4 * a0, a6 * a0**2])
RK = PolynomialRing(KQi, 'z')
z = RK.gen()
roots_g2 = [r for r, _ in (289 * z**2 - 322 * z + 289).roots()]
assert len(roots_g2) == 2
S1 = E1K(0)
S2 = E2K(0)
for b in roots_g2:
    S1 += E1K([a6 * b**2, 0, 1])
    S2 += E2K([a0 / b**2, 0, 1])
assert S1 == E1K([-a0, 0, 1]) and S2 == E2K([-a6, 0, 1])
assert M1([-a0 + r1shift, 0, 1]) == T1 and M2([-a6 + r2shift, 0, 1]) == T2
say("[6a] класс 2-кручения [W(289x^2-322x+289) - inf+ - inf-] in Jac(D)(Q)")
say("[6a]   Psi этого класса = (T1, T2) = (%s, %s)  -- ПРОВЕРЕНО над Q(i)" % (T1, T2))
say("[6a] => Psi(Jac(D)(Q)) СОДЕРЖИТ 2*(E1(Q)xE2(Q)) + <(T1,T2)>, индекс делит 16/2 = 8")

# 6b. локальные образы Psi(Jac(D)(F_p)) и исключение классов
A1M = [int(z0) for z0 in M1.ainvs()]
A2M = [int(z0) for z0 in M2.ainvs()]


def local_image(p):
    """Возвращает (#D(F_p), |E1xE2(F_p)|, Psi(Jac(D)(F_p)) как множество, локальный индекс)."""
    F = GF(p)
    F2 = GF(p**2, 'aa')
    m1 = M1.change_ring(F)
    m2 = M2.change_ring(F)

    def curves(K):
        return (EllipticCurve(K, [K(z0) for z0 in A1M]),
                EllipticCurve(K, [K(z0) for z0 in A2M]))

    def ph(P, K):
        e1, e2 = curves(K)
        if P[0] == 'inf':
            return e1(0), e2([K(r2shift), K(a0) * K(P[1]) * K(6647), 1])
        xv, yv = K(P[1]), K(P[2])
        A = e1([K(a6) * xv * xv + K(r1shift), K(a6) * yv, 1])
        B = e2(0) if xv == 0 else e2([K(a0) / xv**2 + K(r2shift), K(a0) * yv / xv**3, 1])
        return A, B

    def pts(K):
        L = []
        for xv in K:
            v = K(a6) * xv**6 + K(a4) * xv**4 + K(a2) * xv**2 + K(a0)
            if v.is_square():
                sq = v.sqrt()
                L.append(('fin', xv, sq))
                if sq != 0:
                    L.append(('fin', xv, -sq))
        return L

    Dp = pts(F) + [('inf', 1), ('inf', -1)]
    nD = len(Dp)
    # контроль: #D(F_p) = p + 1 - a_p(E1) - a_p(E2)
    assert nD == p + 1 - m1.trace_of_frobenius() - m2.trace_of_frobenius()
    S = set()
    n = len(Dp)
    for i in range(n):
        for j in range(i, n):
            A1p, B1p = ph(Dp[i], F)
            A2p, B2p = ph(Dp[j], F)
            S.add((A1p + A2p, B1p + B2p))
    Fpset = set(F2(xv) for xv in F)
    seen = set()
    for P in pts(F2):
        if P[1] in Fpset or (P[1], P[2]) in seen:
            continue
        Qc = ('fin', P[1]**p, P[2]**p)
        seen.add((P[1], P[2]))
        seen.add((Qc[1], Qc[2]))
        A1p, B1p = ph(P, F2)
        A2p, B2p = ph(Qc, F2)
        A = A1p + A2p
        B = B1p + B2p
        A = m1(0) if A.is_zero() else m1([F(A[0]), F(A[1]), 1])
        B = m2(0) if B.is_zero() else m2([F(B[0]), F(B[1]), 1])
        S.add((A, B))
    # самоконтроль 1: S -- подгруппа
    L = list(S)
    for _ in range(300):
        u = random.choice(L)
        v = random.choice(L)
        assert (u[0] + v[0], u[1] + v[1]) in S, "S не подгруппа при p=%s" % p
    # самоконтроль 2: 2*(E1 x E2)(F_p) \subset S
    for Q in m1.points()[:40]:
        for Rq in m2.points()[:40]:
            assert (2 * Q, 2 * Rq) in S, "2*(E1xE2)(F_p) не в S при p=%s" % p
    N = m1.cardinality() * m2.cardinality()
    return nD, N, S, ZZ(N) // len(S), m1, m2


classes = [(i, j, k, l) for i in [0, 1] for j in [0, 1] for k in [0, 1] for l in [0, 1]]
survivors = set(classes)
local_rows = []
Dboth = D1 * D2
for p in prime_range(5, 100):
    if Dboth % p == 0:
        continue
    nD, N, S, idx, m1, m2 = local_image(p)
    ok = set()
    for (i, j, k, l) in classes:
        if (m1(i * P1 + j * T1), m2(k * P2 + l * T2)) in S:
            ok.add((i, j, k, l))
    survivors &= ok
    local_rows.append({'p': int(p), 'num_D_Fp': int(nD), 'card_E1xE2_Fp': int(N),
                       'card_image': int(len(S)), 'local_index': int(idx),
                       'surviving_classes': sorted([list(c) for c in ok])})
    say("[6b] p=%-3d #D(F_p)=%-4d |E1xE2(F_p)|=%-7d |Psi(J(F_p))|=%-7d индекс=%s  выжило классов=%d"
        % (p, nD, N, len(S), idx, len(ok)))
    say("[6b]      накопленное пересечение: %s" % sorted(survivors))

# 6c'. НЕЗАВИСИМАЯ теоретическая проверка локальных индексов.
# Вейерштрассовы точки D над Q(i): x = +-alpha (alpha^2=-49/529), +-beta, +-beta_bar
#   (beta, beta_bar -- корни g2 = 289x^2-322x+289, disc(g2) = -480^2).
# И phi1, и phi2 склеивают их в ОДНИ И ТЕ ЖЕ три пары {alpha,-alpha}, {beta,-beta},
# {beta_bar,-beta_bar}. Отсюда K = ker Psi = {0, [W(g1)-Dinf], [{beta,-beta}-Dinf],
# [{bb,-bb}-Dinf]} порядка 4; Фробениус переставляет два последних ровно тогда,
# когда beta не лежит в F_p, то есть когда -1 не квадрат mod p.
# Значит |K(F_p)| = 4 при p = 1 mod 4 и 2 при p = 3 mod 4, а локальный индекс
# [E1xE2(F_p) : Psi(J(F_p))] = |K(F_p)| (так как |J(F_p)| = |E1(F_p)|*|E2(F_p)| по Тейту).
pred_ok = True
for row in local_rows:
    pred = 4 if row['p'] % 4 == 1 else 2
    row['predicted_local_index'] = pred
    if pred != row['local_index']:
        pred_ok = False
    say("[6c'] p=%-3d локальный индекс вычислен=%s, предсказан теорией=%s  %s"
        % (row['p'], row['local_index'], pred, "OK" if pred == row['local_index'] else "РАСХОЖДЕНИЕ"))
say("[6c'] независимое теоретическое предсказание локальных индексов совпало на всех %d простых: %s"
    % (len(local_rows), pred_ok))
assert pred_ok

say("[6c] ИТОГ: классы (i*P1+j*T1, k*P2+l*T2) mod 2(E1(Q)xE2(Q)), совместимые со ВСЕМИ p: %s"
    % sorted(survivors))
assert (0, 0, 0, 0) in survivors and (0, 1, 0, 1) in survivors
exact = (survivors == {(0, 0, 0, 0), (0, 1, 0, 1)})
if exact:
    say("[6c] нижнее и верхнее включения совпали => образ определён ТОЧНО:")
    say("[6c]   Psi(Jac(D)(Q)) = 2*(E1(Q) x E2(Q)) + <(T1,T2)> ,  индекс = 8")
    say("[6c]   т.е. (m1 P1 + e1 T1, m2 P2 + e2 T2) лежит в образе  <=>  m1 = m2 = 0 mod 2 и e1 = e2")
else:
    say("[6c] локальные условия НЕ довели до точного ответа; доказано лишь: индекс делит %s"
        % (16 // 2))

R['image_of_Jac'] = {
    'Psi': '(phi1_*, phi2_*) : Jac(D) -> E1 x E2',
    'lower_containment': {
        'statement': 'Psi(Jac(D)(Q)) содержит 2*(E1(Q)xE2(Q)) + <(T1,T2)>',
        'reason_1': 'Psi o Phi = [2], крест-члены обнулены тождественно',
        'reason_2': 'класс 2-кручения [W(289x^2-322x+289) - inf+ - inf-] in Jac(D)(Q) даёт (T1,T2)',
        'verified_over': 'Q(i)',
        'status': 'proved',
    },
    'upper_containment': {
        'method': 'редукция mod p: z in Jac(D)(Q) => Psi(z) mod p in Psi(Jac(D)(F_p))',
        'primes_used': [row['p'] for row in local_rows],
        'per_prime': local_rows,
        'self_checks': ['S -- подгруппа (300 случайных пар на каждое p)',
                        '2*(E1xE2)(F_p) содержится в S',
                        '#D(F_p) = p + 1 - a_p(E1) - a_p(E2)'],
        'surviving_classes': sorted([list(c) for c in survivors]),
        'status': 'proved_software',
    },
    'kernel_of_Psi': {
        'order': 4,
        'description': ('K = {0, [W(529x^2+49)-Dinf], [{beta,-beta}-Dinf], [{bbar,-bbar}-Dinf]}, '
                        'beta -- корень 289x^2-322x+289'),
        'K_over_Q': 'Z/2 (только [W(529x^2+49)-Dinf] рационален)',
        'local_index_formula': '|K(F_p)| = 4 если p = 1 mod 4, иначе 2',
        'formula_matches_all_computed_local_indices': bool(pred_ok),
        'status': 'proved (сверено с брутфорсом на всех использованных p)',
    },
    'exact_image_determined': bool(exact),
    'index_in_product': 8 if exact else None,
    'index_divides': 8,
    'index_divides_unconditionally': 16,
    'membership_rule': ('(m1*P1 + e1*T1, m2*P2 + e2*T2) in Psi(Jac(D)(Q))  <=>  '
                        'm1 = 0 mod 2, m2 = 0 mod 2, e1 = e2') if exact else None,
    'warning': ('E1(Q) x E2(Q) -- НАДмножество образа. Брать всё произведение как область решета '
                'безопасно (можно недоисключить, но нельзя исключить лишнее). Обратное '
                '(считать произведение равным образу) НЕ обосновано и здесь опровергнуто: индекс 8.'),
}

# ====================================== 7. КОНТРОЛЬ: четыре известные точки D
# базовая точка Абеля-Якоби: inf+ ;  z |-> [z - inf+]
A1pt = M1([r1shift, a6 * 2023, 1])            # phi1((0, 2023))
B2pt = M2([r2shift, a0 * 6647, 1])            # phi2(inf+)
say("[7] phi1((0,2023)) = %s ; это -2*P1 ? %s" % (A1pt, A1pt == -2 * P1))
say("[7] phi2(inf+)     = %s ; это -2*P2 ? %s" % (B2pt, B2pt == -2 * P2))
assert A1pt == -2 * P1 and B2pt == -2 * P2


def coords1(Q):
    """Q = m*P1 + e*T1 -> (m, e); поиск в разумном диапазоне."""
    if Q.is_zero():
        return (0, 0)
    for m in range(-8, 9):
        for e in [0, 1]:
            if m * P1 + e * T1 == Q:
                return (m, e)
    return None


def coords2(Q):
    if Q.is_zero():
        return (0, 0)
    for m in range(-8, 9):
        for e in [0, 1]:
            if m * P2 + e * T2 == Q:
                return (m, e)
    return None


known = []
for (xv, yv, tag) in [(QQ(0), QQ(2023), 't=1'), (QQ(0), QQ(-2023), 't=1')]:
    Q1 = M1([r1shift, a6 * yv, 1])
    Q2 = M2(0)
    c1, c2 = coords1(Q1), coords2(Q2)
    # условие принадлежности: Psi([z - inf+]) = (Q1, Q2 - B2pt)
    d2 = coords2(Q2 - B2pt)
    passes = (c1[0] % 2 == 0) and (d2[0] % 2 == 0) and (c1[1] == d2[1])
    known.append({'point': ['0', str(yv)], 't': tag, 'phi1': str(Q1), 'phi2': 'O',
                  'phi1_coords_(m,e)': list(c1), 'Psi_minus_base_phi2_coords': list(d2),
                  'passes_sieve_condition': bool(passes)})
    say("[7] z=(0,%s): phi1 = %s*P1 + %s*T1 ; [z-inf+] -> phi2-часть = %s*P2 + %s*T2 ; условие решета: %s"
        % (yv, c1[0], c1[1], d2[0], d2[1], passes))
for sgn in [1, -1]:
    Q1 = M1(0)
    Q2 = M2([r2shift, sgn * a0 * 6647, 1])
    c1 = coords1(Q1)
    d2 = coords2(Q2 - B2pt)
    passes = (c1[0] % 2 == 0) and (d2[0] % 2 == 0) and (c1[1] == d2[1])
    known.append({'point': ['infty', '%s*6647' % sgn], 't': 't=-1', 'phi1': 'O', 'phi2': str(Q2),
                  'phi2_coords_(m,e)': list(coords2(Q2)),
                  'Psi_minus_base_phi2_coords': list(d2),
                  'passes_sieve_condition': bool(passes)})
    say("[7] z=inf_%s: phi2 = %s*P2 + %s*T2 ; [z-inf+] -> phi2-часть = %s*P2 + %s*T2 ; условие решета: %s"
        % ('+' if sgn > 0 else '-', coords2(Q2)[0], coords2(Q2)[1], d2[0], d2[1], passes))
allpass = all(k['passes_sieve_condition'] for k in known)
say("[7] КОНТРОЛЬ: все четыре известные точки проходят условие решета: %s" % allpass)
assert allpass, "РЕШЕТО СЛОМАНО: известная точка не проходит -- немедленно сообщить"
R['known_points_control'] = {'points': known, 'all_pass': bool(allpass),
                             'base_point_for_Abel_Jacobi': 'inf+',
                             'status': 'proved_software'}

# ============================ 8. безусловное условие согласования двух факторов
say("[8] Для z=(x,y) in D(Q), x != 0, x != oo:  X1*X2 = a6 x^2 * a0/x^2 = a0*a6 = %s = %s"
    % (a0 * a6, factor(a0 * a6)))
say("[8]   в минимальных координатах: (X1min - %s)*(X2min - %s) = %s" % (r1shift, r2shift, a0 * a6))
say("[8]   и Y1/Y2 = a6 x^3 / a0, откуда x = a0*Y1/(X1raw*Y2)")
say("[8] Это тождество НЕ требует знания индекса образа -- оно верно поточечно.")
say("[8] Вырождения: x=0 <=> phi2(z)=O (и X1raw=0); x=oo <=> phi1(z)=O (и X2raw=0).")

R['compatibility'] = {
    'X1raw_times_X2raw': str(a0 * a6),
    'X1raw_times_X2raw_factored': str(factor(a0 * a6)),
    'in_minimal_coords': '(X1min - %s)*(X2min - %s) = %s' % (r1shift, r2shift, a0 * a6),
    'Y1_over_Y2': 'a6*x^3/a0',
    'x_recovered': 'x = a0*Y1/(X1raw*Y2)',
    'degenerate_cases': {'x=0': 'phi2(z)=O, X1raw=0', 'x=infty': 'phi1(z)=O, X2raw=0'},
    'note': 'поточечное тождество, не зависит от индекса образа Jac(D)(Q)',
    'status': 'proved',
}

R['sieve_input'] = {
    'group_to_sieve_over': 'E1(Q) x E2(Q) = (<P1> (+) <T1>) x (<P2> (+) <T2>)',
    'generators': {'P1': str(P1), 'T1': str(T1), 'P2': str(P2), 'T2': str(T2)},
    'heights': {'hhat_P1': str(hP1), 'hhat_P2': str(hP2)},
    'constraint_A_unconditional': '(X1raw)*(X2raw) = %s' % (a0 * a6),
    'constraint_B_from_index': ('для z in D(Q): phi1(z) = m1 P1 + e1 T1, phi2(z) = m2 P2 + e2 T2 '
                                '=> m1 = 0 mod 2, m2 = 0 mod 2, e1 = e2') if exact else None,
    'constraint_B_status': 'proved_software (локальные препятствия mod p, §6)' if exact else 'не установлено',
    'constraint_B_reduction_factor': 8 if exact else 1,
    'control_passed': bool(allpass),
    'must_not_do': ('нельзя исключать кандидата на основании предположения, что образ Jac(D)(Q) '
                    'совпадает со всем E1(Q)xE2(Q); и нельзя объявлять неудалённого кандидата '
                    'рациональной точкой'),
}

(OUT / 'basis_E2.json').write_text(
    json.dumps(R, indent=2, ensure_ascii=False, default=jdefault) + '\n')
(OUT / 'basis_E2.log').write_text("\n".join(LOG) + "\n")
say("записано: %s , %s" % (OUT / 'basis_E2.json', OUT / 'basis_E2.log'))
