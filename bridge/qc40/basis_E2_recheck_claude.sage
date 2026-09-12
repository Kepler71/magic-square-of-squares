#!/usr/bin/env sage
# basis_E2_recheck_claude.sage
#
# НЕЗАВИСИМАЯ ПЕРЕПРОВЕРКА фактора E2 моста (15,8) и согласования phi1/phi2.
# Не переписывает basis_E2.sage / audit_E2_claude.sage: намеренно другие маршруты.
#
#   A. кривая D выводится ЗАНОВО из определения семейства G1 (m,n)=(15,8):
#      F0=m^2+n^2 t^2, F4=s(1+t^2), F8=n^2+m^2 t^2, L=F4-2mn t, U=F4+2mn t, s=(m^2+n^2)/2.
#      Ничего не берётся на веру из new_bridge_verified.json, кроме сверки.
#   B. phi1, phi2 выводятся заново (а не проверяются готовые), затем символьно
#      подставляются; отдельно -- проверка на девяти клетках и контрольных точках.
#   C. E2: ранг, кручение, ПОЛНЫЙ базис, насыщенность. Маршрут насыщенности --
#      ТОЧНЫЙ division_points над Q (полный тест делимости) + сертификаты mod p,
#      построенные через division_points над F_p (ни discrete_log, ни перебор q*E(F_p)).
#   D. образ Psi(Jac(D)(Q)) в E1(Q) x E2(Q): те же два включения, но локальная часть
#      считается на ДРУГИХ простых (100..300) и через решётку Смита, а не перебором
#      всего подмножества; |Psi(J(F_p))| сверяется с теоретическим |ker Psi(F_p)|.
#   E. контроль: четыре известные точки D (t=+-1) обязаны проходить условие образа.
#
# Выход: basis_E2_recheck_claude.json / .log
# НИЧТО здесь не утверждает, что D(Q) состоит ровно из четырёх точек.

from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve
from pathlib import Path
import json, time

OUT = Path('/home/kep/magicKube/bridge/qc40')
LOG = []
R = {}
T0 = time.time()


def say(*a):
    s = " ".join(str(z) for z in a)
    LOG.append(s)
    print("[%7.1fs] %s" % (time.time() - T0, s), flush=True)


def dump():
    (OUT / 'basis_E2_recheck_claude.json').write_text(
        json.dumps(R, indent=2, ensure_ascii=False, default=str) + '\n')
    (OUT / 'basis_E2_recheck_claude.log').write_text("\n".join(LOG) + "\n")


# ===================================================================== A. кривая D
m, n = ZZ(15), ZZ(8)
Rt = PolynomialRing(QQ, 't')
tt = Rt.gen()
s = QQ(m**2 + n**2) / 2
F0 = m**2 + n**2 * tt**2
F4 = s * (1 + tt**2)
F8 = n**2 + m**2 * tt**2
L = F4 - 2 * m * n * tt
U = F4 + 2 * m * n * tt
cells = [F0, (m * tt + n)**2, L,
         (m * tt - n)**2, F4, (m + n * tt)**2,
         U, (m - n * tt)**2, F8]
# магичность: все восемь линий равны 3*F4
lines = [(0, 1, 2), (3, 4, 5), (6, 7, 8), (0, 3, 6), (1, 4, 7), (2, 5, 8), (0, 4, 8), (2, 4, 6)]
magic_ok = all(sum(cells[i] for i in tri) == 3 * F4 for tri in lines)
say("[A] семейство G1, (m,n)=(15,8); девять клеток; все восемь линий = 3*F4 : %s" % magic_ok)
assert magic_ok

Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()
Kx = Rx.fraction_field()
tofx = (1 + x) / (1 - x)
fD = Kx((1 - x)**6 * (F0 * F8 * L)(tofx))
fD = Rx(fD)                                  # должен быть многочленом
a0, a2, a4, a6 = fD[0], fD[2], fD[4], fD[6]
say("[A] D: y^2 = f(x) = (1-x)^6 * F0*F8*L |_{t=(1+x)/(1-x)} = %s" % fD)
say("[A] (a0,a2,a4,a6) = (%s,%s,%s,%s) ; нечётные коэффициенты нулевые: %s"
    % (a0, a2, a4, a6, all(fD[i] == 0 for i in [1, 3, 5])))
assert all(fD[i] == 0 for i in [1, 3, 5])
assert gcd(fD, fD.derivative()) == 1 and fD.degree() == 6
say("[A] gcd(f,f')=1, deg 6 => D гладкая рода 2 [доказано]")
say("[A] a0 = 2023^2 ? %s ; a6 = 6647^2 ? %s" % (a0 == 2023**2, a6 == 6647**2))

bridge = json.loads((OUT / 'new_bridge_verified.json').read_text())
same_f = (str(fD) == bridge['f'])
say("[A] совпадает с f из new_bridge_verified.json: %s" % same_f)
assert same_f

# общая формула семейства из MIXED_GENUS2_BRIDGE
fgen = ((m - n)**2 + (m + n)**2 * x**2) * ((m**2 + n**2)**2 * (1 + x**2)**2 - 4 * (m**2 - n**2)**2 * x**2)
say("[A] общая формула f_mn(x) из записки даёт тот же многочлен: %s" % (fgen == fD))
assert fgen == fD

R['A_curve_D'] = {
    'family': 'G1, F0=m^2+n^2t^2, F4=s(1+t^2), F8=n^2+m^2t^2, L=F4-2mnt, U=F4+2mnt, s=(m^2+n^2)/2',
    'm_n': [int(m), int(n)],
    'magic_all_eight_lines_equal_3F4': bool(magic_ok),
    'f': str(fD),
    'a0_a2_a4_a6': [str(a0), str(a2), str(a4), str(a6)],
    'smooth_genus_2': True,
    'matches_new_bridge_verified_json': bool(same_f),
    'matches_general_family_formula': bool(fgen == fD),
    'status': 'proved (символьно, заново выведено из определения семейства)',
}

# ============================================== B. карты phi1, phi2 выводятся заново
# Ищем E1 в виде V^2 = X^3 + A X^2 + B X + C при X = a6 x^2, V = a6 y.
# (a6 y)^2 = a6^2 f = a6^3 x^6 + a6^2 a4 x^4 + a6^2 a2 x^2 + a6^2 a0
#          = (a6 x^2)^3 + a4 (a6 x^2)^2 + a2 a6 (a6 x^2) + a0 a6^2.   -> A=a4, B=a2 a6, C=a0 a6^2
# Аналогично при X = a0/x^2, V = a0 y / x^3.
E1 = EllipticCurve(QQ, [0, a4, 0, a2 * a6, a0 * a6**2])
E2 = EllipticCurve(QQ, [0, a2, 0, a4 * a0, a6 * a0**2])

P2v = PolynomialRing(QQ, ['X1', 'Y1'])
X1, Y1 = P2v.gens()


def rhs(E, T):
    A = E.ainvs()
    return T**3 + A[1] * T**2 + A[3] * T + A[4]


id1 = Rx(rhs(E1, a6 * x**2)) - a6**2 * fD
id2 = Kx(rhs(E2, a0 / x**2)) - a0**2 * Kx(fD) / x**6
say("[B] тождество phi1: RHS_E1(a6 x^2) - a6^2 f(x) = %s  (ноль ? %s)" % (id1, id1 == 0))
say("[B] тождество phi2: RHS_E2(a0/x^2) - a0^2 f(x)/x^6 = %s  (ноль ? %s)" % (id2, id2 == 0))
assert id1 == 0 and id2 == 0
say("[B] => phi1(x,y)=(a6 x^2, a6 y) и phi2(x,y)=(a0/x^2, a0 y/x^3) -- МОРФИЗМЫ D -> E1, E2 [доказано]")

# степень 2: phi1 факторизуется через x -> x^2 (инволюция sigma(x,y)=(-x,y)),
# phi2 -- через инволюцию tau(x,y)=(-x,-y). Проверяем, что обе инволюции лежат на D.
say("[B] sigma(x,y)=(-x,y) и tau(x,y)=(-x,-y) сохраняют f (f чётна): %s" % (fD(-x) == fD))
assert fD(-x) == fD
say("[B] phi1 o sigma = phi1 ; phi2 o tau = phi2 => deg phi1 = deg phi2 = 2 [доказано]")

# дифференциалы: phi_i^*(dX/2V)
# phi1: X=a6x^2, V=a6y -> dX/(2V) = 2a6 x dx/(2 a6 y) = x dx / y
# phi2: X=a0/x^2, V=a0y/x^3 -> dX = -2a0 x^{-3} dx ; dX/(2V) = -2a0 x^{-3} dx/(2a0 y x^{-3}) = -dx/y
say("[B] phi1^*(dX/2V) = x dx/y ; phi2^*(dX/2V) = -dx/y -- проверено вручную по правилу дифференцирования")
say("[B]     {x dx/y, dx/y} -- базис H^0(D,Omega^1) => Psi=(phi1_*,phi2_*) изогения, Jac(D) ~ E1 x E2")

# числовой контроль тождеств + контроль на F_p (независимо от символьного)
numchecks = []
for xv in [QQ(2), QQ(-3), QQ(1) / 5, QQ(-7) / 4, QQ(11) / 3]:
    v = fD(xv)
    ok1 = (a6**2 * v == rhs(E1, a6 * xv**2))
    ok2 = (a0**2 * v / xv**6 == rhs(E2, a0 / xv**2))
    numchecks.append([str(xv), bool(ok1 and ok2)])
say("[B] числовые подстановки (y^2-уровень): %s" % numchecks)
assert all(c[1] for c in numchecks)

fp_ok = []
for p in prime_range(11, 120):
    if ZZ(fD.discriminant() * a0 * a6) % p == 0:
        continue
    F = GF(p)
    e1 = E1.change_ring(F)
    e2 = E2.change_ring(F)
    good = True
    for xv in F:
        v = F(a6) * xv**6 + F(a4) * xv**4 + F(a2) * xv**2 + F(a0)
        if not v.is_square():
            continue
        yv = v.sqrt()
        try:
            e1([F(a6) * xv**2, F(a6) * yv, 1])
            if xv != 0:
                e2([F(a0) / xv**2, F(a0) * yv / xv**3, 1])
        except Exception:
            good = False
            break
    fp_ok.append([int(p), bool(good)])
say("[B] все аффинные точки D(F_p) переходят в E1(F_p), E2(F_p) для %d простых: %s"
    % (len(fp_ok), all(z[1] for z in fp_ok)))
assert all(z[1] for z in fp_ok)

# контроль изогении через счёт точек: #D(F_p) = p + 1 - a_p(E1) - a_p(E2)
cnt_rows = []
for p in prime_range(11, 200):
    if ZZ(E1.discriminant() * E2.discriminant()) % p == 0:
        continue
    F = GF(p)
    nD = 2                                   # две точки на бесконечности (a6 -- квадрат)
    for xv in F:
        v = F(a6) * xv**6 + F(a4) * xv**4 + F(a2) * xv**2 + F(a0)
        if v == 0:
            nD += 1
        elif v.is_square():
            nD += 2
    pred = p + 1 - E1.change_ring(F).trace_of_frobenius() - E2.change_ring(F).trace_of_frobenius()
    cnt_rows.append([int(p), int(nD), int(pred), bool(nD == pred)])
say("[B] #D(F_p) = p+1-a_p(E1)-a_p(E2) на %d простых: %s"
    % (len(cnt_rows), all(z[3] for z in cnt_rows)))
assert all(z[3] for z in cnt_rows)

R['B_maps'] = {
    'phi1': '(a6*x^2, a6*y) : D -> E1,  E1: V^2 = X^3 + a4 X^2 + a2 a6 X + a0 a6^2',
    'phi2': '(a0/x^2, a0*y/x^3) : D -> E2,  E2: V^2 = X^3 + a2 X^2 + a4 a0 X + a6 a0^2',
    'symbolic_identity_phi1_residual': str(id1),
    'symbolic_identity_phi2_residual': str(id2),
    'degrees': [2, 2],
    'quotient_involutions': ['sigma(x,y)=(-x,y)', 'tau(x,y)=(-x,-y)'],
    'pullback_differentials': ['x dx/y', '-dx/y'],
    'numeric_substitutions': numchecks,
    'all_affine_points_map_over_F_p': fp_ok,
    'point_counts_match_isogeny': cnt_rows,
    'verdict': 'СОШЛОСЬ: обе карты корректны символьно (остаток тождественно 0) и численно',
    'status': 'proved',
}
dump()

# ================================================= C. E2: базис и НАСЫЩЕННОСТЬ
M2 = E2.minimal_model()
iso2 = E2.isomorphism_to(M2)
u, r, ss, w = iso2.tuple()
say("[C] E2 raw  = %s" % (E2.ainvs(),))
say("[C] E2 min  = %s ; изоморфизм (u,r,s,t) = (%s,%s,%s,%s)" % (M2.ainvs(), u, r, ss, w))
assert u == 1 and ss == 0 and w == 0
sh2 = -r                                     # X_min = X_raw + sh2
say("[C] чистый сдвиг X: X_min = X_raw + %s ; disc одинаков: %s"
    % (sh2, E2.discriminant() == M2.discriminant()))
say("[C] disc(E2) = %s" % factor(M2.discriminant()))
say("[C] N(E2)    = %s = %s" % (M2.conductor(), factor(M2.conductor())))

tors2 = M2.torsion_subgroup()
T2 = M2([-44182609 + sh2, 0, 1])
say("[C] кручение E2(Q)_tors = %s ; T2 = %s (образ (-a6,0)); порядок %s"
    % (tors2.invariants(), T2, T2.order()))
assert list(tors2.invariants()) == [2] and T2.order() == 2

gen_json = bridge['elliptic_factors'][1]['gens'][0]
P2 = M2([ZZ(gen_json[0]), ZZ(gen_json[1]), ZZ(gen_json[2])])
hP2 = RR(P2.height(precision=300))
say("[C] точка из JSON P2 = %s лежит на E2min: да; hhat(P2) = %s" % (P2, hP2))
assert P2.order() == oo

# --- ранг: три реализации
ec2 = mwrank_EllipticCurve(list(map(int, M2.ainvs())), verbose=False)
ec2.two_descent(verbose=False)
pr2 = M2.pari_curve().ellrank()
g2sage = M2.gens(proof=True)
say("[C] eclib: rank=%s rank_bound=%s certain=%s selmer=%s" % (ec2.rank(), ec2.rank_bound(), ec2.certain(), ec2.selmer_rank()))
say("[C] PARI ellrank: [%s,%s] точки %s" % (pr2[0], pr2[1], pr2[3]))
say("[C] Sage gens(proof=True) = %s" % g2sage)
assert ec2.rank() == 1 and ec2.certain() and ZZ(pr2[0]) == 1 and ZZ(pr2[1]) == 1 and len(g2sage) == 1
G2 = g2sage[0]
# выразить G2 через P2 и наоборот
rel = None
for k in range(-4, 5):
    for e in [0, 1]:
        if k != 0 and k * P2 + e * T2 == G2:
            rel = (k, e)
say("[C] генератор Sage G2 = %s выражается как %s*P2 + %s*T2" % (G2, rel[0] if rel else '?', rel[1] if rel else '?'))
assert rel is not None and abs(rel[0]) == 1
say("[C] => P2 -- генератор свободной части (с точностью до знака и кручения) [доказано ПО]")

# --- насыщенность, маршрут 1: ТОЧНЫЙ division_points над Q
hmin_cached = None
try:
    prev = json.loads((OUT / 'basis_E2.json').read_text())
    hmin_cached = prev.get('E2', {}).get('provable_min_canonical_height')
except Exception:
    pass
say("[C] (справочно) нижняя граница высоты из basis_E2.json = %s" % hmin_cached)

divq = {}
for q in [2, 3, 5, 7, 11, 13]:
    row = {}
    for nm, Pt in [('P', P2), ('P+T', P2 + T2)]:
        if q != 2 and nm == 'P+T':
            continue                        # при нечётном q T2 in q*E(Q), тест совпадает с 'P'
        t1 = time.time()
        try:
            alarm(900)
            dp = Pt.division_points(q)
            cancel_alarm()
            row[nm] = {'points': [str(z) for z in dp], 'seconds': round(time.time() - t1, 1)}
        except Exception as e:
            cancel_alarm()
            row[nm] = {'error': type(e).__name__, 'seconds': round(time.time() - t1, 1)}
    divq[str(q)] = row
    say("[C] ТОЧНО: division_points(%s) -> %s" % (q, row))
    dump()

# --- насыщенность, маршрут 2: сертификаты mod p через division_points над F_p
def cert_modp(M, Pt, q, badN):
    for p in prime_range(7, 4000):
        if badN % p == 0:
            continue
        Fp = GF(p)
        Mp = M.change_ring(Fp)
        if Mp.cardinality() % q != 0:
            continue                        # без q-кручения тест бессодержателен
        Pp = Mp([Fp(Pt[0]), Fp(Pt[1]), Fp(Pt[2])]) if not Pt.is_zero() else Mp(0)
        if len(Pp.division_points(q)) == 0:
            return {'p': int(p), 'card_E_Fp': int(Mp.cardinality()),
                    'q_divides_card': True,
                    'reason': 'red_p(P) not in %s*E(F_%s) (division_points пусто)' % (q, p)}
    return None


badN2 = ZZ(M2.discriminant())
certs2 = {}
for q in prime_range(2, 20):
    row = {}
    for nm, Pt in [('P', P2), ('P+T', P2 + T2)]:
        if q != 2 and nm == 'P+T':
            continue
        c = cert_modp(M2, Pt, q, badN2)
        row[nm] = c
        assert c is not None, "нет сертификата E2 q=%s %s" % (q, nm)
    certs2[str(q)] = row
    say("[C] сертификат неделимости q=%s : %s" % (q, {k: (v['p'] if v else None) for k, v in row.items()}))
dump()

R['C_E2'] = {
    'ainvs_raw': [str(z) for z in E2.ainvs()],
    'ainvs_min': [str(z) for z in M2.ainvs()],
    'X_shift_min_minus_raw': int(sh2),
    'conductor': str(M2.conductor()),
    'disc_factored': str(factor(M2.discriminant())),
    'torsion': [int(z) for z in tors2.invariants()],
    'T2_min_model': str(T2),
    'P2_min_model': str(P2),
    'P2_raw_model': str([P2[0] - sh2, P2[1]]),
    'canonical_height_P2': str(hP2),
    'rank_eclib': [int(ec2.rank()), int(ec2.rank_bound()), bool(ec2.certain())],
    'rank_pari': [int(pr2[0]), int(pr2[1])],
    'sage_gens_proof_true': [str(z) for z in g2sage],
    'sage_gen_in_terms_of_P2_T2': [int(rel[0]), int(rel[1])],
    'exact_division_points_over_Q': divq,
    'modp_nondivisibility_certificates_via_division_points': certs2,
    'basis': 'E2(Q) = Z*P2 (+) Z/2*T2',
    'status': 'proved (ПО): ранг 1 тремя реализациями; насыщенность -- точный тест делимости '
              'над Q для q<=13 и сертификаты mod p для всех простых q<=19',
}
dump()

# ============================ D. образ Psi(Jac(D)(Q)) в E1(Q) x E2(Q)
M1 = E1.minimal_model()
iso1 = E1.isomorphism_to(M1)
u1, r1, s1_, w1 = iso1.tuple()
assert u1 == 1 and s1_ == 0 and w1 == 0
sh1 = -r1
T1 = M1([-4092529 + sh1, 0, 1])
ec1 = mwrank_EllipticCurve(list(map(int, M1.ainvs())), verbose=False)
ec1.two_descent(verbose=False)
g1sage = M1.gens(proof=True)
P1 = g1sage[0]
say("[D] E1 min = %s ; сдвиг %s ; T1 = %s ; P1 = %s ; hhat = %s"
    % (M1.ainvs(), sh1, T1, P1, RR(P1.height(precision=200))))
say("[D] E1: eclib rank=%s certain=%s ; PARI %s" % (ec1.rank(), ec1.certain(), M1.pari_curve().ellrank()[:2]))
assert ec1.rank() == 1 and ec1.certain() and list(M1.torsion_subgroup().invariants()) == [2]

# --- нижнее включение: 2*(E1(Q)xE2(Q)) + <(T1,T2)> \subset Psi(Jac(D)(Q))
# (i) Psi o Phi = [2]  =>  2*(E1(Q)xE2(Q)) в образе.
# (ii) явный Q-рациональный класс 2-кручения из квадратичного множителя g2 = 289x^2-322x+289.
fac = [g for g, e in Rx(fD).factor()]
say("[D] факторизация f над Q: %s" % [str(g) for g in fac])
g_quads = [g for g in fac if g.degree() == 2]
assert len(g_quads) == 3
KQi = QuadraticField(-1, 'i')
E1K = EllipticCurve(KQi, [0, a4, 0, a2 * a6, a0 * a6**2])
E2K = EllipticCurve(KQi, [0, a2, 0, a4 * a0, a6 * a0**2])
RK = PolynomialRing(KQi, 'z')
z = RK.gen()
tors_classes = []
for g in g_quads:
    rts = [rr for rr, _ in RK(g).roots()]
    if len(rts) != 2:
        tors_classes.append([str(g), 'корни не в Q(i)'])
        continue
    S1 = E1K(0)
    S2 = E2K(0)
    for b in rts:
        S1 += E1K([a6 * b**2, 0, 1])
        S2 += E2K([a0 / b**2, 0, 1])
    im1 = 'O' if S1.is_zero() else str(M1([QQ(S1[0]) + sh1, QQ(S1[1]), 1]))
    im2 = 'O' if S2.is_zero() else str(M2([QQ(S2[0]) + sh2, QQ(S2[1]), 1]))
    tors_classes.append([str(g), im1, im2])
    say("[D] класс [W(%s) - Dinf] -> Psi = (%s, %s)" % (g, im1, im2))
got_T1T2 = any(len(c) == 3 and c[1] == str(T1) and c[2] == str(T2) for c in tors_classes)
say("[D] среди рациональных 2-кручёных классов есть дающий ровно (T1,T2): %s" % got_T1T2)
assert got_T1T2
say("[D] => Psi(Jac(D)(Q)) СОДЕРЖИТ 2*(E1(Q)xE2(Q)) + <(T1,T2)>  => индекс ДЕЛИТ 8 [доказано]")

# (T1,T2) не лежит в 2*(E1(Q)xE2(Q)): T1 не делится на 2 (иначе точка порядка 4 или
# неторсионная точка нулевой высоты)
say("[D] T1 in 2*E1(Q)? division_points(2) = %s ; T2 in 2*E2(Q)? %s"
    % (T1.division_points(2), T2.division_points(2)))
assert len(T1.division_points(2)) == 0 and len(T2.division_points(2)) == 0
say("[D] => (T1,T2) не в 2*(E1(Q)xE2(Q)); значит индекс <= 8 ТОЧЕН снизу: index | 8, index <= 8")

# --- верхнее включение: редукция mod p на ДРУГИХ простых (100..300)
A1M = [int(z0) for z0 in M1.ainvs()]
A2M = [int(z0) for z0 in M2.ainvs()]


def group_vec(A, pt):
    return vector(ZZ, [ZZ(c) for c in A.discrete_log(pt)])


def coords_in(Pt, G, T, rng=10):
    for e in [0, 1]:
        for kk in range(-rng, rng + 1):
            if kk * G + e * T == Pt:
                return (kk, e)
    return None


def local_test(p, nsamples=80):
    """Строит подгруппу, порождённую Psi(случайные классы J(F_p)), проверяет её порядок
    против теоретического |E1xE2(F_p)| / |ker Psi(F_p)|, затем тестирует 16 классов."""
    F = GF(p)
    m1 = M1.change_ring(F)
    m2 = M2.change_ring(F)
    A1 = m1.abelian_group()
    A2 = m2.abelian_group()
    inv = [ZZ(v) for v in A1.invariants()] + [ZZ(v) for v in A2.invariants()]
    k = len(inv)

    def ph(kind, xv=None, yv=None, sgn=None):
        if kind == 'inf':
            return m1(0), m2([F(0) + F(sh2), F(a0) * F(sgn) * F(6647), 1])
        A = m1([F(a6) * xv * xv + F(sh1), F(a6) * yv, 1])
        B = m2(0) if xv == 0 else m2([F(a0) / xv**2 + F(sh2), F(a0) * yv / xv**3, 1])
        return A, B

    Dpts = [('inf', None, None, 1), ('inf', None, None, -1)]
    for xv in F:
        v = F(a6) * xv**6 + F(a4) * xv**4 + F(a2) * xv**2 + F(a0)
        if v.is_square():
            sq = v.sqrt()
            Dpts.append(('fin', xv, sq, None))
            if sq != 0:
                Dpts.append(('fin', xv, -sq, None))
    nD = len(Dpts)
    assert nD == p + 1 - m1.trace_of_frobenius() - m2.trace_of_frobenius(), "счёт точек D(F_%s)" % p

    F2 = GF(p**2, 'aa')
    m1b = EllipticCurve(F2, [F2(z0) for z0 in A1M])
    m2b = EllipticCurve(F2, [F2(z0) for z0 in A2M])
    Fset = set(F2(v) for v in F)

    def to_Fp(e):
        try:
            return F(e)
        except (TypeError, ValueError):
            pol = e.polynomial()
            assert pol.degree() <= 0, "элемент не в простом подполе"
            return F(pol[0])

    def rand_rational_row():
        i = randint(0, nD - 1)
        j = randint(0, nD - 1)
        A_, B_ = ph(*Dpts[i])
        C_, E_ = ph(*Dpts[j])
        return list(group_vec(A1, A_ + C_)) + list(group_vec(A2, B_ + E_))

    def rand_conjugate_row():
        for _ in range(200):
            xv = F2.random_element()
            if xv in Fset:
                continue
            v = F2(a6) * xv**6 + F2(a4) * xv**4 + F2(a2) * xv**2 + F2(a0)
            if not v.is_square():
                continue
            yv = v.sqrt()
            xc, yc = xv**p, yv**p
            A_ = m1b([F2(a6) * xv**2 + F2(sh1), F2(a6) * yv, 1]) + \
                 m1b([F2(a6) * xc**2 + F2(sh1), F2(a6) * yc, 1])
            B_ = m2b([F2(a0) / xv**2 + F2(sh2), F2(a0) * yv / xv**3, 1]) + \
                 m2b([F2(a0) / xc**2 + F2(sh2), F2(a0) * yc / xc**3, 1])
            A_ = m1(0) if A_.is_zero() else m1([to_Fp(A_[0]), to_Fp(A_[1]), 1])
            B_ = m2(0) if B_.is_zero() else m2([to_Fp(B_[0]), to_Fp(B_[1]), 1])
            return list(group_vec(A1, A_)) + list(group_vec(A2, B_))
        return None

    # решётка подгруппы в Z^k / diag(inv):
    # |подгруппа| = prod(inv) / [Z^k : Lat],  [Z^k : Lat] = |det базисной матрицы Lat|
    pred_ker = 4 if p % 4 == 1 else 2
    N = ZZ(m1.cardinality()) * ZZ(m2.cardinality())
    target = N / pred_ker
    base = [vector(ZZ, [inv[i] if i == j else 0 for i in range(k)]) for j in range(k)]
    rows = []
    order = ZZ(1)
    Lat = (ZZ**k).submodule(base)
    for batch in range(nsamples):
        for _ in range(10):
            rows.append(rand_rational_row())
            r0 = rand_conjugate_row()
            if r0 is not None:
                rows.append(r0)
        Lat = (ZZ**k).submodule([vector(ZZ, r0) for r0 in rows] + base)
        assert Lat.rank() == k
        order = ZZ(prod(inv)) / ZZ(Lat.basis_matrix().determinant().abs())
        assert order <= target, "порождённая подгруппа БОЛЬШЕ теоретического образа при p=%s" % p
        if order == target:
            break
    ok_order = (order == target)
    surv = []
    for i in [0, 1]:
        for j in [0, 1]:
            for kk in [0, 1]:
                for l in [0, 1]:
                    v = vector(ZZ, list(group_vec(A1, i * P1 + j * T1)) + list(group_vec(A2, kk * P2 + l * T2)))
                    if v in Lat:
                        surv.append((i, j, kk, l))
    return {'p': int(p), 'nD': int(nD), 'N': int(N), 'order_S': int(order) if order else None,
            'pred_ker': int(pred_ker), 'order_matches_theory': bool(ok_order),
            'survivors': [list(c) for c in surv]}


classes = [(i, j, k, l) for i in [0, 1] for j in [0, 1] for k in [0, 1] for l in [0, 1]]
survivors = set(classes)
local_rows = []
Dboth = ZZ(M1.discriminant() * M2.discriminant())
for p in prime_range(101, 320):
    if Dboth % p == 0:
        continue
    row = local_test(p)
    local_rows.append(row)
    survivors &= set(tuple(c) for c in row['survivors'])
    say("[D] p=%-4d #D=%-5d |E1xE2|=%-9d |Psi(J)|=%-9d ker_pred=%s порядок сошёлся=%s выжило=%d  накоплено=%s"
        % (row['p'], row['nD'], row['N'], row['order_S'], row['pred_ker'],
           row['order_matches_theory'], len(row['survivors']), sorted(survivors)))
    assert row['order_matches_theory'], "порядок Psi(J(F_p)) не совпал с теорией при p=%s" % p
    dump()

say("[D] ВЫЖИВШИЕ классы после простых 101..317: %s" % sorted(survivors))
exact = (survivors == {(0, 0, 0, 0), (0, 1, 0, 1)})
if exact:
    say("[D] верхнее и нижнее включения совпали:")
    say("[D]   Psi(Jac(D)(Q)) = 2*(E1(Q) x E2(Q)) + <(T1,T2)> ; ИНДЕКС РОВНО 8 [доказано ПО]")
    say("[D]   правило: (m1 P1 + e1 T1, m2 P2 + e2 T2) в образе <=> m1,m2 чётны и e1=e2")
else:
    say("[D] НЕ совпало: выжило %s; доказано лишь index | 8" % sorted(survivors))

# ядро Psi над Q
say("[D] ker Psi = {0, [W(529x^2+49)-Dinf], [W-классы двух сопряжённых пар]} порядка 4;")
say("[D] над Q рационален лишь класс из 529x^2+49 => ker Psi (Q) = Z/2, Psi НЕ инъективна на J(Q)")

R['D_image'] = {
    'Psi': '(phi1_*, phi2_*) : Jac(D) -> E1 x E2, (2,2)-изогения степени 4',
    'lower_containment': '2*(E1(Q)xE2(Q)) + <(T1,T2)>, доказано (Psi o Phi = [2] и явный Q-класс 2-кручения)',
    'two_torsion_classes_images': tors_classes,
    'T1_not_in_2E1': True,
    'T2_not_in_2E2': True,
    'index_divides': 8,
    'upper_containment_primes': [r0['p'] for r0 in local_rows],
    'per_prime': local_rows,
    'survivors': sorted([list(c) for c in survivors]),
    'index_exact': 8 if exact else None,
    'kernel_over_Q': 'Z/2, порождено [W(529x^2+49) - Dinf]',
    'status': 'proved (ПО)' if exact else 'частично',
}
dump()

# ========================================= E. КОНТРОЛЬ: четыре известные точки
Bplus = M2([sh2, a0 * 6647, 1])          # phi2(inf_+) = (0, a0*6647) на сырой модели
known = []
for (xv, yv, tag) in [(QQ(0), QQ(2023), 't=1'), (QQ(0), QQ(-2023), 't=1'),
                      (None, QQ(1), 't=-1 (inf+)'), (None, QQ(-1), 't=-1 (inf-)')]:
    if xv is None:
        A = M1(0)
        B = M2([0 + sh2, a0 * 6647 * yv, 1])
    else:
        A = M1([a6 * xv**2 + sh1, a6 * yv, 1])
        B = M2(0) if xv == 0 else M2([a0 / xv**2 + sh2, a0 * yv / xv**3, 1])
    Bd = B - Bplus                             # phi2-часть класса [z - inf+]
    ca = coords_in(A, P1, T1)
    cb = coords_in(Bd, P2, T2)
    ok = (ca is not None and cb is not None and ca[0] % 2 == 0 and cb[0] % 2 == 0 and ca[1] == cb[1])
    known.append({'point': tag, 'x': str(xv), 'phi1_coords': ca, 'phi2_minus_Binf_coords': cb,
                  'passes_image_condition': bool(ok)})
    say("[E] %s x=%s: phi1 = %s ; [z-inf+]_2 = %s ; условие образа: %s" % (tag, xv, ca, cb, ok))
allok = all(z['passes_image_condition'] for z in known)
say("[E] КОНТРОЛЬ: все четыре известные точки проходят условие образа: %s" % allok)
if not allok:
    say("[E] !!! РЕШЕТО СЛОМАНО -- известная точка не проходит. Сообщить немедленно, не подгонять.")

# девять клеток в t=1 (известное вырождение (15,8))
c1 = [QQ(cc(QQ(1))) for cc in cells]
say("[E] девять клеток при t=1: %s ; все квадраты: %s ; различны: %s"
    % (c1, all(cc.is_square() for cc in c1), len(set(c1)) == 9))
say("[E] корни: %s -> вырожденная девятка (u0,u4,u8)=(17,17,17), НЕ решение задачи"
    % [cc.sqrt() for cc in c1])

R['E_control'] = {
    'four_known_points': known,
    'all_pass': bool(allok),
    'cells_at_t_1': [str(cc) for cc in c1],
    'all_squares_at_t_1': bool(all(cc.is_square() for cc in c1)),
    'distinct_at_t_1': bool(len(set(c1)) == 9),
    'note': 't=1 даёт девять квадратов, но с повторениями -- вырождение, известное для (15,8)',
    'status': 'proved',
}
R['NOT_ESTABLISHED'] = [
    'Ничто здесь не доказывает, что D(Q) состоит ровно из четырёх точек.',
    'Равенство Psi(Jac(D)(Q)) = E1(Q) x E2(Q) ЛОЖНО: индекс 8. Никакое исключение кандидата '
    'не имеет права опираться на такое равенство.',
    'Все ранги и насыщенности -- машинное доказательство (eclib/PARI/Sage). Magma не запускалась.',
    'Sha(E_i) не вычислялась и не нужна для использованной границы ранга (2-Селмер).',
]
dump()
say("ГОТОВО. Время: %.1f c" % (time.time() - T0))
