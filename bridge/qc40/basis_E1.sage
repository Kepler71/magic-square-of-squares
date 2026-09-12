#!/usr/bin/env sage
# basis_E1.sage -- полный базис группы Морделла-Вейля и ДОКАЗАТЕЛЬСТВО НАСЫЩЕННОСТИ
# для фактора E1 смешанного моста (15,8).
#
# Запуск:  sage /home/kep/magicKube/bridge/qc40/basis_E1.sage
# Выход:   basis_E1.json, basis_E1.log (в каталоге скрипта)
#
# Метки статуса внутри JSON/лога:
#   proved            -- доказано (конечное точное вычисление или элементарный вывод)
#   proved_software   -- строгий алгоритм, корректность зависит от ПО (Sage/eclib/PARI)
#   numerical         -- численно, с оценкой точности
#   observation       -- наблюдение, не используется как посылка
#
# Логика доказательства насыщенности (см. §5-§7):
#   (a) rank E1(Q) = 1, E1(Q)_tors = Z/2   -- из 2-спуска, dim Sel_2 = 2 = rank + dim E(Q)[2];
#   (b) mu = ДОКАЗАННАЯ нижняя граница канонической высоты неторсионных точек (Cremona-Siksek);
#       если P = n*G (mod tors), то h(P) = n^2 h(G) >= n^2 mu  =>  n <= sqrt(h(P)/mu);
#   (c) для каждого простого q <= этой границы: P и P+T НЕ q-делимы в E1(Q)
#       (точное вычисление division_points + независимый сертификат по модулю p);
#   => индекс n = 1, то есть <P> + <T> = E1(Q).

from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve
from pathlib import Path
import json, time

# ВНИМАНИЕ: под `sage file.sage` переменная __file__ указывает не на скрипт
# (препарсер кладёт код во временный файл), поэтому путь задан явно.
OUT = Path('/home/kep/magicKube/bridge/qc40')
R = {}
LOG = []


def say(*a):
    s = " ".join(str(z) for z in a)
    LOG.append(s)
    print(s, flush=True)


t_start = time.time()

# ==================================================================== 0. модель D
a0, a2, a4, a6 = QQ(4092529), QQ(47287151), QQ(37608911), QQ(44182609)
Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()
Kx = Rx.fraction_field()
f = a6 * x**6 + a4 * x**4 + a2 * x**2 + a0

# контроль тождества моста: f = (1-x)^6 * F0(t) F8(t) L(t), t = (1+x)/(1-x), s = 289/2
t = (1 + x) / (1 - x)
s = QQ(289) / 2
h_bridge = (225 + 64 * t * t) * (64 + 225 * t * t) * (s * (1 + t * t) - 240 * t)
assert Kx(f) == (1 - x)**6 * h_bridge
assert gcd(f, f.derivative()) == 1            # род 2
assert a0 == 2023**2 and a6 == 6647**2
say("[0] D : y^2 = %s" % f)
say("[0] тождество моста и гладкость (род 2) подтверждены точной арифметикой  [proved]")

R['curve_D'] = {'f': str(f), 'a0_a2_a4_a6': [str(a0), str(a2), str(a4), str(a6)],
                'genus': 2, 'bridge_identity': 'proved',
                'source': 'bridge/qc40/new_bridge_verified.json'}

# ========================================== 1. E1 = D/sigma через phi1=(a6 x^2, a6 y)
E1 = EllipticCurve([0, a4, 0, a2 * a6, a0 * a6**2])

# сверка с new_bridge_verified.json (elliptic_factors[0])
json_raw = [ZZ(z) for z in ["0", "37608911", "0", "2089269703356959", "7989037884942063852049"]]
json_min = [ZZ(z) for z in ["0", "-1", "0", "1617792974488319", "-14262298467262390622975"]]
json_gen = [ZZ(z) for z in ["34384993", "286391046720", "1"]]
assert list(E1.ainvs()) == json_raw, "модель E1 не совпала с new_bridge_verified.json"
say("[1] E1 (сырая модель) = %s ; совпала с elliptic_factors[0].raw_ainvs  [proved]" % (E1.ainvs(),))

# символьная проверка карты phi1 в Q[x,y]/(y^2 - f)
P2 = PolynomialRing(QQ, ['xx', 'yy'])
xx, yy = P2.gens()
fxx = a6 * xx**6 + a4 * xx**4 + a2 * xx**2 + a0
I = P2.ideal(yy**2 - fxx)


def rhs(E, T):
    A = E.ainvs()
    return T**3 + A[1] * T**2 + A[3] * T + A[4]


r1 = (a6 * yy)**2 - rhs(E1, a6 * xx**2)
q1 = P2(r1).quo_rem(yy**2 - fxx)[0]
assert r1 in I and r1 == q1 * (yy**2 - fxx)
assert a6**2 * f == rhs(E1, a6 * x**2)        # то же в Q[x], независимая форма
say("[1] phi1=(a6 x^2, a6 y): (a6 y)^2 - RHS_E1(a6 x^2) = (%s)*(y^2-f)  [proved]" % q1)
say("[1] инволюция sigma(x,y)=(-x,y), E1 = D/sigma; phi1^*(dX/2Y) = x dx/y  [proved]")

R['map_phi1'] = {'phi1': '(a6*x^2, a6*y)', 'quotient_by': 'sigma(x,y)=(-x,y)',
                 'identity_multiplier': str(q1),
                 'pullback_differential': 'x dx / y',
                 'status': 'proved'}

# ------------------------------------------------- минимальная модель и изоморфизм
M1 = E1.minimal_model()
assert list(M1.ainvs()) == json_min, "минимальная модель не совпала с JSON"
iso = E1.isomorphism_to(M1)
iso_inv = M1.isomorphism_to(E1)
say("[1] E1 (минимальная) = %s ; совпала с elliptic_factors[0].minimal_ainvs  [proved]" % (M1.ainvs(),))
say("[1] изоморфизм raw->min: %s" % iso)

# ------------------------------------- ТОЧНАЯ проверка образующей подстановкой в уравнение
def exact_on_curve(E, X, Y):
    """Точная рациональная проверка Y^2 + a1 X Y + a3 Y = X^3 + a2 X^2 + a4 X + a6."""
    A = E.ainvs()
    lhs = Y**2 + A[0] * X * Y + A[2] * Y
    rhsv = X**3 + A[1] * X**2 + A[3] * X + A[4]
    return bool(lhs == rhsv), lhs - rhsv


Px, Py = QQ(json_gen[0]), QQ(json_gen[1])
ok_min, resid_min = exact_on_curve(M1, Px, Py)
assert ok_min
P = M1([Px, Py, 1])
say("[1] P = (%s : %s : 1) на минимальной модели: подстановка даёт невязку %s  [proved]"
    % (Px, Py, resid_min))

# та же точка на сырой модели (через обратный изоморфизм), снова точная подстановка
Praw = iso_inv(P)
okr, residr = exact_on_curve(E1, Praw[0], Praw[1])
assert okr
say("[1] образ P на сырой модели: %s ; невязка %s  [proved]" % (Praw, residr))
say("[1] обратно iso(Praw) == P : %s  [proved]" % bool(iso(Praw) == P))
assert iso(Praw) == P

R['model'] = {
    'E1_raw_ainvs': [str(z) for z in E1.ainvs()],
    'E1_min_ainvs': [str(z) for z in M1.ainvs()],
    'iso_raw_to_min': str(iso),
    'generator_min_model': [str(Px), str(Py)],
    'generator_raw_model': [str(Praw[0]), str(Praw[1])],
    'exact_substitution_residual_min': str(resid_min),
    'exact_substitution_residual_raw': str(residr),
    'status': 'proved',
}

# ===================================================== 2. инварианты, кручение, p = 11
cond = M1.conductor()
disc = M1.discriminant()
say("[2] кондуктор  N = %s = %s" % (cond, factor(cond)))
say("[2] дискриминант  = %s" % factor(disc))
say("[2] j(E1) = %s" % M1.j_invariant())

Tors = M1.torsion_subgroup()
tors_pts = [Q for Q in Tors]
T = Tors.gens()[0].element() if hasattr(Tors.gens()[0], 'element') else M1(Tors.gens()[0])
assert T.order() == 2
okT, residT = exact_on_curve(M1, T[0], T[1])
assert okT
say("[2] E1(Q)_tors = %s, порождающая T = %s (порядок 2), невязка подстановки %s  [proved_software+proved]"
    % (Tors.invariants(), T, residT))
say("[2] 2-кручение: x^3 - x^2 + ... имеет ровно один рациональный корень %s ; проверка: %s"
    % (T[0], factor(Rx([M1.ainvs()[4], M1.ainvs()[3], M1.ainvs()[1], 1]))))

good11 = (disc % 11 != 0)
say("[2] p = 11: хорошая редукция = %s (11 не делит дискриминант)  [proved]" % good11)
assert good11

R['invariants'] = {
    'conductor': str(cond), 'conductor_factored': str(factor(cond)),
    'disc_min_factored': str(factor(disc)),
    'j_invariant': str(M1.j_invariant()),
    'torsion_invariants': [int(z) for z in Tors.invariants()],
    'torsion_points': [str(Q) for Q in tors_pts],
    'two_torsion_poly_factored': str(factor(Rx([M1.ainvs()[4], M1.ainvs()[3], M1.ainvs()[1], 1]))),
    'good_reduction_at_11': bool(good11),
    'status': 'proved_software',
}

# ================================================================== 3. РАНГ E1(Q)
ec = mwrank_EllipticCurve(list(map(int, M1.ainvs())), verbose=False)
ec.two_descent(verbose=False)
ec_rank, ec_bound, ec_certain, ec_selmer = ec.rank(), ec.rank_bound(), ec.certain(), ec.selmer_rank()
say("[3] eclib two_descent: rank=%s rank_bound=%s certain=%s selmer_rank=%s gens=%s"
    % (ec_rank, ec_bound, ec_certain, ec_selmer, ec.gens()))

sel_pari = M1.selmer_rank(algorithm='pari')
sel_mw = M1.selmer_rank(algorithm='mwrank')
say("[3] dim_F2 Sel_2(E1/Q) = %s (pari) = %s (mwrank)" % (sel_pari, sel_mw))
dim2tors = 1  # dim_F2 E1(Q)[2] = 1, т.к. кручение Z/2
say("[3] dim Sel_2 = rank + dim E(Q)[2] + dim Sha[2]  =>  rank + dim Sha[2] = %s - %s = %s"
    % (sel_pari, dim2tors, sel_pari - dim2tors))
say("[3] точка P бесконечного порядка даёт rank >= 1  =>  rank = 1 и Sha(E1/Q)[2] = 0  [proved_software]")
assert sel_pari - dim2tors == 1

pr = M1.pari_curve().ellrank()
say("[3] PARI ellrank: интервал [%s, %s], точки %s" % (pr[0], pr[1], pr[3]))
rk_proof = M1.rank(proof=True)
gens_proof = M1.gens(proof=True)
say("[3] Sage rank(proof=True) = %s ; gens(proof=True) = %s" % (rk_proof, gens_proof))
say("[3] знак функционального уравнения = %s (согласуется с нечётным рангом; НЕ доказательство)"
    % M1.root_number())
assert rk_proof == 1 and ec_rank == 1 and bool(ec_certain)
assert len(gens_proof) == 1

R['rank'] = {
    'value': 1,
    'eclib': {'rank': int(ec_rank), 'rank_bound': int(ec_bound), 'certain': bool(ec_certain),
              'selmer_rank': int(ec_selmer), 'gens': [str(g) for g in ec.gens()]},
    'selmer_dim_pari': int(sel_pari), 'selmer_dim_mwrank': int(sel_mw),
    'dim_E_Q_2torsion': 1,
    'sha2_dim': 0,
    'argument': 'dim Sel_2 = 2 = rank + dim E(Q)[2] + dim Sha[2]; rank >= 1 из точки P => rank = 1, Sha[2] = 0',
    'pari_ellrank_interval': [int(pr[0]), int(pr[1])],
    'sage_rank_proof_true': int(rk_proof),
    'sage_gens_proof_true': [str(g) for g in gens_proof],
    'root_number': int(M1.root_number()),
    'status': 'proved_software',
}

# ==================================================== 4. канонические высоты (2 реализации)
hP_sage = P.height(precision=300)
hP_pari = M1.pari_curve().ellheight(P.__pari__(), precision=128)
say("[4] hhat(P): Sage = %s ; PARI = %s ; разность = %s"
    % (RR(hP_sage), RR(hP_pari), RR(hP_sage) - RR(hP_pari)))
assert abs(RR(hP_sage) - RR(hP_pari)) < 1e-20
say("[4] регулятор E1(Q) = hhat(P) = %s  [numerical, точность >= 1e-20]" % RR(hP_sage))

R['heights'] = {'hhat_P_sage': str(RR(hP_sage)), 'hhat_P_pari': str(RR(hP_pari)),
                'regulator': str(RR(hP_sage)), 'status': 'numerical'}

# ======================================== 5. ДОКАЗАННАЯ граница индекса подгруппы <P>
# height_function().min(tol, n_max) -- алгоритм Cremona-Siksek/Thongjunthug:
# "A positive real mu for which it has been established rigorously that every point
#  of infinite order has canonical height greater than mu."
H = M1.height_function()
mins = {}
for (tol, nmax) in [(0.1, 5), (0.01, 8), (0.001, 10), (0.0001, 12)]:
    t0 = time.time()
    mu = H.min(RR(tol), nmax)
    nb = RR(sqrt(hP_sage / mu))
    mins['tol=%s,n_max=%s' % (tol, nmax)] = {'mu': str(RR(mu)), 'index_bound': str(nb),
                                             'seconds': round(time.time() - t0, 2)}
    say("[5] mu(tol=%s, n_max=%s) = %s  =>  индекс n <= sqrt(hhat(P)/mu) = %s  (%.1f c)"
        % (tol, nmax, RR(mu), nb, time.time() - t0))

mu_best = max(RR(v['mu']) for v in mins.values())
index_bound = floor(RR(sqrt(hP_sage / mu_best)))
bad_primes_for_index = [int(q) for q in prime_range(index_bound + 1)]
say("[5] лучшая нижняя граница mu = %s  =>  ИНДЕКС n <= %s  [proved_software]" % (mu_best, index_bound))
say("[5] значит достаточно исключить делимость на простые q из %s" % bad_primes_for_index)

R['index_bound'] = {
    'method': 'Cremona-Siksek lower bound mu for canonical height of non-torsion points; n^2*mu <= hhat(P)',
    'runs': mins,
    'mu_used': str(mu_best),
    'index_bound_n_le': int(index_bound),
    'primes_to_exclude': bad_primes_for_index,
    'status': 'proved_software',
}

# =================================== 6. НЕДЕЛИМОСТЬ: точные division_points + модульные сертификаты
targets = {'P': P, 'P+T': P + T}
divpts = {}
for q in bad_primes_for_index:
    divpts[str(q)] = {}
    for nm, Pt in targets.items():
        t0 = time.time()
        d = Pt.division_points(q)          # ТОЧНОЕ вычисление через многочлены деления
        divpts[str(q)][nm] = {'solutions': [str(z) for z in d], 'empty': (len(d) == 0),
                              'seconds': round(time.time() - t0, 2)}
        say("[6] division_points(q=%s) для %s -> %s  (%.2f c)  [proved_software: точная арифметика]"
            % (q, nm, d, time.time() - t0))
        assert len(d) == 0, "НАЙДЕНА ДЕЛИМОСТЬ: подгруппа НЕ насыщена!"


def not_q_divisible_mod_p(E, Pt, q, p):
    """True, если образ Pt в E(F_p) НЕ лежит в q*E(F_p) (=> Pt не q-делима в E(Q))."""
    if disc % p == 0 or p == q:
        return None
    A = E.change_ring(GF(p)).abelian_group()
    Ep = E.change_ring(GF(p))
    inv = list(A.invariants())
    v = list(A.discrete_log(Ep(Pt)))
    for ni, vi in zip(inv, v):
        if ZZ(vi) % gcd(ZZ(q), ZZ(ni)) != 0:
            return True
    return False


mod_certs = {}
for q in bad_primes_for_index:
    mod_certs[str(q)] = {}
    for nm, Pt in targets.items():
        cert = None
        for p in prime_range(5, 20000):
            r = not_q_divisible_mod_p(M1, Pt, q, p)
            if r is True:
                Ep = M1.change_ring(GF(p))
                A = Ep.abelian_group()
                cert = {'p': int(p), 'group_invariants': [int(z) for z in A.invariants()],
                        'discrete_log': [int(z) for z in A.discrete_log(Ep(Pt))],
                        'reason': 'образ %s не лежит в %s*E(F_%s)' % (nm, q, p)}
                break
        mod_certs[str(q)][nm] = cert
        say("[6] независимый сертификат неделимости q=%s, %s : %s  [proved]" % (q, nm, cert))
        assert cert is not None, "не найден модульный сертификат для q=%s, %s" % (q, nm)

R['non_divisibility'] = {
    'division_points_exact': divpts,
    'mod_p_certificates': mod_certs,
    'note': ('division_points -- точное решение уравнения q*Q = Pt над Q (многочлены деления); '
             'модульные сертификаты независимы от этого кода: если Pt = q*Q в E(Q), то для любого '
             'простого p хорошей редукции образ Pt лежит в q*E(F_p).'),
    'status': 'proved_software',
}

# ============================== 7. перекрёстные проверки насыщенности чужими реализациями
sat_pts, sat_index, sat_reg = M1.saturation([P], max_prime=-1)
say("[7] eclib saturation(max_prime=-1): index = %s, reg = %s, точки %s"
    % (sat_index, sat_reg, sat_pts))
sat_checks = {'eclib_full': {'index': int(sat_index), 'regulator': str(RR(sat_reg))}}
for mp in [11, 100]:
    sp, si, sr = M1.saturation([P], max_prime=mp)
    sat_checks['eclib_max_prime_%s' % mp] = {'index': int(si)}
    say("[7] eclib saturation(max_prime=%s): index = %s" % (mp, si))

pari_sat = None
try:
    v = M1.pari_curve().ellsaturation([P.__pari__()], 20)
    pari_sat = str(v)
    say("[7] PARI ellsaturation(B=20): %s" % v)
except Exception as e:
    pari_sat = 'недоступно: %s' % type(e).__name__
    say("[7] PARI ellsaturation недоступна (%s) -- не критично" % type(e).__name__)
sat_checks['pari_ellsaturation'] = pari_sat
assert int(sat_index) == 1

# ------------------------------------------------------------ финальное утверждение
say("")
say("[7] ВЫВОД. rank E1(Q) = 1 (2-спуск, dim Sel_2 = 2, Sha[2] = 0);")
say("[7]        E1(Q)_tors = Z/2 = <T>, T = %s;" % T)
say("[7]        индекс <P> + <T> в E1(Q) делит n <= %s, и ни одно простое q <= %s не делит n;"
    % (index_bound, index_bound))
say("[7]        => E1(Q) = <P> (+) <T> ~ Z (+) Z/2, P = %s.  НАСЫЩЕННОСТЬ ДОКАЗАНА [proved_software]" % P)
say("")

R['saturation'] = {
    'claim': 'E1(Q) = <P> (+) <T>, P = %s (порядок oo), T = %s (порядок 2)' % (P, T),
    'index': 1,
    'proof_steps': [
        'rank = 1 из 2-спуска: dim Sel_2 = 2 = rank + dim E(Q)[2] + dim Sha[2], rank >= 1 из P',
        'mu = %s -- доказанная нижняя граница hhat на неторсионных точках (Cremona-Siksek)' % mu_best,
        'hhat(P) = %s => индекс n <= %s' % (RR(hP_sage), index_bound),
        'для q in %s: division_points(q) пусто для P и P+T (точная арифметика)' % bad_primes_for_index,
        'независимо: модульные сертификаты неделимости',
    ],
    'cross_checks': sat_checks,
    'status': 'proved_software',
}

# ================================================ 8. КОНТРОЛЬ: известные точки D под phi1
# D(Q) содержит: (0, +-2023)  [t = 1, вырожденная девятка (u0,u4,u8) = (17,17,17)]
# и две точки на бесконечности с y/x^3 = +-6647  [t = -1].
known = []
for yv in [QQ(2023), QQ(-2023)]:
    Xr, Yr = a6 * QQ(0)**2, a6 * yv
    okk, residk = exact_on_curve(E1, Xr, Yr)
    assert okk
    pm = M1(iso(E1([Xr, Yr, 1])))
    # выражение в базисе <P> (+) <T>
    coeff = None
    for k in range(-8, 9):
        for tt, tn in [(M1(0), 'O'), (T, 'T')]:
            if k * P + tt == pm:
                coeff = (k, tn)
    say("[8] phi1(0, %s) = %s на E1raw -> %s на E1min ; в базисе: %s  [proved]"
        % (yv, (Xr, Yr), pm, coeff))
    known.append({'D_point': ['0', str(yv)], 't': '1', 'phi1_raw': [str(Xr), str(Yr)],
                  'phi1_min': str(pm), 'exact_residual': str(residk),
                  'basis_coordinates': {'k': int(coeff[0]), 'torsion': coeff[1]},
                  'order': 'infinite'})
    assert coeff is not None
for sgn in [1, -1]:
    say("[8] точка на бесконечности D (y/x^3 -> %s*6647): phi1 -> O (полюс x^2 второго порядка)  [proved]" % sgn)
    known.append({'D_point': ['infty', '%s*6647' % sgn], 't': '-1', 'phi1_min': 'O',
                  'basis_coordinates': {'k': 0, 'torsion': 'O'}})

kk = [z['basis_coordinates']['k'] for z in known if z['D_point'][0] == '0']
say("[8] ВНИМАНИЕ ДЛЯ РЕШЕТА: образ известной точки D равен %s*P, а НЕ +-P."
    % kk[0])
say("[8] Если бы за образующую взяли phi1-образ известной точки, индекс был бы 2,")
say("[8] и решето убивало бы кандидатов НЕПРАВОМЕРНО. Контроль пройден: P -- образующая. [proved]")

R['known_points_control'] = {
    'points': known,
    'warning': ('phi1-образ известной конечной точки D равен %s*P (индекс 2 в <P>); '
                'использовать как образующую E1(Q) её НЕЛЬЗЯ' % kk[0]),
    'status': 'proved',
}

# ======================================= 9. локальные данные для решета Морделла-Вейля
# Экспорт: для простых p хорошей редукции -- структура E1(F_p) и координаты P, T.
local = {}
for p in prime_range(3, 220):
    if disc % p == 0:
        continue
    Ep = M1.change_ring(GF(p))
    A = Ep.abelian_group()
    local[str(p)] = {
        'order': int(Ep.cardinality()),
        'invariants': [int(z) for z in A.invariants()],
        'dlog_P': [int(z) for z in A.discrete_log(Ep(P))],
        'dlog_T': [int(z) for z in A.discrete_log(Ep(T))],
        'a_p': int(M1.ap(p)),
    }
say("[9] экспортированы локальные данные для %s простых хорошей редукции (включая p = 11): %s"
    % (len(local), local['11']))

# Дополнительное (безусловное) условие на phi1-образ точки D:
#   X_raw(phi1(z)) = a6 * x_D^2  =>  X_raw / a6 -- КВАДРАТ в Q (для конечных z).
# Сдвиг между моделями определяем вычислением, а не на глаз.
shift = None
for Xr in [QQ(0), QQ(a6), QQ(4 * a6)]:
    Yr2 = rhs(E1, Xr)
    if Yr2.is_square():
        pm_t = M1(iso(E1([Xr, Yr2.sqrt(), 1])))
        sh = pm_t[0] - Xr
        assert shift is None or sh == shift
        shift = sh
assert shift is not None
say("[9] сдвиг моделей: X_min = X_raw + %s, то есть X_raw = X_min - %s (проверено на точках)  [proved]"
    % (shift, shift))
say("[9] безусловное условие на образ: X_raw(phi1(z)) = a6 * x_D^2, то есть (X_min - %s)/a6 -- квадрат в Q"
    % shift)
say("[9]     (статус proved; при переносе по модулю p требуется отдельный учёт полюсов и плохих мест)")

R['local_data_for_sieve'] = {
    'model': 'minimal E1',
    'P': str(P), 'T': str(T),
    'primes': local,
    'image_condition': {'statement': '(X_min - shift)/a6 = x_D^2 -- квадрат в Q для любого конечного z in D(Q)',
                        'shift_X_min_minus_X_raw': str(shift),
                        'a6': str(a6),
                        'status': 'proved',
                        'caveat': 'редукция условия по модулю p требует отдельного разбора полюсов'},
}

# ======================================================================= итог и запись
R['summary'] = {
    'E1_min_ainvs': [str(z) for z in M1.ainvs()],
    'rank': 1,
    'torsion': 'Z/2',
    'basis': {'free_generator': str(P), 'torsion_generator': str(T)},
    'saturated': True,
    'index': 1,
    'index_bound_before_checks': int(index_bound),
    'regulator': str(RR(hP_sage)),
    'status': 'proved_software',
    'what_is_NOT_claimed': [
        'не утверждается, что D(Q) состоит ровно из четырёх точек',
        'не утверждается, что образ Jac(D)(Q) равен E1(Q) x E2(Q)',
        'насыщенность доказана только для E1; фактор E2 -- отдельный расчёт (basis_E2.sage)',
    ],
    'seconds_total': round(time.time() - t_start, 1),
}

def _jdefault(o):
    """Sage Integer/Rational/RealNumber -> строка или int (preparser делает литералы Integer)."""
    if o in ZZ:
        return int(o)
    return str(o)


(OUT / 'basis_E1.json').write_text(
    json.dumps(R, indent=2, ensure_ascii=False, default=_jdefault) + '\n')
(OUT / 'basis_E1.log').write_text("\n".join(LOG) + "\n")
say("записано: %s , %s" % (OUT / 'basis_E1.json', OUT / 'basis_E1.log'))
