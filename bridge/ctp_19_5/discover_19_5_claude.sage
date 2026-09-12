# -*- coding: utf-8 -*-
# ОТКРЫТИЕ сертификата Касселса-Тейта для G1 (m,n) = (19,5).   Claude, 2026-09-12.
#
# Ранг кривой E НИГДЕ не вычисляется и НИГДЕ не используется.
# Из PARI берётся ТОЛЬКО ell2cover (список всюду локально разрешимых 2-накрытий = квартик).
# Локальная разрешимость всех используемых квартик затем перепроверяется явными свидетелями,
# так что даже вывод ell2cover не входит в цепочку доверия сертификата.
#
# Выход: certificate_19_5_claude.json  -> проверяется bridge/ctp_cert_19_5_фишер.py на чистом Python.

from sage.all import *
from cysignals.alarm import alarm, cancel_alarm
from cysignals.signals import AlarmInterrupt
from pathlib import Path
import json, random, time

OUT = Path('/home/kep/magicKube/bridge/ctp_19_5')
load('/home/kep/magicKube/descent/ctp_quartic.sage')
random.seed(int(19052026)); set_random_seed(int(19052026))

S = lambda v: str(v)
SL = lambda v: [str(c) for c in v]

# ---------------------------------------------------------------- 1. семейство G1 и класс delta
m, n = 19, 5
assert gcd(m, n) == 1 and (m * m + n * n) % 2 == 0
s = QQ(m ** 2 + n ** 2) / 2                       # 193
b = s * m ** 2 * n ** 2                           # 1741825
assert s == 193 and b == 1741825
assert not QQ(s).is_square()                      # исключение бесконечно удалённых точек C

# ТОЖДЕСТВА, задающие класс delta: проверяем как многочлены от t (точно, а не численно)
Rt = PolynomialRing(QQ, 't'); t = Rt.gen()
F0 = m ** 2 + n ** 2 * t ** 2
F4 = s * (1 + t ** 2)
F8 = n ** 2 + m ** 2 * t ** 2
X = b * t ** 2
e = [-b, -s * m ** 4, -s * n ** 4]
assert X - e[0] == (m * n) ** 2 * F4, 'X-e1 != (mn)^2 F4'
assert X - e[1] == s * m ** 2 * F0, 'X-e2 != s m^2 F0'
assert X - e[2] == s * n ** 2 * F8, 'X-e3 != s n^2 F8'
# на C:  F4 = u4^2, F0 = u0^2, F8 = u8^2  =>  классы (X-e1, X-e2, X-e3) = (1, s, s)
TRIP = [QQ(1), QQ(s), QQ(s)]

# ---------------------------------------------------------------- 2. кривая и минимальная модель
Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
f = prod(x - ee for ee in e)
E = EllipticCurve([0, f[2], 0, f[1], f[0]])
assert all(E.defining_polynomial()(ee, 0, 1) == 0 for ee in e)
M = E.minimal_model()
iso = E.isomorphism_to(M)
mr = [iso(E([ee, 0]))[0] for ee in e]
u = iso.u

# точная связь: X - e_i = lam * (x_M - mr_i) с ОДНИМ квадратным lam  =>  квадратный класс сохраняется
LAM = (e[0] - e[1]) / (mr[0] - mr[1])
assert QQ(LAM).is_square(), 'масштаб между моделями не квадрат — класс мог сдвинуться'
SH = e[0] - LAM * mr[0]
for i in range(3):
    assert e[i] == LAM * mr[i] + SH, 'сдвиг моделей не согласован'

I, J = IJ_of_curve(M)
F = FisherCTP(QQ, I, J, verbose=True)
phis = [phi_of_root(M, ee) for ee in mr]
assert len(F.E.comps) == 3 and all(c['deg'] == 1 for c in F.E.comps), \
    'L не расщепилась на три компоненты — полное 2-кручение над Q ожидалось'
order = _comp_order(F, phis)
delta = F.E.from_comps([TRIP[i] for i in order])

rec = {'family': 'G1', 'm': int(m), 'n': int(n), 's': S(s), 'b': S(b),
       'orig_roots': SL(e), 'lam': S(LAM), 'shift': S(SH),
       'M_ainvs': SL(M.ainvs()), 'minimal_roots': SL(mr), 'phi_roots': SL(phis),
       'I': S(I), 'J': S(J), 'discriminant': S(F.disc), 'order': [int(o) for o in order],
       'delta_trip': SL(TRIP), 'delta': S(delta), 'pairs': []}
print('SETUP', rec['M_ainvs'], rec['minimal_roots'], 'order', order, flush=True)

# ---------------------------------------------------------------- 3. квартика класса delta
t0 = time.time()
g1 = F.quartic_from_delta(delta, tag='target')
rec['target_quartic'] = SL(g1); rec['target_seconds'] = time.time() - t0
print('TARGET', g1, flush=True)

# контроль связи с исходной C: z(g1) лежит в классе (1, 193, 193) ПОКОМПОНЕНТНО в порядке (e1,e2,e3)
z1 = [QQ(F.z_inv(g1).lift()(ph)) for ph in phis]
assert all((z1[i] / TRIP[i]).is_square() for i in range(3)), 'z(g1) не в классе (1,193,193)'
print('z(g1) класс подтверждён', z1, flush=True)

# ---------------------------------------------------------------- 4. партнёры из ell2cover
cover = pari(M).ell2cover()
rec['n_cover'] = len(cover)
found = None
for i, c in enumerate(cover):
    t1 = time.time()
    try:
        alarm(1800)
        q = Rx(c[0]); g2 = [QQ(q[4 - j]) for j in range(5)]
        scale = QQ(J * quartic_I(g2) / (quartic_J(g2) * I))
        assert scale.is_square(), 'scale %s не квадрат' % scale
        g2 = [scale * a0 for a0 in g2]
        F.check_quartic(g2)
        g2 = F._make_z_unit(g2)
        g3 = F.quartic_from_delta(F.z_inv(g1) * F.z_inv(g2), tag='sum%d' % i)
        val, npl = F.pair(g1, g2, g3, reps=3)
        cancel_alarm()
        row = {'i': i, 'scale': S(scale), 'g2': SL(g2), 'g3': SL(g3),
               'pair': int(val), 'places': int(npl), 'seconds': time.time() - t1}
        rec['pairs'].append(row); print('PAIR', row['i'], row['pair'], flush=True)
        if val and found is None:
            found = row
    except (Exception, AlarmInterrupt) as ex:
        cancel_alarm()
        rec['pairs'].append({'i': i, 'error': repr(ex), 'seconds': time.time() - t1})
        print('ERROR', i, repr(ex), flush=True)
    (OUT / 'ctp_probe_19_5_claude.json').write_text(json.dumps(rec, indent=2) + '\n')

assert found is not None, 'ненулевого спаривания не найдено — сертификат НЕ построен'
rec['chosen'] = found['i']

# ---------------------------------------------------------------- 5. экспорт сертификата
gs = [list(map(QQ, g)) for g in [rec['target_quartic'], found['g2'], found['g3']]]
for g in gs:
    F.check_quartic(g)
zs = [[QQ(F.z_inv(g).lift()(ph)) for ph in phis] for g in gs]
assert all((zs[0][i] / TRIP[i]).is_square() for i in range(3))

rawgam, mm = F.gamma1(*gs)
gam = list(rawgam)
den = lcm([QQ(c).denominator() for c in gam])
gam = [c * den for c in gam]
num = gcd([ZZ(QQ(c).numerator()) for c in gam if c != 0])
gam = [c / num for c in gam]
gscale = QQ(den) / num

a = gs[1][0]                    # ДОСЛОВНО g2(1,0) из Теоремы 3.1 (никакого small_rep)
assert a != 0
pls = F.places_for(gs[0], gam, a)

rec['quartics'] = [SL(g) for g in gs]
rec['z_values'] = [SL(z) for z in zs]
rec['m_values'] = SL([mm.lift()(ph) for ph in phis])
rec['gamma_raw'] = SL(rawgam); rec['gamma'] = SL(gam); rec['gamma_scale'] = S(gscale)
rec['a'] = S(a)
rec['pair_runs'] = []; rec['local_solubility'] = []

for rep in range(3):
    F._lpcache = {}
    loc = []
    for pl in pls:
        xx, zz = F.local_point(gs[0], gam, pl)
        qv = quartic_eval(gs[0], xx, zz)
        gv = gam[0] * xx ** 2 + gam[1] * xx * zz + gam[2] * zz ** 2
        sg = hilb(QQ, a, gv, pl)
        loc.append({'place': pl.name, 'x': S(xx), 'z': S(zz), 'g': S(qv),
                    'gamma': S(gv), 'hilbert': int(sg)})
    assert prod(v['hilbert'] for v in loc) == -1, [(v['place'], v['hilbert']) for v in loc]
    rec['pair_runs'].append(loc)
    print('PAIR_RUN', rep, [(v['place'], v['hilbert']) for v in loc if v['hilbert'] == -1], flush=True)

# симметрия и Замечание 3.2(v)
rev = F.pair(gs[1], gs[0], gs[2], reps=3); rec['reverse_pair'] = int(rev[0])
alt13 = F.pair(gs[0], gs[2], gs[1], reps=3); rec['pair_g1_g3'] = int(alt13[0])
alt31 = F.pair(gs[2], gs[0], gs[1], reps=3); rec['pair_g3_g1'] = int(alt31[0])
diag = F.pair(gs[0], gs[0], F.quartic_from_delta(F.z_inv(gs[0]) ** 2, tag='diag'), reps=2)
rec['pair_g1_g1'] = int(diag[0])
print('REVERSE', rev[0], 'g1g3', alt13[0], 'g3g1', alt31[0], 'diag', diag[0], flush=True)
assert rec['reverse_pair'] == 1, 'СИММЕТРИЯ НАРУШЕНА'

# ELS всех трёх квартик явными свидетелями
for i, g in enumerate(gs):
    ps = set(ZZ(F.disc.numerator()).prime_divisors()) | set(ZZ(F.disc.denominator()).prime_divisors()) | {ZZ(2), ZZ(3)}
    for a0 in g:
        if a0 != 0:
            ps |= set(QQ(a0).denominator().prime_divisors())
    local = []
    for pl in [LocSq(QQ, p) for p in sorted(ps)] + [RealPlace(QQ)]:
        xx, zz = F.local_point(g, [QQ(1), QQ(0), QQ(1)], pl)
        local.append({'place': pl.name, 'x': S(xx), 'z': S(zz), 'value': S(quartic_eval(g, xx, zz))})
    rec['local_solubility'].append(local)
    print('ELS', i, [v['place'] for v in local], flush=True)

(OUT / 'certificate_19_5_claude.json').write_text(json.dumps(rec, indent=2) + '\n')
print('CERTIFICATE EXPORTED', flush=True)
