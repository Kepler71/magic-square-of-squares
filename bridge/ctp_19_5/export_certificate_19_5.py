"""Экспорт сертификата CTP для G1 (19,5) в certificate_19_5.json.
Использует Sage только на этапе ОТКРЫТИЯ; сам сертификат проверяется чистым Python
(ctp_cert_19_5_фишер.py).  Ранг кривой нигде не вычисляется."""
from sage.all import *
from pathlib import Path
import json, random

out = Path('/home/kep/magicKube/bridge/ctp_19_5')
load(str(out / 'ctp_quartic_snapshot.sage'))
random.seed(12090195); set_random_seed(12090195)

r = json.loads((out / 'ctp_probe_19_5.json').read_text())
p = next(q for q in r['pairs'] if q.get('pair') == 1)
I, J = QQ(r['I']), QQ(r['J'])
F = FisherCTP(QQ, I, J)
gs = [list(map(QQ, g)) for g in [r['target_quartic'], p['g2'], p['g3']]]
phis = [QQ(v) for v in r['phi_roots']]
mr = [QQ(v) for v in r['minimal_roots']]
assert phis == [-12 * e for e in mr]

for g in gs:
    F.check_quartic(g)
zs = [[QQ(F.z_inv(g).lift()(phi)) for phi in phis] for g in gs]
trip = [QQ(1), QQ(193), QQ(193)]
assert all((zs[0][i] / d).is_square() for i, d in enumerate(trip)), 'z(g1) не в классе (1,193,193)'

rawgam, mm = F.gamma1(*gs)
gam = list(rawgam)
den = lcm(c.denominator() for c in gam)
gam = [c * den for c in gam]
num = gcd([ZZ(c) for c in gam if c])
gam = [c / num for c in gam]

a = gs[1][0]                      # ДОСЛОВНО g2(1,0) из Теоремы 3.1, без small_rep
assert a != 0
pls = F.places_for(gs[0], gam, a)

# собственный дискриминант бинарной квартики (произведение квадратов разностей корней)
Rx = PolynomialRing(QQ, 'x'); xx = Rx.gen()
poly = sum(gs[0][j] * xx ** (4 - j) for j in range(5))
disc_form = poly.discriminant() / gs[0][0] ** 2 * gs[0][0] ** 2  # = disc бинарной формы при deg 4

record = {
    'family': 'G1', 'm': 19, 'n': 5, 's': '193', 'b': '1741825',
    'orig_roots': r['orig_roots'], 'M_ainvs': r['M_ainvs'],
    'minimal_roots': r['minimal_roots'], 'phi_roots': r['phi_roots'],
    'I': str(I), 'J': str(J), 'discriminant': str(F.disc),
    'delta_trip': ['1', '193', '193'],
    'quartics': [list(map(str, g)) for g in gs],
    'z_values': [[str(v) for v in z] for z in zs],
    'm_values': [str(mm.lift()(ph)) for ph in phis],
    'gamma_raw': list(map(str, rawgam)), 'gamma': list(map(str, gam)),
    'gamma_scale': str(QQ(den) / num), 'a': str(a),
    'pair_runs': [], 'local_solubility': [],
}

for rep in range(3):
    F._lpcache = {}
    loc = []
    for pl in pls:
        x, z = F.local_point(gs[0], gam, pl)
        q = quartic_eval(gs[0], x, z)
        v = gam[0] * x * x + gam[1] * x * z + gam[2] * z * z
        sg = hilb(QQ, a, v, pl)
        loc.append({'place': pl.name, 'x': str(x), 'z': str(z), 'g': str(q),
                    'gamma': str(v), 'hilbert': int(sg)})
    assert prod(v['hilbert'] for v in loc) == -1, [ (v['place'],v['hilbert']) for v in loc ]
    record['pair_runs'].append(loc)
    print('PAIR_RUN', rep, [(v['place'], v['hilbert']) for v in loc], flush=True)

# симметрия: обратный порядок аргументов -> та же величина
rev = F.pair(gs[1], gs[0], gs[2], reps=2)
assert rev[0] == 1, rev
record['reverse_pair'] = int(rev[0])
# и третья пара (Замечание 3.2(v)): <g1,g2> = <g1,g3>
alt = F.pair(gs[0], gs[2], gs[1], reps=2)
record['pair_g1_g3'] = int(alt[0])
print('REVERSE', rev, 'ALT <g1,g3>', alt, flush=True)

# ELS всех трёх квартик: явные свидетели во всех плохих местах
for i, g in enumerate(gs):
    ps = set(ZZ(F.disc.numerator()).prime_divisors()) | set(ZZ(F.disc.denominator()).prime_divisors()) | {ZZ(2), ZZ(3)}
    for a0 in g:
        if a0 != 0:
            ps.update(QQ(a0).denominator().prime_divisors())
    local = []
    for pl in [LocSq(QQ, pp) for pp in sorted(ps)] + [RealPlace(QQ)]:
        x, z = F.local_point(g, [QQ(1), QQ(0), QQ(1)], pl)
        local.append({'place': pl.name, 'x': str(x), 'z': str(z), 'value': str(quartic_eval(g, x, z))})
    record['local_solubility'].append(local)
    print('ELS', i, [v['place'] for v in local], flush=True)

(out / 'certificate_19_5.json').write_text(json.dumps(record, indent=2) + '\n')
print('CERTIFICATE EXPORTED', flush=True)
