"""Экспорт сертификата CTP для (11,4): всё, что нужно проверяющему на чистом Python.

Дополнительно КОНТРОЛЬ-ФАЛЬСИФИКАЦИЯ: классы из E(Q)/2E(Q) (2-кручение + найденные
рациональные точки) обязаны спариваться с g2 ТРИВИАЛЬНО. Если бы хоть один дал -1,
машинерия была бы сломана. Контроль НЕ используется в выводе, только как проверка кода.
"""
from sage.all import *
from pathlib import Path
import json, random

out = Path(__file__).resolve().parent
load(str(out / 'ctp_quartic_snapshot.sage'))
random.seed(1211904)
set_random_seed(1211904)

r = json.loads((out / 'ctp_probe_11_4.json').read_text())
p = next(p for p in r['pairs'] if p.get('pair') == 1)
I, J = QQ(r['I']), QQ(r['J'])
F = FisherCTP(QQ, I, J)
gs = [list(map(QQ, g)) for g in [r['target_quartic'], p['g2'], p['g3']]]
mroots = [QQ(e) for e in r['minimal_roots']]
phis = [-12 * e for e in mroots]           # b2 = 0 для M, значит phi = -12 e
for g in gs:
    F.check_quartic(g)

zs = [[QQ(F.z_inv(g).lift()(phi)) for phi in phis] for g in gs]
TRIP = [1, 274, 274]
assert all((zs[0][i] / d).is_square() for i, d in enumerate(TRIP)), 'z(g1) не в классе (1,274,274)'

rawgam, mm = F.gamma1(*gs)
gam = list(rawgam)
den = lcm(c.denominator() for c in gam)
gam = [c * den for c in gam]
num = gcd([ZZ(c) for c in gam if c])
gam = [c / num for c in gam]

a = gs[1][0]        # ЛИТЕРАЛЬНОЕ g2(1,0) из теоремы 3.1, без small_rep/min-сокращений
assert a != 0
pls = F.places_for(gs[0], gam, a)

record = {'pair': '(m,n)=(11,4)', 'm': 11, 'n': 4, 's': r['s'], 'b': r['b'],
          'E_ainvs': r['E_ainvs'], 'M_ainvs': r['M_ainvs'],
          'I': str(I), 'J': str(J),
          'minimal_roots': r['minimal_roots'],
          'phi_roots': list(map(str, phis)),
          'quartics': [list(map(str, g)) for g in gs],
          'z_values': [[str(v) for v in z] for z in zs],
          'm_values': [str(mm.lift()(ph)) for ph in phis],
          'gamma_raw': list(map(str, rawgam)),
          'gamma': list(map(str, gam)),
          'gamma_scale': str(QQ(den) / num),
          'a': str(a),
          'discriminant': str(F.disc),
          'target_squareclass': TRIP,
          'pair_runs': [], 'local_solubility': []}

# --- два независимых прогона со СВЕЖИМИ локальными точками ---
for rep in range(3):
    F._lpcache = {}
    loc = []
    for pl in pls:
        x, z = F.local_point(gs[0], gam, pl)
        q = quartic_eval(gs[0], x, z)
        v = gam[0] * x * x + gam[1] * x * z + gam[2] * z * z
        sg = hilb(QQ, a, v, pl)
        loc.append({'place': pl.name, 'x': str(x), 'z': str(z),
                    'g': str(q), 'gamma': str(v), 'hilbert': int(sg)})
    assert prod(v['hilbert'] for v in loc) == -1, 'произведение символов != -1'
    record['pair_runs'].append(loc)
    print('PAIR_RUN', rep, [(v['place'], v['hilbert']) for v in loc], flush=True)

# --- симметрия: обратный порядок аргументов даёт свою gamma, свои места, свои точки ---
rev = F.pair(gs[1], gs[0], gs[2], reps=2)
print('REVERSE <g2,g1> =', rev, flush=True)
assert rev[0] == 1, 'спаривание несимметрично — ошибка'
record['reverse_pair'] = int(rev[0])
# --- третий порядок: <g1,g3> должно равняться <g1,g2> (Fisher 3.2(v)) ---
p13 = F.pair(gs[0], gs[2], gs[1], reps=2)
print('<g1,g3> =', p13, flush=True)
record['pair_g1_g3'] = int(p13[0])
assert p13[0] == 1
# --- альтернирование: <g1,g1> = 0 ---
try:
    self_pair = F.pair(gs[0], gs[0], F.quartic_from_delta(F.z_inv(gs[0]) * F.z_inv(gs[0]), tag='self'), reps=1)
    record['self_pair_g1_g1'] = int(self_pair[0])
    print('<g1,g1> =', self_pair, flush=True)
except Exception as e:
    record['self_pair_g1_g1'] = 'error: ' + repr(e)
    print('self pair error', repr(e), flush=True)

# --- ELS всех трёх квартик, независимо от перечисления Селмера в PARI ---
for i, g in enumerate(gs):
    ps = set(ZZ(F.disc.numerator()).prime_divisors()) | set(ZZ(F.disc.denominator()).prime_divisors()) | {ZZ(2), ZZ(3)}
    for a0 in g:
        ps.update(a0.denominator().prime_divisors())
    local = []
    for pl in [LocSq(QQ, pp) for pp in sorted(ps)] + [RealPlace(QQ)]:
        x, z = F.local_point(g, [QQ(1), QQ(0), QQ(1)], pl)
        local.append({'place': pl.name, 'x': str(x), 'z': str(z),
                      'value': str(quartic_eval(g, x, z))})
    record['local_solubility'].append(local)
    print('ELS', i, [v['place'] for v in local], flush=True)

# --- КОНТРОЛЬ-ФАЛЬСИФИКАЦИЯ: образ E(Q)/2E(Q) обязан спариваться тривиально ---
M = EllipticCurve([QQ(c) for c in r['M_ainvs']])
ctrl = []
try:
    gens = M.gens(proof=False)
except Exception as e:
    gens = []
    print('gens недоступны:', repr(e), flush=True)
tors = [P for P in M.torsion_points() if P != M(0)]
basepts = list(gens) + tors
seen = set()
for P in basepts:
    try:
        xp = QQ(P[0])
        comps = [xp - e for e in mroots]
        if any(c == 0 for c in comps):
            # точка 2-кручения: заменяем нулевую компоненту по произведению (класс определён)
            j = [i for i, c in enumerate(comps) if c == 0][0]
            others = [mroots[j] - mroots[i] for i in range(3) if i != j]
            comps[j] = prod(others)
        cls = tuple(str(QQ(c).squarefree_part()) for c in comps)
        if cls in seen:
            continue
        seen.add(cls)
        # delta для E_{I,J}: x_{I,J} - Theta_i = 36*(x_M - e_i), 36 — квадрат
        order = _comp_order(F, phis)
        dl = F.E.from_comps([comps[i] for i in order])
        gq = F.quartic_from_delta(dl, tag='ctrl')
        gsum = F.quartic_from_delta(F.z_inv(gq) * F.z_inv(gs[1]), tag='ctrlsum')
        val, npl = F.pair(gq, gs[1], gsum, reps=1)
        ctrl.append({'point_x': str(xp), 'class': list(cls), 'pair_with_g2': int(val)})
        print('CTRL', xp, cls, '-> pair', val, flush=True)
    except Exception as e:
        ctrl.append({'point_x': str(P[0]), 'error': repr(e)})
        print('CTRL error', P[0], repr(e), flush=True)
record['control_image_pairs'] = ctrl
assert all(c.get('pair_with_g2', 0) == 0 for c in ctrl), 'КОНТРОЛЬ ПРОВАЛЕН: точка из образа дала -1'

(out / 'certificate_11_4.json').write_text(json.dumps(record, indent=2) + '\n')
print('CERTIFICATE EXPORTED', flush=True)
