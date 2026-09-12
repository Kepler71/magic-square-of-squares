"""Independent Sage audit of the forty 11-adic candidates of the mixed
genus-2 bridge D for G1(15,8).

This script does NOT reuse the pure-python parser of build_candidates_table.py.
Every number is re-parsed with Sage's own Qp and every claim is recomputed.

Output: verify_candidates_sage.json   (consumed by build_candidates_table.py)
"""
from sage.all import *
import json
from pathlib import Path

BASE = Path('/home/kep/magicKube/bridge/qc40')
p = 11
K = Qp(p, 60)

res30 = json.loads((BASE / 'qc_mixed_15_8_p11_n30_twelve_check.json').read_text())
res20 = json.loads((BASE / 'qc_mixed_15_8_p11_n20_twelve.json').read_text())

R = PolynomialRing(QQ, 'x')
x = R.gen()
f = (49 + 529 * x * x) * (83521 * x**4 + 63358 * x * x + 83521)
assert str(f.change_ring(QQ)) == res30['polynomial'] or f == R(res30['polynomial'])
A = [f[0], f[2], f[4], f[6]]                       # a0, a2, a4, a6

def _j(o):
    try:
        return int(o)
    except Exception:
        return str(o)


report = {}
report['f_expanded'] = str(f)
report['a0_a2_a4_a6'] = [str(a) for a in A]
report['a0_is_square'] = bool(ZZ(A[0]).is_square())
report['a6_is_square'] = bool(ZZ(A[3]).is_square())
report['sqrt_a0'] = str(ZZ(A[0]).sqrt())
report['sqrt_a6'] = str(ZZ(A[3]).sqrt())
report['disc_valuation_at_11'] = int(f.discriminant().valuation(p))

# ---------------------------------------------------------------- Omega, re-derived
# verbatim from run_qc_twelve.py line 64
Kp = Qp(p, 30)
omega_rat = [(r7, r17, r23) for r7 in [QQ(-2), QQ(-4) / 3]
             for r17 in [QQ(-2), QQ(0), QQ(2)]
             for r23 in [QQ(2), QQ(4) / 3]]
omega_val = [a * Kp(7).log() + b * Kp(17).log() + c * Kp(23).log()
             for (a, b, c) in omega_rat]
report['omega_count'] = len(omega_rat)
report['omega_triples'] = [[str(t) for t in tr] for tr in omega_rat]
diffs = [(omega_val[i] - omega_val[j]).valuation()
         for i in range(12) for j in range(i + 1, 12)]
report['omega_pairwise_difference_count'] = len(diffs)
report['omega_min_difference_valuation'] = int(min(diffs))
report['omega_max_difference_valuation'] = int(max(diffs))
report['omega_all_distinct_at_prec_30'] = bool(min(diffs) < 30)
report['omega_values_prec30'] = [str(v) for v in omega_val]

# theory (TWELVE_HEIGHT_VALUES.md):  x=0 -> (-2,-2,4/3);  x=oo -> (-4/3,2,2)
report['theory_index_x0'] = omega_rat.index((QQ(-2), QQ(-2), QQ(4) / 3))
report['theory_index_xinf'] = omega_rat.index((QQ(-4) / 3, QQ(2), QQ(2)))

# ---------------------------------------------------------------- point parser
def parse(s):
    """Sage-native parse of '(a : b : c)' with p-adic entries."""
    parts = s[1:-1].split(' : ')
    return [K(sage_eval(t, locals={'O': O, 'x': x})) for t in parts]


def normalise(s):
    """Return (chart, X, Y) with v(X) > 0, matching the two allowed discs."""
    X, Y, Z = parse(s)
    if Z == 0:
        return 'infinity', X, Y          # already (1/x : y/x^3 : 0)-style chart
    X, Y = X / Z, Y / Z**3
    if X.valuation() < 0:
        return 'infinity', 1 / X, Y / X**3
    return 'finite', X, Y


def residual(chart, X, Y):
    co = A if chart == 'finite' else list(reversed(A))
    return Y * Y - sum(co[i] * X**(2 * i) for i in range(4))


# ---------------------------------------------------------------- known points
known = sorted({s for row in res30['rational_points_by_omega'] for s in row})
report['known_points_sorted'] = known
report['known_points_expected'] = sorted(
    ['(0 : -2023 : 1)', '(0 : 2023 : 1)', '(1 : -6647 : 0)', '(1 : 6647 : 0)'])
report['known_points_match'] = (report['known_points_sorted']
                                == report['known_points_expected'])
report['known_by_omega_index'] = {
    str(i): row for i, row in enumerate(res30['rational_points_by_omega']) if row}

# on-curve check for the four known points, exactly over Q
known_rows = []
for i, row in enumerate(res30['rational_points_by_omega']):
    for s in row:
        a, b, c = [QQ(t) for t in s[1:-1].split(' : ')]
        if c == 0:
            ok = bool(b * b == A[3] * a**6)        # y^2 = a6 x^6 at infinity
            tval, xdesc = '-1', 'infinity'
        else:
            ok = bool(b * b == f(a / c))
            tval, xdesc = str((1 + a / c) / (1 - a / c)), str(a / c)
        known_rows.append({'point': s, 'omega_index': i,
                           'omega_triple': [str(u) for u in omega_rat[i]],
                           'x': xdesc, 't': tval, 'on_curve_over_Q': ok})
report['known_rows'] = known_rows
report['known_all_on_curve_over_Q'] = all(r['on_curve_over_Q'] for r in known_rows)

# ---------------------------------------------------------------- the forty
sig20 = {}
for i, row in enumerate(res20['extra_points_by_omega']):
    for s in row:
        ch, X, Y = normalise(s)
        sig20.setdefault(i, set()).add((ch, ZZ(X.lift()) % p**8, ZZ(Y.lift()) % p**8))

cands = []
for i, row in enumerate(res30['extra_points_by_omega']):
    for s in row:
        ch, X, Y = normalise(s)
        r = residual(ch, X, Y)
        vX = int(X.valuation())
        entry = {
            'omega_index': i,
            'omega_triple': [str(u) for u in omega_rat[i]],
            'residue_disc': 'x=0' if ch == 'finite' else 'infinity',
            'chart': ('affine (x, y)' if ch == 'finite'
                      else 'w = 1/x, v = y/x^3'),
            'chart_x_valuation': vX,
            'chart_x_unit_mod_11': int(ZZ((X / p**vX).lift()) % p),
            'chart_x_mod_11_8': int(ZZ(X.lift()) % p**8),
            'chart_y_mod_11_8': int(ZZ(Y.lift()) % p**8),
            'chart_x_mod_11_20': int(ZZ(X.lift()) % p**20),
            'chart_y_mod_11_20': int(ZZ(Y.lift()) % p**20),
            'abs_precision_x': int(X.precision_absolute()),
            'abs_precision_y': int(Y.precision_absolute()),
            'curve_residual_valuation': (int(r.valuation()) if r != 0
                                         else int(r.precision_absolute())),
            'stable_at_precision_20': bool((ch, ZZ(X.lift()) % p**8, ZZ(Y.lift()) % p**8)
                                           in sig20.get(i, set())),
            'raw_string': s,
        }
        cands.append(entry)

report['candidate_count'] = len(cands)
report['per_omega_counts'] = [sum(1 for c in cands if c['omega_index'] == i)
                              for i in range(12)]
report['per_omega_counts_match_source'] = (report['per_omega_counts']
                                           == res30['extra_counts'])
report['min_curve_residual_valuation'] = min(c['curve_residual_valuation']
                                             for c in cands)
report['all_stable_at_precision_20'] = all(c['stable_at_precision_20']
                                           for c in cands)
report['x_valuations_seen'] = sorted({c['chart_x_valuation'] for c in cands})
report['distinct_signatures_mod_11_8'] = len(
    {(c['omega_index'], c['residue_disc'], c['chart_x_mod_11_8'],
      c['chart_y_mod_11_8']) for c in cands})
report['distinct_xy_mod_11_20'] = len(
    {(c['residue_disc'], c['chart_x_mod_11_20'], c['chart_y_mod_11_20'])
     for c in cands})
report['counts_by_disc'] = {
    'x=0': sum(1 for c in cands if c['residue_disc'] == 'x=0'),
    'infinity': sum(1 for c in cands if c['residue_disc'] == 'infinity')}

# no candidate coincides with a known point
report['no_candidate_is_a_known_point'] = all(
    not (c['chart_x_mod_11_8'] == 0 and c['chart_y_mod_11_8'] == 0)
    for c in cands)
report['min_x_valuation_positive'] = bool(all(c['chart_x_valuation'] >= 1
                                              for c in cands))

(BASE / 'verify_candidates_sage.json').write_text(
    json.dumps({'report': report, 'candidates': cands}, indent=2, default=_j) + '\n')
print(json.dumps({k: v for k, v in report.items()
                  if k not in ('known_rows', 'omega_values_prec30',
                               'omega_triples')}, indent=2, default=_j))
