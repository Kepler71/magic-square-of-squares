"""Structure of the forty candidates over E1 x E2.

phi1 = (a6 x^2, a6 y), phi2 = (a0/x^2, a0 y/x^3).
Both X-coordinates depend on x only through x^2, and X1*X2 = a0*a6.
So each orbit of <sigma: x->-x, iota: y->-y> has ONE pair (X1, X2) and four
sign choices (eps1, eps2).  This is checked numerically here, not assumed.
"""
from sage.all import *
import json
from pathlib import Path

BASE = Path('/home/kep/magicKube/bridge/qc40')
p = 11
K = Qp(p, 60)
M = p**18


def _j(o):
    try:
        return int(o)
    except Exception:
        return str(o)


tab = json.loads((BASE / 'candidates_table.json').read_text())
R = PolynomialRing(QQ, 'x')
xv = R.gen()
f = (49 + 529 * xv * xv) * (83521 * xv**4 + 63358 * xv * xv + 83521)
a0, a2, a4, a6 = f[0], f[2], f[4], f[6]


def parse(s):
    return [K(sage_eval(t, locals={'O': O, 'x': xv})) for t in s[1:-1].split(' : ')]


rows = []
for c in tab['candidates']:
    X, Y, Z = parse(c['raw_point_string_from_run'])
    X, Y = X / Z, Y / Z**3
    X1, Y1 = a6 * X**2, a6 * Y
    X2, Y2 = a0 / X**2, a0 * Y / X**3
    rows.append({
        'id': c['id'], 'omega_index': c['omega_index'],
        'symmetry_orbit': c['symmetry_orbit'],
        'residue_disc': c['residue_disc'],
        'X1_val': int(X1.valuation()), 'X2_val': int(X2.valuation()),
        'X1_key': '%d:%d' % (X1.valuation(), ZZ((X1 / p**X1.valuation()).lift()) % M),
        'X2_key': '%d:%d' % (X2.valuation(), ZZ((X2 / p**X2.valuation()).lift()) % M),
        'Y1_key': '%d:%d' % (Y1.valuation(), ZZ((Y1 / p**Y1.valuation()).lift()) % M),
        'Y2_key': '%d:%d' % (Y2.valuation(), ZZ((Y2 / p**Y2.valuation()).lift()) % M),
        'X1X2_minus_a0a6_val': int((X1 * X2 - a0 * a6).valuation())
        if X1 * X2 - a0 * a6 != 0 else 999,
    })

rep = {}
rep['distinct_X1'] = len({r['X1_key'] for r in rows})
rep['distinct_X2'] = len({r['X2_key'] for r in rows})
rep['distinct_(X1,X2)'] = len({(r['X1_key'], r['X2_key']) for r in rows})
rep['X1X2_equals_a0a6_min_valuation'] = min(r['X1X2_minus_a0a6_val'] for r in rows)
rep['a0a6'] = str(a0 * a6)
# each symmetry orbit: one (X1,X2), four points
byorb = {}
for r in rows:
    byorb.setdefault(r['symmetry_orbit'], []).append(r)
rep['orbit_count'] = len(byorb)
rep['orbit_sizes'] = sorted(len(v) for v in byorb.values())
rep['each_orbit_has_single_X1'] = all(len({r['X1_key'] for r in v}) == 1
                                      for v in byorb.values())
rep['each_orbit_has_single_X2'] = all(len({r['X2_key'] for r in v}) == 1
                                      for v in byorb.values())
rep['each_orbit_has_four_distinct_(Y1,Y2)'] = all(
    len({(r['Y1_key'], r['Y2_key']) for r in v}) == 4 for v in byorb.values())
rep['orbits_have_single_omega'] = all(len({r['omega_index'] for r in v}) == 1
                                      for v in byorb.values())
rep['orbits_have_single_disc'] = all(len({r['residue_disc'] for r in v}) == 1
                                     for v in byorb.values())
rep['omega_per_orbit'] = {str(k): int(v[0]['omega_index']) for k, v in sorted(byorb.items())}
rep['disc_per_orbit'] = {str(k): v[0]['residue_disc'] for k, v in sorted(byorb.items())}

(BASE / 'orbit_structure_sage.json').write_text(
    json.dumps({'report': rep, 'rows': rows}, indent=2, default=_j) + '\n')
print(json.dumps(rep, indent=2, default=_j))
