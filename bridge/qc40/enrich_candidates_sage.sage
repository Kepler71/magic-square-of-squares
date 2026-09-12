"""Enrich the forty candidates: symmetry structure, t = (1+x)/(1-x),
and the images under phi1 = (a6 x^2, a6 y), phi2 = (a0/x^2, a0 y/x^3).

Everything recomputed in Sage from the raw point strings.
Output: enrich_candidates_sage.json
"""
from sage.all import *
import json
from pathlib import Path

BASE = Path('/home/kep/magicKube/bridge/qc40')
p = 11
K = Qp(p, 60)


def _j(o):
    try:
        return int(o)
    except Exception:
        return str(o)


res30 = json.loads((BASE / 'qc_mixed_15_8_p11_n30_twelve_check.json').read_text())
R = PolynomialRing(QQ, 'x')
xv = R.gen()
f = (49 + 529 * xv * xv) * (83521 * xv**4 + 63358 * xv * xv + 83521)
a0, a2, a4, a6 = f[0], f[2], f[4], f[6]

E1 = EllipticCurve([0, a4, 0, a2 * a6, a0 * a6**2])
E2 = EllipticCurve([0, a2, 0, a0 * a4, a0**2 * a6])
rep = {'E1_ainvs': [str(t) for t in E1.a_invariants()],
       'E2_ainvs': [str(t) for t in E2.a_invariants()]}

# sanity: phi1, phi2 land on E1, E2 as formal identities in x
Fr = FractionField(R)
yy = Fr(f)                                    # y^2
X1, Y2_1 = Fr(a6 * xv**2), Fr(a6**2) * yy     # (a6 x^2, (a6 y)^2)
rep['phi1_identity'] = bool(Y2_1 == X1**3 + a4 * X1**2 + a2 * a6 * X1 + a0 * a6**2)
X2, Y2_2 = Fr(a0) / xv**2, Fr(a0**2) * yy / xv**6
rep['phi2_identity'] = bool(Y2_2 == X2**3 + a2 * X2**2 + a0 * a4 * X2 + a0**2 * a6)

rep['E1_rank_pari'] = str(pari(E1).ellrank())
rep['E2_rank_pari'] = str(pari(E2).ellrank())
rep['E1_torsion'] = str(E1.torsion_subgroup().invariants())
rep['E2_torsion'] = str(E2.torsion_subgroup().invariants())
rep['E1_good_ordinary_at_11'] = bool(E1.has_good_reduction(p) and E1.is_ordinary(p))
rep['E2_good_ordinary_at_11'] = bool(E2.has_good_reduction(p) and E2.is_ordinary(p))
rep['E1_bad_primes'] = [int(q) for q in E1.conductor().prime_factors()]
rep['E2_bad_primes'] = [int(q) for q in E2.conductor().prime_factors()]


def parse(s):
    return [K(sage_eval(t, locals={'O': O, 'x': xv})) for t in s[1:-1].split(' : ')]


def affine_x_y(s):
    """Return the affine (x, y) in Qp, possibly with v(x) < 0."""
    X, Y, Z = parse(s)
    if Z == 0:
        # projective chart (1 : y/x^3 : 0)-style: the point lies over x = infinity
        return None, None
    X, Y = X / Z, Y / Z**3
    return X, Y


rows = []
for i, row in enumerate(res30['extra_points_by_omega']):
    for s in row:
        X, Y, Z = parse(s)
        if Z != 0:
            X, Y = X / Z, Y / Z**3
        vx = int(X.valuation())
        # phi1 and phi2 in the affine coordinates of E1, E2
        P1 = (a6 * X**2, a6 * Y)
        P2 = (a0 / X**2, a0 * Y / X**3)
        e1_res = 'nonsingular' if P1[0].valuation() >= 0 else 'infinity'
        e2_res = 'nonsingular' if P2[0].valuation() >= 0 else 'infinity'
        d1 = P1[1]**2 - (P1[0]**3 + a4 * P1[0]**2 + a2 * a6 * P1[0] + a0 * a6**2)
        d2 = P2[1]**2 - (P2[0]**3 + a2 * P2[0]**2 + a0 * a4 * P2[0] + a0**2 * a6)
        t = (1 + X) / (1 - X)
        rows.append({
            'omega_index': i,
            'raw_string': s,
            'affine_x_valuation': vx,
            'affine_y_valuation': int(Y.valuation()),
            'phi1_X_valuation': int(P1[0].valuation()),
            'phi1_Y_valuation': int(P1[1].valuation()),
            'phi2_X_valuation': int(P2[0].valuation()),
            'phi2_Y_valuation': int(P2[1].valuation()),
            'phi1_chart': e1_res,
            'phi2_chart': e2_res,
            'phi1_residual_valuation': (int(d1.valuation()) if d1 != 0
                                        else int(d1.precision_absolute())),
            'phi2_residual_valuation': (int(d2.valuation()) if d2 != 0
                                        else int(d2.precision_absolute())),
            't_valuation': int(t.valuation()),
            't_minus_1_valuation': int((t - 1).valuation()) if t != 1 else None,
            't_plus_1_valuation': int((t + 1).valuation()) if t != -1 else None,
            't_mod_11_12': int(ZZ(t.lift()) % p**12) if t.valuation() >= 0 else None,
        })

rep['min_phi1_residual'] = min(r['phi1_residual_valuation'] for r in rows)
rep['min_phi2_residual'] = min(r['phi2_residual_valuation'] for r in rows)
rep['phi1_charts'] = sorted({r['phi1_chart'] for r in rows})
rep['phi2_charts'] = sorted({r['phi2_chart'] for r in rows})
rep['affine_x_valuations'] = sorted({r['affine_x_valuation'] for r in rows})
rep['t_valuations'] = sorted({r['t_valuation'] for r in rows})

(BASE / 'enrich_candidates_sage.json').write_text(
    json.dumps({'report': rep, 'rows': rows}, indent=2, default=_j) + '\n')
print(json.dumps(rep, indent=2, default=_j))
