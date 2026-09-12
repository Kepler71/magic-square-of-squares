"""Audit of the E1(Q), E2(Q) data in new_bridge_verified.json.

Point of this script: the maps phi1 = (a6 x^2, a6 y), phi2 = (a0/x^2, a0 y/x^3)
land on the RAW models
    E1 : Y^2 = X^3 + a4 X^2 + a2 a6 X + a0 a6^2
    E2 : Y^2 = X^3 + a2 X^2 + a0 a4 X + a0^2 a6
but the 'gens' recorded in new_bridge_verified.json are points of the MINIMAL
models.  They are the same points, related by an integral translation of X;
this script pins that down so nobody feeds a minimal-model point to phi_i.
"""
from sage.all import *
import json
from pathlib import Path

BASE = Path('/home/kep/magicKube/bridge/qc40')
br = json.loads((BASE / 'new_bridge_verified.json').read_text())
R = PolynomialRing(QQ, 'x')
x = R.gen()
f = (49 + 529 * x * x) * (83521 * x**4 + 63358 * x * x + 83521)
a0, a2, a4, a6 = f[0], f[2], f[4], f[6]

raw = {'E1': EllipticCurve([0, a4, 0, a2 * a6, a0 * a6**2]),
       'E2': EllipticCurve([0, a2, 0, a0 * a4, a0**2 * a6])}
pari_pt = {'E1': (21848689, 286391046720), 'E2': (-18241391, 87162492480)}

out = {'note': 'gens in new_bridge_verified.json are on the MINIMAL models; '
               'phi1/phi2 land on the RAW models. Same points, X-translation only.'}
for k, C in raw.items():
    idx = 0 if k == 'E1' else 1
    fac = br['elliptic_factors'][idx]
    Cmin = C.minimal_model()
    iso = C.isomorphism_to(Cmin)
    u, r, s, t = iso.tuple()
    gj = [QQ(v) for v in fac['gens'][0]]
    row = {
        'raw_ainvs': [str(v) for v in C.a_invariants()],
        'minimal_ainvs': [str(v) for v in Cmin.a_invariants()],
        'json_minimal_ainvs': fac['minimal_ainvs'],
        'minimal_ainvs_agree': [str(v) for v in Cmin.a_invariants()] == fac['minimal_ainvs'],
        'iso_raw_to_minimal_urst': [str(v) for v in (u, r, s, t)],
        'iso_is_integral_X_translation': bool(u == 1 and s == 0 and t == 0 and r in ZZ),
        'disc_ratio_raw_over_minimal': str(C.discriminant() / Cmin.discriminant()),
        'raw_model_is_minimal_up_to_translation': bool(
            C.discriminant() == Cmin.discriminant()),
        'json_gen_on_minimal_model': bool(Cmin.is_on_curve(gj[0] / gj[2], gj[1] / gj[2])),
        'json_gen_on_raw_model': bool(C.is_on_curve(gj[0] / gj[2], gj[1] / gj[2])),
        'conductor': str(C.conductor()),
        'conductor_primes': [int(q) for q in C.conductor().prime_factors()],
    }
    Pmin = Cmin(gj[0] / gj[2], gj[1] / gj[2])
    isoback = Cmin.isomorphism_to(C)
    Praw = isoback(Pmin)
    row['json_gen_minimal'] = str(Pmin.xy())
    row['json_gen_pulled_back_to_raw'] = str(Praw.xy())
    row['pari_ellrank_point'] = [str(v) for v in pari_pt[k]]
    Q = C(*[QQ(v) for v in pari_pt[k]])
    row['pari_point_on_raw_model'] = True
    row['pari_point_equals_pulled_back_gen'] = bool(Q == Praw)
    row['canonical_height_raw'] = str(Praw.height())
    row['canonical_height_minimal'] = str(Pmin.height())
    row['heights_agree'] = bool(Praw.height() == Pmin.height())
    row['torsion_invariants'] = [int(v) for v in C.torsion_subgroup().invariants()]
    row['torsion_points_raw'] = [str(P) for P in C.torsion_points()]
    row['rank_pari_ellrank'] = str(pari(C).ellrank())
    row['json_rank'] = fac['rank']
    row['json_rank_bound'] = fac['rank_bound']
    row['json_certain'] = fac['certain']
    # saturation of the single recorded generator
    try:
        sat, index, reg = C.saturation([Praw])
        row['saturation_index'] = str(index)
        row['saturation_index_is_1'] = bool(index == 1)
        row['saturated_generator'] = str(sat[0].xy())
        row['saturation_succeeded'] = True
    except Exception as ex:
        row['saturation_succeeded'] = False
        row['saturation_error'] = repr(ex)
    try:
        G = C.gens()
        row['sage_gens_raw'] = [str(P.xy()) for P in G]
        row['sage_gens_count'] = len(G)
        row['sage_gen_equals_recorded'] = bool(
            len(G) == 1 and (G[0] == Praw or G[0] == -Praw))
    except Exception as ex:
        row['sage_gens_error'] = repr(ex)
    out[k] = row

(BASE / 'check_generators_sage.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
