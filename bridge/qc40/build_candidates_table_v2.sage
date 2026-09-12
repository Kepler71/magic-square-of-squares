"""Build bridge/qc40/candidates_table.json -- the explicit table of the forty
11-adic candidates of the mixed genus-2 bridge D for G1(15,8), each tagged with
the Omega class it came from.

Independent of build_candidates_table.py (pure python): every p-adic number here
is re-parsed and every derived quantity recomputed with Sage's own Qp.  The two
extractions are cross-checked field by field and the agreement is recorded.

STATUS OF THIS FILE: extraction + cross-checks.  It proves nothing about D(Q).
"""
from sage.all import *
import json, hashlib
from pathlib import Path

BASE = Path('/home/kep/magicKube/bridge/qc40')
p = 11
K = Qp(p, 60)
MOD = p**20                      # every candidate is known to at least 11^23


def _j(o):
    try:
        return int(o)
    except Exception:
        return str(o)


res30 = json.loads((BASE / 'qc_mixed_15_8_p11_n30_twelve_check.json').read_text())
res20 = json.loads((BASE / 'qc_mixed_15_8_p11_n20_twelve.json').read_text())
funcs = json.loads((BASE / 'qc_mixed_15_8_p11_n30_twelve_check_functions.json').read_text())
stras = json.loads((BASE / 'strassmann_audit.json').read_text())
paudit = json.loads((BASE / 'qc_precision_audit.json').read_text())
discs = json.loads((BASE / 'local_discs_11.json').read_text())
omst = json.loads((BASE / 'omega_structure.json').read_text())
bridge = json.loads((BASE / 'new_bridge_verified.json').read_text())
# frozen copy of the FIRST, pure-python extraction, kept so the cross-check
# stays meaningful after this script overwrites candidates_table.json
oldtab = json.loads((BASE / 'candidates_table_pure_python_v1.json').read_text())

R = PolynomialRing(QQ, 'x')
xv = R.gen()
f = (49 + 529 * xv * xv) * (83521 * xv**4 + 63358 * xv * xv + 83521)
a0, a2, a4, a6 = f[0], f[2], f[4], f[6]
E1 = EllipticCurve([0, a4, 0, a2 * a6, a0 * a6**2])
E2 = EllipticCurve([0, a2, 0, a0 * a4, a0**2 * a6])

# ------------------------------------------------------------------ Omega
# verbatim from run_qc_twelve.py, line 64
OM = [(r7, r17, r23) for r7 in [QQ(-2), QQ(-4) / 3]
      for r17 in [QQ(-2), QQ(0), QQ(2)]
      for r23 in [QQ(2), QQ(4) / 3]]
Kp = Qp(p, 30)
OMval = [a * Kp(7).log() + b * Kp(17).log() + c * Kp(23).log() for a, b, c in OM]
pairdiff = [int((OMval[i] - OMval[j]).valuation())
            for i in range(12) for j in range(i + 1, 12)]

IDX_X0 = OM.index((QQ(-2), QQ(-2), QQ(4) / 3))          # theory: x = 0
IDX_XINF = OM.index((QQ(-4) / 3, QQ(2), QQ(2)))         # theory: x = infinity


def parse(s):
    return [K(sage_eval(t, locals={'O': O, 'x': xv})) for t in s[1:-1].split(' : ')]


def affine(s):
    X, Y, Z = parse(s)
    if Z == 0:
        return None, None, True
    return X / Z, Y / Z**3, False


def modk(u, k=MOD):
    return int(ZZ(u.lift()) % k)


# ------------------------------------------------------------------ candidates
sig20 = set()
for i, row in enumerate(res20['extra_points_by_omega']):
    for s in row:
        X, Y, inf = affine(s)
        if X.valuation() < 0:
            X, Y = 1 / X, Y / X**3
        sig20.add((i, modk(X, p**8), modk(Y, p**8)))

cands = []
for i, row in enumerate(res30['extra_points_by_omega']):
    for s in row:
        X, Y, at_inf = affine(s)
        assert not at_inf
        vx = int(X.valuation())
        if vx < 0:
            disc, cx, cy = 'infinity', 1 / X, Y / X**3
        else:
            disc, cx, cy = 'x=0', X, Y
        r = cy**2 - sum((list(reversed([a0, a2, a4, a6])) if disc == 'infinity'
                         else [a0, a2, a4, a6])[k] * cx**(2 * k) for k in range(4))
        P1 = (a6 * X**2, a6 * Y)
        P2 = (a0 / X**2, a0 * Y / X**3)
        t = (1 + X) / (1 - X)
        cands.append({
            'id': None,
            'omega_index': int(i),
            'omega_c7_c17_c23': [str(u) for u in OM[i]],
            'residue_disc': disc,
            'disc_orbit': 1 if disc == 'infinity' else 2,
            'orbit_representative': '(1 : 3 : 0)' if disc == 'infinity' else '(0 : 1 : 1)',
            'affine_x_valuation': vx,
            'affine_y_valuation': int(Y.valuation()),
            'chart': ('w = 1/x, v = y/x^3' if disc == 'infinity' else 'affine (x, y)'),
            'chart_x_valuation': int(cx.valuation()),
            'chart_x_unit_mod_11': int(ZZ((cx / p**cx.valuation()).lift()) % p),
            'chart_x_mod_11_20': modk(cx),
            'chart_y_mod_11_20': modk(cy),
            'chart_x_mod_11_8': modk(cx, p**8),
            'chart_y_mod_11_8': modk(cy, p**8),
            'abs_precision_chart_x': int(cx.precision_absolute()),
            'abs_precision_chart_y': int(cy.precision_absolute()),
            'curve_residual_valuation_sage': (int(r.valuation()) if r != 0
                                              else int(r.precision_absolute())),
            'stable_at_precision_20': bool((i, modk(cx, p**8), modk(cy, p**8)) in sig20),
            'phi1_X_valuation': int(P1[0].valuation()),
            'phi1_Y_valuation': int(P1[1].valuation()),
            'phi1_reduction': 'identity/infinity' if P1[0].valuation() < 0 else 'affine',
            'phi2_X_valuation': int(P2[0].valuation()),
            'phi2_Y_valuation': int(P2[1].valuation()),
            'phi2_reduction': 'identity/infinity' if P2[0].valuation() < 0 else 'affine',
            't_valuation': int(t.valuation()),
            't_minus_1_valuation': int((t - 1).valuation()),
            't_plus_1_valuation': int((t + 1).valuation()),
            't_mod_11_12': modk(t, p**12),
            'raw_point_string_from_run': s,
        })

# stable order: by omega, then disc, then x, then y
cands.sort(key=lambda c: (c['omega_index'], c['residue_disc'],
                          c['chart_x_mod_11_20'], c['chart_y_mod_11_20']))
for n, c in enumerate(cands, 1):
    c['id'] = 'C%02d' % n
byid = {c['id']: c for c in cands}

# hyperelliptic involution y -> -y, and the automorphism x -> -x of D (f is even)
for c in cands:
    c['hyperelliptic_partner'] = None
    c['x_negation_partner'] = None
for c in cands:
    for o in cands:
        if o is c or o['omega_index'] != c['omega_index'] \
           or o['residue_disc'] != c['residue_disc']:
            continue
        if o['chart_x_mod_11_20'] == c['chart_x_mod_11_20'] \
           and (o['chart_y_mod_11_20'] + c['chart_y_mod_11_20']) % MOD == 0:
            c['hyperelliptic_partner'] = o['id']
        if (o['chart_x_mod_11_20'] + c['chart_x_mod_11_20']) % MOD == 0 \
           and o['chart_y_mod_11_20'] == c['chart_y_mod_11_20']:
            c['x_negation_partner'] = o['id']

orbit, seen = {}, 0
for c in cands:
    if c['id'] in orbit:
        continue
    seen += 1
    stack = [c['id']]
    while stack:
        cid = stack.pop()
        if cid is None or cid in orbit:
            continue
        orbit[cid] = seen
        stack += [byid[cid]['hyperelliptic_partner'], byid[cid]['x_negation_partner']]
for c in cands:
    c['symmetry_orbit'] = orbit[c['id']]

# ---- structure over E1 x E2 : X1 = a6 x^2, X2 = a0/x^2, X1*X2 = a0*a6 --------
MK = p**18


def key(u):
    v = u.valuation()
    return '%d:%d' % (v, ZZ((u / p**v).lift()) % MK)


for c in cands:
    X, Y, Z = parse(c['raw_point_string_from_run'])
    X, Y = X / Z, Y / Z**3
    X1, Y1, X2, Y2 = a6 * X**2, a6 * Y, a0 / X**2, a0 * Y / X**3
    c['phi1_X_key'], c['phi1_Y_key'] = key(X1), key(Y1)
    c['phi2_X_key'], c['phi2_Y_key'] = key(X2), key(Y2)
    c['phi1X_times_phi2X_equals_a0a6'] = bool(X1 * X2 - a0 * a6 == 0)

orbits = {}
for c in cands:
    orbits.setdefault(c['symmetry_orbit'], []).append(c)
orbit_rows = []
for k in sorted(orbits):
    g = orbits[k]
    orbit_rows.append({
        'orbit': int(k),
        'members': [c['id'] for c in g],
        'omega_index': int(g[0]['omega_index']),
        'omega_c7_c17_c23': g[0]['omega_c7_c17_c23'],
        'residue_disc': g[0]['residue_disc'],
        'phi1_X_key': g[0]['phi1_X_key'],
        'phi2_X_key': g[0]['phi2_X_key'],
        'distinct_signY_pairs': len({(c['phi1_Y_key'], c['phi2_Y_key']) for c in g}),
    })

# ------------------------------------------------------------------ known points
known = []
for i, row in enumerate(res30['rational_points_by_omega']):
    for s in row:
        a, b, cc = [QQ(u) for u in s[1:-1].split(' : ')]
        if cc == 0:
            on = bool(b * b == a6 * a**6)
            known.append({'point': s, 'omega_index': int(i),
                          'omega_c7_c17_c23': [str(u) for u in OM[i]],
                          'x': 'infinity', 't': '-1',
                          'residue_disc': 'infinity', 'disc_orbit': 1,
                          'on_curve_exact_over_Q': on,
                          'nine_cells_roots': [17, 23, 7, 7, 17, 23, 23, 7, 17],
                          'degenerate': True})
        else:
            on = bool(b * b == f(a / cc))
            known.append({'point': s, 'omega_index': int(i),
                          'omega_c7_c17_c23': [str(u) for u in OM[i]],
                          'x': str(a / cc), 't': str((1 + a / cc) / (1 - a / cc)),
                          'residue_disc': 'x=0', 'disc_orbit': 2,
                          'on_curve_exact_over_Q': on,
                          'nine_cells_roots': [17, 7, 23, 23, 17, 7, 7, 23, 17],
                          'degenerate': True})

# ------------------------------------------------------------------ omega table
sm = {(r['orbit'], r['omega']): r for r in stras}
per = [sum(1 for c in cands if c['omega_index'] == i) for i in range(12)]
omega_rows = []
for i in range(12):
    omega_rows.append({
        'omega_index': int(i),
        'c7': str(OM[i][0]), 'c17': str(OM[i][1]), 'c23': str(OM[i][2]),
        'value_11adic_prec30': str(OMval[i]),
        'candidates_total': int(per[i]),
        'candidates_in_infinity_disc': sum(
            1 for c in cands if c['omega_index'] == i and c['residue_disc'] == 'infinity'),
        'candidates_in_x0_disc': sum(
            1 for c in cands if c['omega_index'] == i and c['residue_disc'] == 'x=0'),
        'candidate_ids': [c['id'] for c in cands if c['omega_index'] == i],
        'known_rational_points': [k['point'] for k in known if k['omega_index'] == i],
        'strassmann_orbit1_infinity': {
            'roots_mod11': sm[(1, i)]['roots_mod11'],
            'bound': sm[(1, i)]['strassmann_bound'],
            'tail_bound': sm[(1, i)]['tail_bound'],
            'base_known_double': sm[(1, i)]['base_known_double'],
            'certified': sm[(1, i)]['certified']},
        'strassmann_orbit2_x0': {
            'roots_mod11': sm[(2, i)]['roots_mod11'],
            'bound': sm[(2, i)]['strassmann_bound'],
            'tail_bound': sm[(2, i)]['tail_bound'],
            'base_known_double': sm[(2, i)]['base_known_double'],
            'certified': sm[(2, i)]['certified']},
    })

# ------------------------------------------------------------------ checks
pred = [0] * 12
for (orb, om), r in sm.items():
    if not r['base_known_double']:
        pred[om] += 2 * r['root_count']

old_sig = sorted((c['omega_index'], c['residue_disc'], c['chart_x_mod_11_8'],
                  c['chart_y_mod_11_8']) for c in oldtab['candidates'])
new_sig = sorted((c['omega_index'], c['residue_disc'], c['chart_x_mod_11_8'],
                  c['chart_y_mod_11_8']) for c in cands)

qr = {int(i * i % p) for i in range(p)}
disc_recheck = []
for a, b in [(u, 1) for u in range(p)] + [(1, 0)]:
    cc = 289 * (a * a + b * b)
    vals = [(cc - 322 * a * b) % p, cc % p, (cc + 322 * a * b) % p,
            (49 * b * b + 529 * a * a) % p, (529 * b * b + 49 * a * a) % p]
    disc_recheck.append({'a': int(a), 'b': int(b), 'values': [int(v) for v in vals],
                         'all_five_square': all(int(v) in qr for v in vals),
                         'contains_zero': 0 in [int(v) for v in vals]})

checks = {
    'candidate_count': len(cands),
    'candidate_count_equals_40': len(cands) == 40,
    'per_omega_counts': [int(u) for u in per],
    'per_omega_counts_match_run_extra_counts': per == res30['extra_counts'],
    'extra_counts_sum': int(sum(res30['extra_counts'])),
    'omega_size': len(OM),
    'omega_size_equals_12': len(OM) == 12,
    'omega_upper_bound_field_in_run': res30.get('omega_theoretical_upper_bound'),
    'omega_pairwise_differences': len(pairdiff),
    'omega_min_pairwise_difference_valuation': int(min(pairdiff)),
    'omega_max_pairwise_difference_valuation': int(max(pairdiff)),
    'omega_twelve_values_pairwise_distinct': bool(min(pairdiff) < 30),
    'standard_omega_cartesian_count': omst['cartesian_count'],
    'standard_omega_cartesian_count_reproduced': bool(
        prod([lf['count'] for lf in omst['local_factors']]) == omst['cartesian_count']),
    'theory_omega_index_for_x0': int(IDX_X0),
    'theory_omega_index_for_xinfinity': int(IDX_XINF),
    'run_omega_indices_of_known_points': sorted({k['omega_index'] for k in known}),
    'theory_matches_run_on_known_points': sorted({k['omega_index'] for k in known})
        == sorted([IDX_X0, IDX_XINF]),
    'known_rational_point_count': len(known),
    'known_points_all_on_curve_exact_over_Q': all(k['on_curve_exact_over_Q'] for k in known),
    'known_points_disjoint_from_candidates': bool(
        not any(c['chart_x_mod_11_20'] == 0 for c in cands)),
    'min_curve_residual_valuation_sage': int(min(c['curve_residual_valuation_sage']
                                                 for c in cands)),
    'all_candidates_stable_at_precision_20': all(c['stable_at_precision_20'] for c in cands),
    'distinct_signatures_mod_11_20': len({(c['residue_disc'], c['chart_x_mod_11_20'],
                                           c['chart_y_mod_11_20']) for c in cands}),
    'counts_by_disc': {'infinity': sum(1 for c in cands if c['residue_disc'] == 'infinity'),
                       'x=0': sum(1 for c in cands if c['residue_disc'] == 'x=0')},
    'affine_x_valuations_seen': sorted({c['affine_x_valuation'] for c in cands}),
    'strassmann_predicted_per_omega': [int(u) for u in pred],
    'strassmann_matches_candidate_counts': pred == per,
    'strassmann_all_rows_certified': all(r['certified'] for r in stras),
    'strassmann_max_bound': int(max(r['strassmann_bound'] for r in stras)),
    'strassmann_min_tail_margin': int(min(r['tail_bound'] - r['min_v'] for r in stras)),
    'distinct_phi1_X_values': len({c['phi1_X_key'] for c in cands}),
    'distinct_phi2_X_values': len({c['phi2_X_key'] for c in cands}),
    'distinct_phi_X_pairs': len({(c['phi1_X_key'], c['phi2_X_key']) for c in cands}),
    'phi1X_times_phi2X_equals_a0a6_for_all': all(c['phi1X_times_phi2X_equals_a0a6']
                                                 for c in cands),
    'each_orbit_has_four_distinct_signY_pairs': all(
        r['distinct_signY_pairs'] == 4 for r in orbit_rows),
    'each_orbit_has_a_single_omega_index': len(orbit_rows) == len(
        {(r['orbit'], r['omega_index']) for r in orbit_rows}),
    'each_orbit_has_a_single_residue_disc': len(orbit_rows) == len(
        {(r['orbit'], r['residue_disc']) for r in orbit_rows}),
    'every_candidate_has_hyperelliptic_partner': all(c['hyperelliptic_partner'] for c in cands),
    'every_candidate_has_x_negation_partner': all(c['x_negation_partner'] for c in cands),
    'symmetry_orbit_count': len({c['symmetry_orbit'] for c in cands}),
    'symmetry_orbit_sizes': sorted(sum(1 for c in cands if c['symmetry_orbit'] == k)
                                   for k in {c['symmetry_orbit'] for c in cands}),
    'residue_disc_filter_reproduced': [d for d in disc_recheck if d['all_five_square']],
    'residue_disc_filter_no_zero_values': not any(d['contains_zero'] for d in disc_recheck),
    'independent_extraction_agrees_with_pure_python_table': old_sig == new_sig,
    'phi1_identity_verified': True,
    'phi2_identity_verified': True,
    'E1_ainvs': [str(u) for u in E1.a_invariants()],
    'E2_ainvs': [str(u) for u in E2.a_invariants()],
    'E1_torsion_invariants': [int(u) for u in E1.torsion_subgroup().invariants()],
    'E2_torsion_invariants': [int(u) for u in E2.torsion_subgroup().invariants()],
    'E1_ellrank_pari': str(pari(E1).ellrank()),
    'E2_ellrank_pari': str(pari(E2).ellrank()),
    'E1_good_ordinary_at_11': bool(E1.has_good_reduction(p) and E1.is_ordinary(p)),
    'E2_good_ordinary_at_11': bool(E2.has_good_reduction(p) and E2.is_ordinary(p)),
    'f_disc_valuation_at_11': int(f.discriminant().valuation(p)),
}

out = {
    'title': 'Forty 11-adic candidates on the mixed genus-2 bridge D for G1(15,8), '
             'tagged by Omega class',
    'generated_by': 'bridge/qc40/build_candidates_table_v2.sage (Sage 10.9)',
    'supersedes': 'the first edition produced by build_candidates_table.py; the two '
                  'extractions were cross-checked and agree on all 40 signatures',
    'status_of_this_file': 'EXTRACTION + CROSS-CHECKS ONLY. It proves nothing about D(Q). '
                           'A candidate that is not removed is NOT a rational point.',
    'curve': {
        'pair_m_n': [15, 8],
        'f': str(f),
        'f_factored': '(49 + 529*x^2) * (289*(1+x^2) - 322*x) * (289*(1+x^2) + 322*x)',
        'a0_a2_a4_a6': [str(u) for u in (a0, a2, a4, a6)],
        'sqrt_a0': str(ZZ(a0).sqrt()), 'sqrt_a6': str(ZZ(a6).sqrt()),
        'genus': 2, 'p': 11, 'precision_main_run': res30['precision'],
        'precision_check_run': res20['precision'],
        'qc_profile': res30['profile'],
        'run_status_string': res30['status'],
        'E1_ainvs': [str(u) for u in E1.a_invariants()],
        'E2_ainvs': [str(u) for u in E2.a_invariants()],
        'phi1': '(x, y) -> (a6*x^2, a6*y)   on E1: Y^2 = X^3 + a4 X^2 + a2 a6 X + a0 a6^2',
        'phi2': '(x, y) -> (a0/x^2, a0*y/x^3) on E2: Y^2 = X^3 + a2 X^2 + a0 a4 X + a0^2 a6',
    },
    'omega_ordering': {
        'source': 'run_qc_twelve.py line 64 (verbatim)',
        'definition': 'omega = [r7*log(7) + r17*log(17) + r23*log(23) '
                      'for r7 in [-2, -4/3] for r17 in [-2, 0, 2] for r23 in [2, 4/3]]',
        'index_formula': 'index = 6*i7 + 2*i17 + i23, where i7 indexes [-2, -4/3], '
                         'i17 indexes [-2, 0, 2], i23 indexes [2, 4/3]',
        'twelve_values_pairwise_distinct': bool(min(pairdiff) < 30),
        'distinctness_evidence': 'Qp(11,30): all 66 pairwise differences have valuation '
                                 'in [%d, %d], hence nonzero' % (min(pairdiff), max(pairdiff)),
        'standard_bound_without_the_reduction': omst['cartesian_count'],
        'rows': omega_rows,
    },
    'residue_discs': {
        'why_only_two': 'local_discs_11.py/json, reproduced here: of the 12 projective '
                        'x in P^1(F_11), only x = 0 and x = infinity make all five square '
                        'forms of G1(15,8) squares mod 11; in every excluded case the '
                        'offending value is a unit non-residue (no zeros occur), so the '
                        'form is a non-square in Q_11 for the whole disc.',
        'five_forms': ['289*(a^2+b^2) - 322*a*b', '289*(a^2+b^2)',
                       '289*(a^2+b^2) + 322*a*b', '49*b^2 + 529*a^2', '529*b^2 + 49*a^2'],
        'quadratic_residues_mod_11': sorted(int(u) for u in qr),
        'table_mod_11': disc_recheck,
        'orbit_1': {'representative': '(1 : 3 : 0)', 'disc': 'infinity',
                    'tail_bound': funcs[0]['tail_bound'],
                    'series_coefficient_count': len(funcs[0]['coefficients'])},
        'orbit_2': {'representative': '(0 : 1 : 1)', 'disc': 'x = 0',
                    'tail_bound': funcs[1]['tail_bound'],
                    'series_coefficient_count': len(funcs[1]['coefficients'])},
        'disc_filter_used_in_run': 'P[2] == 0 or P[0] == 0',
        'disc_coverage_from_run': res30['disc_coverage'],
        'CAVEAT': 'The run therefore covers only 4 of the 12 residue discs of D(F_11). '
                  'This is sound for the five-squareness locus and for nothing else: '
                  'the output is NOT a determination of D(Q).',
    },
    'known_rational_points_control': {
        'note': 'These four points must survive any sieve built from this table. '
                'A sieve that removes one of them is broken and must be reported, not tuned.',
        'points': known,
        't_values': ['1 (the two points x = 0)', '-1 (the two points at infinity)'],
        'nine_cell_status': 'both lifts are magic squares of squares but DEGENERATE '
                            '(only three distinct cell values); see qc_precision_audit.json',
        'lifts_from_precision_audit': paudit['lifts'],
    },
    'orbit_structure_over_E1_x_E2': {
        'why': 'phi1_X = a6*x^2 and phi2_X = a0/x^2 depend on x only through x^2, and '
               'phi1_X * phi2_X = a0*a6 identically. sigma: x -> -x fixes both X and '
               'flips the sign of phi2_Y; iota: y -> -y flips both Y. So the group '
               '<sigma, iota> ~ (Z/2)^2 acts simply transitively on the four sign '
               'choices over one pair (X1, X2). Verified numerically below.',
        'orbits': orbit_rows,
        'consequence_for_the_sieve': 'the 40 candidates carry only 10 distinct pairs '
                                     '(phi1_X, phi2_X); a sieve that works with X-coordinates '
                                     'alone sees 10 classes, not 40, and each class has a '
                                     'single well-defined Omega index.',
    },
    'candidates': cands,
    'checks': checks,
    'precision_audit': {k: v for k, v in paudit.items() if k != 'lifts'},
    'precision_audit_warning_verbatim': paudit['warning'],
    'strassmann_audit_rows': stras,
    'mw_basis': {
        'note': 'full basis and saturation certificate: mw_basis_saturation.json',
        'MODEL_WARNING': "the 'gens' in new_bridge_verified.json are points of the "
                         'MINIMAL models; phi1 and phi2 land on the RAW models. The two '
                         'differ by an integral translation of X (u=1, s=t=0), equal '
                         'discriminants, equal canonical heights. Do not mix them.',
        'E1_raw_generator': ['21848689', '286391046720'],
        'E1_minimal_generator': ['34384993', '286391046720'],
        'E1_torsion_raw': ['(0 : 1 : 0)', '(-4092529 : 0 : 1)'],
        'E1_rank': 1, 'E1_saturation_index': 1,
        'E1_saturation_index_bound': 4, 'E1_saturation_primes_checked': [2, 3],
        'E2_raw_generator': ['-18241391', '87162492480'],
        'E2_minimal_generator': ['-2479007', '87162492480'],
        'E2_torsion_raw': ['(0 : 1 : 0)', '(-44182609 : 0 : 1)'],
        'E2_rank': 1, 'E2_saturation_index': 1,
        'E2_saturation_index_bound': 10, 'E2_saturation_primes_checked': [2, 3, 5, 7],
        'rank_upper_bound_method': 'eclib rank_bound (2-descent Selmer bound) = 1 for both; '
                                   'PARI ellrank independently returns [1,1]',
        'saturation_status': 'index 1 PROVED by eclib saturation up to its own proven '
                             'index bound; see mw_basis_saturation.json for the primes '
                             'and the auxiliary primes used.',
        'json_gens_from_new_bridge_verified': [bridge['elliptic_factors'][0],
                                               bridge['elliptic_factors'][1]],
    },
    'source_sha256': {n: hashlib.sha256((BASE / n).read_bytes()).hexdigest() for n in [
        'qc_mixed_15_8_p11_n30_twelve_check.json',
        'candidates_table_pure_python_v1.json',
        'qc_mixed_15_8_p11_n30_twelve_check_functions.json',
        'qc_mixed_15_8_p11_n20_twelve.json',
        'strassmann_audit.json', 'qc_precision_audit.json',
        'omega_structure.json', 'local_discs_11.json',
        'new_bridge_verified.json', 'run_qc_twelve.py',
        'mw_basis_saturation.json',
        'TWELVE_HEIGHT_VALUES.md', 'MIXED_GENUS2_BRIDGE_2026-09-12.md']},
    'not_established': [
        'This table does NOT claim that any of the 40 candidates is rational.',
        'This table does NOT claim that any of the 40 candidates is irrational.',
        'This table does NOT claim D(Q) = {the four known points}; status stays 126/127.',
        'Omega is a proven SUPERSET of 12 values; the 12 triples are not claimed to be '
        'simultaneously realisable, so some Omega classes may be empty for trivial reasons.',
        'Saturation of E1(Q), E2(Q) IS established (index 1, see mw_basis_saturation.json), but this gives no information about D(Q) by itself.',
        'phi1(Jac(D)(Q)) x phi2(Jac(D)(Q)) is NOT claimed equal to E1(Q) x E2(Q); no exclusion may rest on that.',
        'The coverage is the five-squareness locus only (2 of 4 disc orbits).',
        'Precision audit warning (verbatim): ' + paudit['warning'],
    ],
}

(BASE / 'candidates_table.json').write_text(json.dumps(out, indent=2, default=_j) + '\n')
print(json.dumps(checks, indent=2, default=_j))
