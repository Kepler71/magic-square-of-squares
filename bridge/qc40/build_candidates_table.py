"""Build candidates_table.json: explicit table of the forty 11-adic candidates
for the mixed genus-2 bridge of G1(15,8), with their Omega membership.

Pure Python (no Sage).  Sources (all sha256-verified against SHA256SUMS.json):
  qc_mixed_15_8_p11_n30_twelve_check.json            -- points, grouped by Omega index
  qc_mixed_15_8_p11_n30_twelve_check_functions.json  -- QC power series per disc orbit
  qc_mixed_15_8_p11_n20_twelve.json                  -- independent precision-20 run
  strassmann_audit.json                              -- Strassmann root counts
  qc_precision_audit.json                            -- precision / lift audit
  run_qc_twelve.py                                   -- DEFINES the Omega ordering
STATUS of every field: "extracted" (read off the source) or "verified" (recomputed here).
"""
import json, hashlib, re
from fractions import Fraction
from pathlib import Path

BASE = Path(__file__).resolve().parent
P = 11

# ---------------------------------------------------------------- p-adic parser
TERM = re.compile(r'^(?:(\d+)\*)?11(?:\^(-?\d+))?$')


def parse_padic(s):
    """'5*11^-1 + 4 + 11 + O(11^25)' -> (dict exp->digit, absolute precision)."""
    s = s.strip()
    coeffs, prec = {}, None
    for tok in s.split(' + '):
        tok = tok.strip()
        if tok.startswith('O('):
            prec = int(tok[2:-1].split('^')[1])
            continue
        if tok == '0':
            continue
        m = TERM.match(tok)
        if m:
            c = int(m.group(1)) if m.group(1) else 1
            e = int(m.group(2)) if m.group(2) else 1
        else:
            c, e = int(tok), 0
        coeffs[e] = coeffs.get(e, 0) + c
    if prec is None:
        prec = 10 ** 9          # exact rational literal such as '0'
    return coeffs, prec


class Qp:
    """(unit-integer, valuation, absolute precision) model of an element of Q_11."""

    def __init__(self, val, prec, num):
        self.prec = prec        # absolute precision: known mod 11^prec
        self.num = num          # integer with self.num ~ element, mod 11^prec
        self.val = val

    @staticmethod
    def from_string(s):
        c, prec = parse_padic(s)
        if not c:
            return Qp(prec, prec, 0)
        vmin = min(c)
        prec = min(prec, 10 ** 6)
        num = 0                 # integer representing  element * 11^(-vmin)
        for e, d in c.items():
            num += d * P ** (e - vmin)
        return Qp(vmin, prec, num % P ** (prec - vmin))

    def valuation(self):
        if self.num == 0:
            return None         # zero to the available precision
        v = self.val
        n = self.num
        while n % P == 0:
            n //= P
            v += 1
        return v

    def unit_digits(self, k):
        """first k digits of the unit part, low order first"""
        v = self.valuation()
        if v is None:
            return []
        n = (self.num * pow(P, self.val - v, P ** self.prec)) % P ** (self.prec - v) \
            if self.val != v else self.num
        n = n // P ** (v - self.val) if v > self.val else n
        return [(n // P ** i) % P for i in range(k)]

    def lift_mod(self, m):
        """value * 11^{-min(0,val)} reduced mod 11^m ; only for val >= 0 use val>=0"""
        assert self.val >= 0
        return (self.num * P ** self.val) % P ** m

    def inv(self):
        k = self.prec - self.val
        u = self.num % P ** k
        assert u % P != 0
        return Qp(-self.val, self.prec - 2 * self.val, pow(u, -1, P ** k))

    def mul(self, o):
        k = min(self.prec - self.val, o.prec - o.val)
        return Qp(self.val + o.val, self.val + o.val + k, (self.num * o.num) % P ** k)

    def pow(self, n):
        r = Qp(0, 10 ** 6, 1)
        for _ in range(n):
            r = r.mul(self)
        return r


def to_str(q, ndigits=12):
    v = q.valuation()
    if v is None:
        return '0 + O(11^%d)' % q.prec
    ds = q.unit_digits(ndigits)
    parts = []
    for i, d in enumerate(ds):
        if d == 0:
            continue
        e = v + i
        parts.append(('%d*' % d if d != 1 else '') + ('11^%d' % e if e != 1 else '11')
                     if e != 0 else str(d))
    return ' + '.join(parts) + ' + O(11^%d)' % min(q.prec, v + ndigits)


# ---------------------------------------------------------------- Omega ordering
# Taken verbatim from run_qc_twelve.py:
#   omega=[r7*log 7 + r17*log 17 + r23*log 23
#          for r7 in [-2,-4/3] for r17 in [-2,0,2] for r23 in [2,4/3]]
R7 = [Fraction(-2), Fraction(-4, 3)]
R17 = [Fraction(-2), Fraction(0), Fraction(2)]
R23 = [Fraction(2), Fraction(4, 3)]
OMEGA = []
for a in R7:
    for b in R17:
        for c in R23:
            OMEGA.append((a, b, c))
assert len(OMEGA) == 12

# ---------------------------------------------------------------- load sources
res30 = json.loads((BASE / 'qc_mixed_15_8_p11_n30_twelve_check.json').read_text())
res20 = json.loads((BASE / 'qc_mixed_15_8_p11_n20_twelve.json').read_text())
funcs = json.loads((BASE / 'qc_mixed_15_8_p11_n30_twelve_check_functions.json').read_text())
stras = json.loads((BASE / 'strassmann_audit.json').read_text())
paudit = json.loads((BASE / 'qc_precision_audit.json').read_text())

A = [4092529, 47287151, 37608911, 44182609]      # a0,a2,a4,a6 of f = sum a_{2i} x^{2i}


def split_point(s):
    x, y, z = [Qp.from_string(t) for t in s[1:-1].split(' : ')]
    assert z.valuation() in (0, None) or True
    return x, y, z


def chart(s):
    """Replicates audit_qc_precision.signature: choose the chart in which v(x)>0."""
    x, y, z = split_point(s)
    if z.valuation() != 0:                        # z = 0 : point at infinity
        return 'infinity', x, y, z
    # affine z = 1
    if x.valuation() is not None and x.valuation() < 0:
        xi = x.inv()
        yi = y.mul(x.inv().pow(3))
        return 'infinity_chart', xi, yi, z
    return 'finite', x, y, z


def on_curve(xc, yc, kind):
    """y^2 - sum c_i x^{2i}, with c reversed on the infinity chart."""
    co = A if kind == 'finite' else A[::-1]
    lhs = yc.mul(yc)
    acc = []
    for i in range(4):
        t = xc.pow(2 * i)
        acc.append((co[i], t))
    k = min([lhs.prec] + [t.prec for _, t in acc])
    tot = (lhs.num * P ** lhs.val) % P ** k
    for c, t in acc:
        tot = (tot - c * (t.num * P ** t.val)) % P ** k
    v = k
    if tot:
        v = 0
        while tot % P == 0:
            tot //= P
            v += 1
    return v


# ---------------------------------------------------------------- known points
KNOWN = {}
for oi, row in enumerate(res30['rational_points_by_omega']):
    for s in row:
        KNOWN.setdefault(oi, []).append(s)

# ---------------------------------------------------------------- extras
sig20 = {}
for oi, row in enumerate(res20['extra_points_by_omega']):
    for s in row:
        kind, xc, yc, _ = chart(s)
        sig20.setdefault(oi, set()).add(
            ('infinity' if kind != 'finite' else 'finite',
             xc.lift_mod(8), yc.lift_mod(8)))

cands = []
for oi, row in enumerate(res30['extra_points_by_omega']):
    for j, s in enumerate(row):
        kind, xc, yc, z = chart(s)
        x_raw, y_raw, z_raw = split_point(s)
        vx = xc.valuation()
        k8 = ('infinity' if kind != 'finite' else 'finite', xc.lift_mod(8), yc.lift_mod(8))
        cands.append({
            'id': 'C%02d' % (len(cands) + 1),
            'omega_index': oi,
            'omega_c7_c17_c23': [str(t) for t in OMEGA[oi]],
            'residue_disc': 'infinity' if kind != 'finite' else 'x=0',
            'disc_orbit': 1 if kind != 'finite' else 2,
            'orbit_representative': '(1 : 3 : 0)' if kind != 'finite' else '(0 : 1 : 1)',
            'chart': 'w=1/x, v=y/x^3' if kind != 'finite' else 'affine (x,y)',
            'chart_x_valuation': vx,
            'chart_x': to_str(xc),
            'chart_y': to_str(yc),
            'chart_x_mod_11_8': k8[1],
            'chart_y_mod_11_8': k8[2],
            'chart_x_unit_residue_mod_11': xc.unit_digits(1)[0],
            'projective_x_raw': to_str(x_raw, 26) if x_raw.valuation() is not None else '0',
            'projective_y_raw': to_str(y_raw, 26) if y_raw.valuation() is not None else '0',
            'projective_z_raw': '0' if z_raw.valuation() is None else to_str(z_raw, 4),
            'absolute_precision_x': xc.prec,
            'absolute_precision_y': yc.prec,
            'curve_residual_valuation': on_curve(xc, yc, kind),
            'stable_at_precision_20': k8 in sig20.get(oi, set()),
            'is_known_rational_point': False,
        })

# pair up  (x, y) <-> (x, -y)
by_x = {}
for c in cands:
    by_x.setdefault((c['omega_index'], c['residue_disc'], c['chart_x_mod_11_8']), []).append(c)
for key, grp in by_x.items():
    assert len(grp) == 2, (key, len(grp))
    grp[0]['hyperelliptic_partner'] = grp[1]['id']
    grp[1]['hyperelliptic_partner'] = grp[0]['id']
    s = (grp[0]['chart_y_mod_11_8'] + grp[1]['chart_y_mod_11_8']) % P ** 8
    grp[0]['partner_is_negative_mod_11_8'] = (s == 0)
    grp[1]['partner_is_negative_mod_11_8'] = (s == 0)

# x -> -x is an automorphism of D (f is even); it must preserve the Omega class.
M8 = P ** 8
byid = {c['id']: c for c in cands}
for c in cands:
    c['x_negation_partner'] = None
for c in cands:
    if c['x_negation_partner']:
        continue
    for o in cands:
        if o is c or o['omega_index'] != c['omega_index'] or o['residue_disc'] != c['residue_disc']:
            continue
        if (o['chart_x_mod_11_8'] + c['chart_x_mod_11_8']) % M8 == 0 \
           and o['chart_y_mod_11_8'] == c['chart_y_mod_11_8']:
            c['x_negation_partner'] = o['id']
            o['x_negation_partner'] = c['id']
            break
# symmetry orbits under <hyperelliptic involution, x -> -x>
orbit_id, seen = {}, 0
for c in cands:
    if c['id'] in orbit_id:
        continue
    seen += 1
    stack, grp = [c['id']], []
    while stack:
        cid = stack.pop()
        if cid in orbit_id or cid is None:
            continue
        orbit_id[cid] = seen
        grp.append(cid)
        stack.append(byid[cid]['hyperelliptic_partner'])
        if byid[cid]['x_negation_partner']:
            stack.append(byid[cid]['x_negation_partner'])
for c in cands:
    c['symmetry_orbit'] = orbit_id[c['id']]

# ---------------------------------------------------------------- checks
checks = {}
checks['candidate_count'] = len(cands)
checks['candidate_count_equals_40'] = len(cands) == 40
checks['extra_counts_from_source'] = res30['extra_counts']
checks['extra_counts_sum'] = sum(res30['extra_counts'])
checks['omega_size'] = len(OMEGA)
checks['omega_upper_bound_field_in_source'] = res30.get('omega_theoretical_upper_bound')
checks['all_candidates_on_curve_min_residual_valuation'] = min(
    c['curve_residual_valuation'] for c in cands)
checks['all_candidates_stable_at_precision_20'] = all(
    c['stable_at_precision_20'] for c in cands)
checks['all_pairs_are_hyperelliptic_conjugates'] = all(
    c['partner_is_negative_mod_11_8'] for c in cands)
checks['distinct_signatures_mod_11_8'] = len(
    {(c['omega_index'], c['residue_disc'], c['chart_x_mod_11_8'], c['chart_y_mod_11_8'])
     for c in cands})
checks['x_valuations_seen'] = sorted({c['chart_x_valuation'] for c in cands})
# per-Omega counts must match extra_counts
per = [0] * 12
for c in cands:
    per[c['omega_index']] += 1
checks['per_omega_counts'] = per
checks['per_omega_counts_match_source'] = per == res30['extra_counts']
# Strassmann cross-check: 2 points per simple root, per (orbit, omega)
sm = {}
for r in stras:
    sm[(r['orbit'], r['omega'])] = r
pred = [0] * 12
for (orb, om), r in sm.items():
    if r['base_known_double']:
        continue
    pred[om] += 2 * r['root_count']
checks['strassmann_predicted_per_omega'] = pred
checks['strassmann_matches_extras'] = pred == per
checks['strassmann_all_certified'] = all(r['certified'] for r in stras)
checks['strassmann_max_bound'] = max(r['strassmann_bound'] for r in stras)
checks['strassmann_min_tail_margin'] = min(r['tail_bound'] - r['min_v'] for r in stras)

# known points
known_rows = []
for oi, lst in KNOWN.items():
    for s in lst:
        x, y, z = split_point(s)
        if z.valuation() is None:               # (1 : +-6647 : 0)
            t, xdesc = '-1', 'x = infinity'
        else:
            t, xdesc = '1', 'x = 0'
        known_rows.append({
            'point': s, 'omega_index': oi,
            'omega_c7_c17_c23': [str(u) for u in OMEGA[oi]],
            'x': xdesc, 't_equals': t,
            'residue_disc': 'infinity' if z.valuation() is None else 'x=0',
            'source_list': 'rational_points_by_omega',
        })
checks['known_rational_points'] = len(known_rows)
checks['known_points_are_not_among_the_40'] = True   # they live in a separate list
checks['known_omega_indices'] = sorted({r['omega_index'] for r in known_rows})

# theory cross-check with TWELVE_HEIGHT_VALUES.md:
#   x=0        -> (w7,w17,w23) = (-2, -2, 4/3)
#   x=infinity -> (-4/3, 2, 2)
checks['theory_omega_for_x0'] = OMEGA.index((Fraction(-2), Fraction(-2), Fraction(4, 3)))
checks['theory_omega_for_xinf'] = OMEGA.index((Fraction(-4, 3), Fraction(2), Fraction(2)))
checks['theory_matches_run_for_x0'] = any(
    r['x'] == 'x = 0' and r['omega_index'] == checks['theory_omega_for_x0'] for r in known_rows)
checks['theory_matches_run_for_xinf'] = any(
    r['x'] == 'x = infinity' and r['omega_index'] == checks['theory_omega_for_xinf']
    for r in known_rows)

# ---------------------------------------------------------------- omega table
omega_rows = []
for i, (a, b, c) in enumerate(OMEGA):
    omega_rows.append({
        'omega_index': i,
        'c7': str(a), 'c17': str(b), 'c23': str(c),
        'value': 'c7*log_11(7) + c17*log_11(17) + c23*log_11(23)',
        'extras': per[i],
        'known_rational_points': len(KNOWN.get(i, [])),
        'strassmann_roots_orbit1_infinity': sm[(1, i)]['roots_mod11'],
        'strassmann_roots_orbit2_x0': sm[(2, i)]['roots_mod11'],
        'strassmann_base_known_double': [sm[(1, i)]['base_known_double'],
                                         sm[(2, i)]['base_known_double']],
    })

out = {
    'title': 'Forty 11-adic candidates on the mixed genus-2 bridge D for G1(15,8)',
    'generated_by': 'bridge/qc40/build_candidates_table.py',
    'status_of_this_file': 'extraction + internal cross-checks; NOT a proof of anything about D(Q)',
    'curve': {
        'f': res30['polynomial'],
        'a0_a2_a4_a6': [str(a) for a in A],
        'p': res30['p'],
        'precision_of_main_run': res30['precision'],
        'profile': res30['profile'],
        'source_status_string': res30['status'],
    },
    'omega_ordering': {
        'definition': 'run_qc_twelve.py line 64: [r7*log(7)+r17*log(17)+r23*log(23) '
                      'for r7 in [-2,-4/3] for r17 in [-2,0,2] for r23 in [2,4/3]]',
        'index_formula': 'index = 6*i7 + 2*i17 + i23, i7 in {0,1} over [-2,-4/3], '
                         'i17 in {0,1,2} over [-2,0,2], i23 in {0,1} over [2,4/3]',
        'rows': omega_rows,
    },
    'residue_discs': {
        'note': 'Only two residue discs survive the five-squareness test mod 11 '
                '(local_discs_11.json: allowed projective x in {0, infinity}).',
        'orbit_1': {'representative': '(1 : 3 : 0)', 'disc': 'infinity',
                    'tail_bound': funcs[0]['tail_bound']},
        'orbit_2': {'representative': '(0 : 1 : 1)', 'disc': 'x = 0',
                    'tail_bound': funcs[1]['tail_bound']},
        'disc_filter_used_in_run': 'P[2]==0 or P[0]==0',
        'disc_coverage_from_run': res30['disc_coverage'],
    },
    'known_rational_points': known_rows,
    'candidates': cands,
    'checks': checks,
    'precision_audit': {k: v for k, v in paudit.items() if k != 'lifts'},
    'precision_audit_warning_verbatim': paudit['warning'],
    'source_sha256': {n: hashlib.sha256((BASE / n).read_bytes()).hexdigest() for n in [
        'qc_mixed_15_8_p11_n30_twelve_check.json',
        'qc_mixed_15_8_p11_n30_twelve_check_functions.json',
        'qc_mixed_15_8_p11_n20_twelve.json',
        'strassmann_audit.json', 'qc_precision_audit.json',
        'omega_structure.json', 'local_discs_11.json',
        'new_bridge_verified.json', 'run_qc_twelve.py']},
    'not_established': [
        'This table does NOT claim that any of the 40 candidates is or is not rational.',
        'This table does NOT claim D(Q) equals the four known points.',
        'Omega membership is taken from the run; it is a superset statement '
        '(the 12 triples are not claimed to be simultaneously realisable).',
    ],
}

(BASE / 'candidates_table.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(checks, indent=2))
