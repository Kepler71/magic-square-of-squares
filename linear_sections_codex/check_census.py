"""Check that every reduced slope 0<k<1 of denominator <=48 has a witness."""
from collections import Counter
from fractions import Fraction
from math import gcd
from pathlib import Path
import json

root = Path(__file__).resolve().parent
expected = {Fraction(a, b) for b in range(2, 49) for a in range(1, b) if gcd(a, b) == 1}
witnesses = {}

def add(k, method, file, **extra):
    k = Fraction(k)
    assert k not in witnesses, (k, "duplicate exclusion")
    witnesses[k] = dict(k=str(k), method=method, file=file, **extra)

name = 'small_slopes_verified.json'
small = json.loads((root / name).read_text())
for k in small['four_AP_exclusions']:
    add(k, 'four_square_AP', name)
add(small['repeated_cell_slope'], 'repeated_cells', name)
for row in small['rank_zero_exclusions']:
    assert row['pari_rank_bounds'] == [0, 0]
    assert row['eclib'] == dict(rank=0, upper_bound=0, certain=True)
    add(row['k'], 'cubic', name)
assert set(map(Fraction, small['all_slopes'])) == set(witnesses)

for name in ['extended_slopes_7_12.json', 'extended_slopes_13_24.json',
             'four_cell_followup.json', 'extended_slopes_25_48.json',
             'four_cell_25_48.json', 'exception_40_43_certified.json']:
    data = json.loads((root / name).read_text())
    for row in data['slopes']:
        if row['status'] != 'excluded':
            continue
        witness = row['tried'][row['witness_index']]
        assert data['rank_cache'][witness['curve_key']]['rank_bounds'] == [0, 0]
        cert = witness['certificate']
        special = cert.get('rank_proof') == 'eclib_selmer_plus_fisher_pairing'
        if special:
            assert row['k'] == '40/43'
            audit = json.loads((root / 'pairing_40_43_audit.json').read_text())
            assert audit['all_passed'] and audit['rank_upper_bound'] == 0
        else:
            assert witness['eclib'] == dict(rank=0, upper_bound=0, certain=True)
        assert not any(t['positive'] and t['distinct'] and len(t['square_indices']) == 9
                       for t in cert['all_p_tests'])
        method = 'quartic_fisher' if special else ('cubic' if len(witness['coefficients']) == 3 else 'quartic')
        add(row['k'], method, name, witness_index=row['witness_index'], curve_key=witness['curve_key'])

assert expected == set(witnesses), (expected-set(witnesses), set(witnesses)-expected)
counts = dict(Counter(w['method'] for w in witnesses.values()))
assert counts == dict(four_square_AP=4, repeated_cells=1, cubic=650, quartic=55, quartic_fisher=1)
report = dict(all_passed=True, denominator_max=48, total_slopes=len(expected), counts=counts,
              scope='Reduced k=q/p after D8 normalization 0<q<p; all rational finite p, no height bound.',
              limitation='Coverage audit links existing certificates; it does not recompute their rank proofs.',
              witnesses=[witnesses[k] for k in sorted(witnesses)])
(root / 'census_through_48.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k != 'witnesses'}, indent=2))
