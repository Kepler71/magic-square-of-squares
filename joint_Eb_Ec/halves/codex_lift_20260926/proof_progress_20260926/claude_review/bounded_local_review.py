from pathlib import Path
import json
import random

src = Path('s3_padic_local.py').read_text()
ns = {'__name__': 'claude_local_functions'}
exec(src.split('# ------------------------------------------------------------------ основной прогон')[0], ns)
random.seed(20260926)
primes = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 73, 89, 97, 101, 113)
out = {}
for p in primes:
    qr = {x*x % p for x in range(1, p)}
    triples = [(r, c) for r in range(1, p) for c in range(1, p)
               if all((r*r + j*c) % p in qr for j in (-1, 0, 1))]
    five = [(r, b) for r in range(1, p) for b in range(1, p)
            if all((r*r + j*b) % p in qr for j in (-2, -1, 0, 1, 2))]
    rec = {'three_AP_seed_count': len(triples), 'five_AP_seed_count': len(five),
           'type_b': {'samples': 0, 'failures': []},
           'type_c': {'samples': 0, 'failures': []},
           'type_bc': {'samples': 0, 'failures': []},
           'type_bpc': {'samples': 0, 'failures': []}}
    if triples:
        for k in range(30):
            r0, c = random.choice(triples)
            a = r0*r0
            b = p**random.choice((1, 3)) * random.randrange(1, p)
            roots = ns['make_roots'](a, b, c, p)
            assert roots is not None
            bad, _ = ns['check_type_b'](a, b, c, roots, p)
            rec['type_b']['samples'] += 1
            if bad: rec['type_b']['failures'].append(bad)
            transposed_problem_roots = ns['make_roots'](a, c, b, p)
            assert transposed_problem_roots is not None
            bad, _ = ns['check_type_b'](a, b, c, ns['transpose'](transposed_problem_roots), p)
            rec['type_c']['samples'] += 1
            if bad: rec['type_c']['failures'].append(bad)
    for k in range(30):
        a = random.randrange(1, p)**2
        b = p * random.randrange(1, p)
        c = p*p * random.randrange(1, p)
        roots = ns['make_roots'](a, b, c, p)
        assert roots is not None
        pred, _ = ns['predict_common'](a, b, c, roots, p)
        got = ns['classes'](roots, p)
        bad = [(nm, j, pred[nm][j], got[nm]['coords'][j])
               for nm in pred for j in range(3)
               if pred[nm][j] != tuple(got[nm]['coords'][j])]
        rec['type_bc']['samples'] += 1
        if bad: rec['type_bc']['failures'].append(bad)
    if five:
        for k in range(30):
            r0, b = random.choice(five)
            a = r0*r0
            c = -b + p * random.randrange(1, p)
            roots = ns['make_roots'](a, b, c, p)
            assert roots is not None
            exact, allowed, _, _, _ = ns['check_type_bpc'](a, b, c, roots, p)
            rec['type_bpc']['samples'] += 1
            if not (exact and allowed): rec['type_bpc']['failures'].append([exact, allowed])
    out[p] = rec
    print(p, {k: v if not isinstance(v, dict) else (v['samples'], len(v['failures'])) for k, v in rec.items()}, flush=True)
assert not out[5]['three_AP_seed_count']
assert not out[7]['three_AP_seed_count']
assert all(not r[k]['failures'] for r in out.values() for k in ('type_b', 'type_c', 'type_bc', 'type_bpc'))
Path('bounded_local_review.json').write_text(json.dumps(out, indent=2))
