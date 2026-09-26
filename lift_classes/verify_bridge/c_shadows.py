# (в) Какая «тень» Claude исключена новым критерием и корректно ли.
# Также: (i) исключала ли её уже лемма B Claude при p=3 (код делал исключение для 3|N);
#        (ii) что даёт тот же критерий, применённый ко ВСЕМ оценкам (числители тоже), а не только к полюсам.
from sage.all import *
import json
from pathlib import Path
OUT = Path(__file__).resolve().parent
src = Path('/home/kep/magicKube/joint_Eb_Ec/lattice2d/equal_area_rank_M1500.json')
data = json.loads(src.read_text())

def rule(X, mid, p, only_poles):
    """None если прогрессия X[e1]+X[e2]=2X[mid] не исключается простым p; иначе причина."""
    v = [x.valuation(p) for x in X]
    if only_poles:
        if max(-vi for vi in v) <= 0: return None
    mu = min(v); at = [i for i in range(3) if v[i] == mu]
    if len(at) == 3: return None
    if len(at) == 1: return 'unique min valuation'
    if p == 2: return 'two min valuations at 2'
    high = [i for i in range(3) if i not in at][0]
    need = -1 if high == mid else 2      # исчез средний -> крайние: отношение -1; исчез крайний -> 2
    return None if legendre_symbol(need, p) == 1 else 'needs (%d/%d)=1' % (need, p)

out = []
for rec in data['R1_triples']:
    if not rec['ap_shadow']: continue
    N = ZZ(rec['N']); E = EllipticCurve([-N**2, 0])
    X = []
    for m, n, s in rec['tri']:
        Z = QQ(m*m + n*n)/s; x = Z**2/4
        assert all((x + dd).is_square() for dd in (-N, 0, N))        # точка из 2E_N(Q)
        X.append(x)
    assert len(set(X)) == 3
    defects = [X[(k+1)%3] + X[(k+2)%3] - 2*X[k] for k in range(3)]
    assert all(d != 0 for d in defects)
    rootden = [ZZ(x.denominator()).sqrt() for x in X]
    primes = sorted(set(sum([prime_divisors(x.numerator()) + prime_divisors(x.denominator()) for x in X], [])))
    poles_only = {k: {int(p): rule(X, k, p, True) for p in primes if rule(X, k, p, True)} for k in range(3)}
    allval = {k: {int(p): rule(X, k, p, False) for p in primes if rule(X, k, p, False)} for k in range(3)}
    lemmaB3 = [x.valuation(3) for x in X]
    out.append(dict(N=int(N), N_mod3=int(N % 3), tri=rec['tri'], root_den=list(map(int, rootden)),
                    root_den_factored=[str(factor(d)) for d in rootden],
                    v3_of_X=list(map(int, lemmaB3)), lemmaB_p3_excludes=len(set(lemmaB3)) > 1,
                    old_code_eq3_skipped_because_3_divides_N=(N % 3 == 0),
                    survivors_poles_only=[k for k in range(3) if not poles_only[k]],
                    witnesses_poles_only={k: v for k, v in poles_only.items()},
                    survivors_all_valuations=[k for k in range(3) if not allval[k]],
                    witnesses_all_valuations={k: v for k, v in allval.items()},
                    numerators_factored=[str(factor(x.numerator())) for x in X]))
# положительный контроль: настоящая прогрессия квадратов с полюсами (1,25,49)/121 и (7^2,13^2,17^2)/p^2...
for trip in [(1, 25, 49), (49, 169, 289), (1, 841, 1681), (529, 1369, 2209)]:
    for sc in (1, 11**2, 3**2*19**2):
        Xc = [QQ(t)/sc for t in trip]
        assert Xc[0] + Xc[2] == 2*Xc[1]
        pr = sorted(set(sum([prime_divisors(x.numerator()) + prime_divisors(x.denominator()) for x in Xc], [])))
        assert all(rule(Xc, 1, p, False) is None for p in pr), (trip, sc)
res = dict(shadows=out, positive_controls='4 настоящие прогрессии квадратов x 3 масштаба проходят критерий по всем оценкам')
(OUT / 'c_shadows.json').write_text(json.dumps(res, indent=1, default=str))
for o in out:
    print(o['N'], 'N mod 3 =', o['N_mod3'], 'rootden', o['root_den'], o['root_den_factored'])
    print('   v3(X) =', o['v3_of_X'], ' lemma B (p=3) excludes:', o['lemmaB_p3_excludes'])
    print('   poles-only survivors:', o['survivors_poles_only'], o['witnesses_poles_only'])
    print('   all-valuation survivors:', o['survivors_all_valuations'], o['witnesses_all_valuations'])
    print('   numerators:', o['numerators_factored'])
