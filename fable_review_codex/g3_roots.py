# Fable, 14.09.2026. Символьная проверка множеств вырождения в A, B (THEOREM_INFINITY_CONSTRAINTS) и в
# теореме о фиксированных шести (PROOF_FIXED_SIX_COMPLETE): при каких значениях параметра нули линейных форм совпадают.
import json
from sage.all import *
out = {}
p = var('p'); k = var('k'); a = var('a'); b = var('b')

def coincidences(roots, v):
    bad = set()
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            eq = SR(roots[i]) - SR(roots[j])
            if eq.is_zero():
                bad.add('тождественно')
                continue
            for s in solve(eq == 0, v):
                bad.add(str(s.rhs()))
    return sorted(bad)

# A: четыре клетки 1+q, 1-q, 1+p+q, 1+p-q -> нули по q
A4 = coincidences([-1, 1, -1 - p, 1 + p], p)
# A: все шесть меняющихся клеток 1±q, 1+p±q, 1-p±q -> нули ±1, ±1-p, ±1+p
A6 = coincidences([-1, 1, -1 - p, 1 + p, -1 + p, 1 - p], p)
# B: восемь клеток 1+lambda p, lambda in {±1, ±k, ±(1+k), ±(1-k)}: нули -1/lambda; совпадение <=> lambda_i = lambda_j; нуль lambda = 0
lams = [1, -1, k, -k, 1 + k, -1 - k, 1 - k, -1 + k]
B8 = coincidences(lams, k)
Bzero = sorted(set(str(s.rhs()) for l in lams if not SR(l).is_constant() for s in solve(SR(l) == 0, k)))
# фиксированные шесть: нули 0, b, -b, -2a
F4 = coincidences([0, b, -b, -2 * a], a) + coincidences([0, b, -b, -2 * a], b)
out = {'A_four_cells_bad_p0': A4, 'A_six_cells_bad_p0': A6, 'B_coincidence_k': B8, 'B_zero_lambda_k': Bzero,
       'fixed_six_bad': F4}
print(out)
# положительность: max(|b+c|, |b-c|) = |b|+|c|
import itertools, random
ok = True
for _ in range(2000):
    bb = QQ(random.randint(-50, 50)) / random.randint(1, 9); cc = QQ(random.randint(-50, 50)) / random.randint(1, 9)
    ok &= max(abs(bb + cc), abs(bb - cc)) == abs(bb) + abs(cc)
out['max_abs_identity_random'] = bool(ok)
# ромб |p|+|q|<1 <=> все восемь клеток 1±p, 1±q, 1±(p+q), 1±(p-q) положительны
ok2 = True
for _ in range(5000):
    pp = QQ(random.randint(-30, 30)) / 20; qq = QQ(random.randint(-30, 30)) / 20
    pos = all(1 + s * v > 0 for v in [pp, qq, pp + qq, pp - qq] for s in [1, -1])
    ok2 &= (pos == (abs(pp) + abs(qq) < 1))
out['rhombus_iff_positive_random'] = bool(ok2)
print(out['max_abs_identity_random'], out['rhombus_iff_positive_random'])
json.dump(out, open('/home/kep/magicKube/fable_review_codex/g3_roots.json', 'w'), indent=1, ensure_ascii=False)
print('OK')
