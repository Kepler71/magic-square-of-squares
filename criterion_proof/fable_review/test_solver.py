"""Sanity controls for fr_selmer local solver (Fable review, 2026-09-13)."""
import random, math, sys, time
sys.path.insert(0, '/home/kep/magicKube/criterion_proof/fable_review')
from fr_selmer import *

random.seed(20260913)


def naive_soluble(a, b, c, p, N):
    """Brute force over z = x/p^k (k=0..3), x mod p^N; positive if f(z)p^{4k} is a square
    with the decision stable within precision.  Also returns 'maybe' if only undecidable classes."""
    maybe = False
    if a == 0:
        return True
    # infinity
    v = vp(a, p)
    if v % 2 == 0 and is_unit_square(a // p**v, p):
        return True
    for k in range(0, 4):
        pk = p ** k
        for x in range(p ** N):
            if k > 0 and x % p == 0:
                continue
            F = a * x**4 + b * x * x * pk * pk + c * pk**4
            if F == 0:
                return True
            v = vp(F, p)
            margin = 3 if p == 2 else 1
            if v + margin > N:      # F mod p^N doesn't determine unit part
                maybe = True
                continue
            if v % 2 == 0 and is_unit_square(F // p**v, p):
                return True
    return 'maybe' if maybe else False


t0 = time.time()
# --- 1. random quartics vs naive search
mism = 0; tested = 0; maybes = 0
for p, N in [(2, 11), (3, 6), (5, 5), (7, 4)]:
    for _ in range(300):
        a = random.choice([1, -1, 2, -2, 3, 5, 6, -3, 7, 10, -6, 15, 21, -15, 30]) * random.choice([1, 1, 1, p, p*p])
        b = random.randint(-60, 60)
        c = random.choice([1, -1, 2, -2, 3, 5, 6, -3, 7, 10, -6, 15, 21, 4, 9, 25, -5, 12]) * random.choice([1, 1, p, p*p, p**3])
        # discriminant of a u^2 + b u + c nonzero and a,c nonzero to keep squarefree-ish
        if a == 0 or c == 0 or b*b - 4*a*c == 0:
            continue
        mine = quartic_locally_soluble(a, b, c, p)
        nv = naive_soluble(a, b, c, p, N)
        tested += 1
        if nv == 'maybe':
            maybes += 1
            if mine is False:
                pass   # naive undecided, mine says no: cannot compare
            continue
        if mine != nv:
            mism += 1
            print('MISMATCH', p, (a, b, c), 'mine', mine, 'naive', nv)
print(f'[1] random quartics vs naive: tested {tested}, mismatches {mism}, naive-undecided {maybes}, {time.time()-t0:.1f}s')

# --- 2. positive controls on pairs
t0 = time.time()
bad = 0; cnt = 0
for _ in range(400):
    m = random.randint(1, 400); n = random.randint(1, 400)
    if m == n or math.gcd(m, n) != 1:
        continue
    S1, S2 = selmer_sets(m, n)
    cnt += 1
    for d in (1, -1, m*n, -m*n):
        # squarefree part of d
        if d not in S1 and d != 0:
            # d may be non-squarefree: reduce
            pass
    def sqf(d):
        s = 1 if d > 0 else -1; d = abs(d); out = 1
        for q in prime_factors(d):
            if vp(d, q) % 2 == 1:
                out *= q
        return s * out
    for d in (1, -1, m*n, -m*n):
        if sqf(d) not in S1:
            bad += 1; print('POSITIVE CONTROL FAIL Sel_phi', m, n, d, S1)
    if 1 not in S2:
        bad += 1; print('POSITIVE CONTROL FAIL Sel_phihat', m, n, S2)
    if len(S1) & (len(S1) - 1) or len(S2) & (len(S2) - 1):
        bad += 1; print('NOT POWER OF 2', m, n, S1, S2)
    # subgroup check (Selmer sets are groups mod squares)
    def sqfmul(x, y):
        return sqf(x * y)
    for x in S1:
        for y in S1:
            if sqfmul(x, y) not in S1:
                bad += 1; print('NOT A GROUP Sel_phi', m, n, x, y); break
    for x in S2:
        for y in S2:
            if sqfmul(x, y) not in S2:
                bad += 1; print('NOT A GROUP Sel_phihat', m, n, x, y); break
print(f'[2] positive/group controls on {cnt} pairs: failures {bad}, undetermined {UNDETERMINED[0]}, {time.time()-t0:.1f}s')

# --- 3. primes outside 2mn(m^2-n^2) never obstruct (sample)
t0 = time.time()
viol = 0; checks = 0
for _ in range(60):
    m = random.randint(1, 60); n = random.randint(1, 60)
    if m == n or math.gcd(m, n) != 1:
        continue
    S1, S2 = selmer_sets(m, n)
    extra = [p for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] if (2 * m * n * (m*m - n*n)) % p != 0]
    T1, T2 = selmer_sets(m, n, extra_primes=extra)
    checks += 1
    if S1 != T1 or S2 != T2:
        viol += 1; print('EXTRA PRIME OBSTRUCTS', m, n, S1, T1, S2, T2)
print(f'[3] good primes never obstruct: {checks} pairs, violations {viol}, {time.time()-t0:.1f}s')
