# Независимая точная локальная проверка шагов 2-5 теоремы об остаточных классах.
# Только целые числа; критерий квадратности в Q_p: чётная оценка и
#   p=2: единичная часть ≡1 mod 8;  p нечётно: символ Лежандра единичной части = 1.
# Выборка: случайные (r,s) из семейств с вынужденной делимостью p^k | r, s, r-s, r+s;
#   z = p^n u/v, u,v случайны и взаимно просты с p, |u|,|v| до 10^6, n в [-7,7].
# Все циклы с явной верхней границей NS. Прогресс печатается.
import sys, random, json, time
from math import gcd

def vpu(n, p):
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v, n

def issq(n, p):          # n целое ненулевое: квадрат ли в Q_p
    v, u = vpu(n, p)
    if v & 1: return False
    if p == 2: return u % 8 == 1
    return pow(u % p, (p - 1) // 2, p) == 1

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
LINES = [[(i, j) for j in (-1, 0, 1)] for i in (-1, 0, 1)] + \
        [[(i, j) for i in (-1, 0, 1)] for j in (-1, 0, 1)] + \
        [[(t, t) for t in (-1, 0, 1)], [(t, -t) for t in (-1, 0, 1)]]

def rand_unit(p, B, rng):
    while True:
        x = rng.randint(1, B)
        if x % p: return x * rng.choice((1, -1))

def rand_rs(p, rng):
    fam = rng.randrange(6)
    for _ in range(1000):
        k = rng.randint(1, 4)
        m = rng.randint(1, 3000) * rng.choice((1, -1))
        s = rng.randint(1, 3000) * rng.choice((1, -1))
        if fam == 0:
            r = rng.randint(1, 3000) * rng.choice((1, -1))
        elif fam == 1: r = s + p**k * m          # p^k | r-s
        elif fam == 2: r = -s + p**k * m         # p^k | r+s
        elif fam == 3: r = p**k * m              # p^k | r
        elif fam == 4: r, s = s, p**k * m        # p^k | s
        else:           # 2-адически: оба нечётные / разной чётности
            r = 2 * m + 1; s = 2 * s + rng.choice((0, 1))
        if r and s and r != s and r != -s and gcd(r, s) == 1:
            return r, s, fam
    raise RuntimeError('no rs')

def run(p, NS, seed):
    rng = random.Random(seed)
    st = dict(tested=0, passed=0, violations=[], T_odd=0, L_odd=0,
              pole_pass=0, integral_pass=0, fam_pass=[0]*6)
    t0 = time.time()
    for it in range(NS):
        r, s, fam = rand_rs(p, rng)
        n = rng.randint(-7, 3)
        u = rand_unit(p, 10**6, rng); v = abs(rand_unit(p, 10**6, rng))
        num = u * p**max(n, 0); den = v * p**max(-n, 0)
        N = {(i, j): den + (i*r + j*s) * num for (i, j) in CELLS}
        if any(x == 0 for x in N.values()): continue
        st['tested'] += 1
        ok = all(issq(N[a]*N[b]*N[c]*den, p) for (a, b, c) in LINES)
        if not ok: continue
        st['passed'] += 1; st['fam_pass'][fam] += 1
        vz = vpu(num, p)[0] - vpu(den, p)[0]
        if vz < 0: st['pole_pass'] += 1
        else: st['integral_pass'] += 1
        cv = {k: (vpu(N[k], p)[0] - vpu(den, p)[0]) for k in CELLS}
        T = N[(-1, 0)] * N[(1, 1)] * N[(0, -1)] * den
        L = N[(-1, 0)] * N[(1, -1)] * N[(0, 1)] * den
        vT = vpu(T, p)[0] & 1; vL = vpu(L, p)[0] & 1
        S, U = cv[(-1, -1)] & 1, cv[(-1, 0)] & 1
        bad = []
        # матрица классов по оценкам
        exp = {(-1,-1): S, (1,1): S, (-1,1): S ^ U, (1,-1): S ^ U,
               (-1,0): U, (1,0): U, (0,-1): U, (0,1): U, (0,0): 0}
        if any((cv[k] & 1) != exp[k] for k in CELLS): bad.append('class-matrix parity')
        if vT != S or vL != (S ^ U): bad.append('[T]=S or [L]=S+U parity')
        if p in (2, 3):
            if not all(issq(N[k]*den, p) for k in CELLS): bad.append('not all nine squares')
            vb = vpu(r, p)[0] + vz; vc = vpu(s, p)[0] + vz
            if min(vb, vc) < (3 if p == 2 else 1): bad.append('v(b),v(c) too small')
            if not (issq(T, p) and issq(L, p)): bad.append('T or L not square')
        else:
            if vT and (r - s) % p: bad.append('vT odd but p !| r-s')
            if vL and (r + s) % p: bad.append('vL odd but p !| r+s')
            if vz >= 0 and any(cv[k] & 1 for k in CELLS): bad.append('integral z but odd val')
            if p % 4 == 3:
                if vT or vL: bad.append('p=3mod4 odd val of T/L')
                if vz < 0: bad.append('p=3mod4 but pole (derived claim)')
            if vz < 0:
                nn = -vz
                if r % p == 0 or s % p == 0:
                    if nn & 1 or S or U: bad.append('case p|r or p|s')
                elif nn % 2 == 0:
                    if S or U: bad.append('units, n even')
                else:
                    if U != 1: bad.append('n odd but U=0')
                    if (r - s) % p == 0 and not (S == 1 and (S ^ U) == 0): bad.append('p|r-s case')
                    if (r + s) % p == 0 and not (S == 0 and (S ^ U) == 1): bad.append('p|r+s case')
                    if (r - s) % p and (r + s) % p: bad.append('n odd, both r±s units, passed')
        st['T_odd'] += vT; st['L_odd'] += vL
        if bad and len(st['violations']) < 20:
            st['violations'].append(dict(p=p, r=r, s=s, num=num, den=den, bad=bad))
        if bad: st.setdefault('nviol', 0); st['nviol'] = st.get('nviol', 0) + 1
        if (it + 1) % (NS // 4) == 0:
            print(f'  p={p} {it+1}/{NS} tested={st["tested"]} passed={st["passed"]} '
                  f'Todd={st["T_odd"]} Lodd={st["L_odd"]} viol={st.get("nviol",0)} '
                  f'{time.time()-t0:.0f}s', flush=True)
    return st

if __name__ == '__main__':
    primes = [int(x) for x in sys.argv[1].split(',')]
    NS = int(sys.argv[2]); seed = int(sys.argv[3]); out = sys.argv[4]
    res = {}
    for p in primes:
        res[p] = run(p, NS, seed * 1000 + p)
        print(p, {k: v for k, v in res[p].items() if k != 'violations'},
              'violations:', res[p]['violations'][:3], flush=True)
    json.dump(res, open(out, 'w'), indent=1)
    tot = sum(r.get('nviol', 0) for r in res.values())
    print('TOTAL VIOLATIONS:', tot)
