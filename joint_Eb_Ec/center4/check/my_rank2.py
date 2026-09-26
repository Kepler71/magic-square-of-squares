# Проверяющий: утверждение 12, двухпростая часть. Для p<q<=PMAX (p,q = 1 mod 4) и w = eps*(pi_p/pib_p)^r*(pi_q/pib_q)^s,
# |r|,|s| <= E, eps in mu_4: в каждом классе по модулю квадратов ищем связи |Im|-значений типов (1,1,1)/(1,1,2) (со знаками).
# Своя реализация: гауссовы целые как пары int, Im w = Im(rho^2)/N(rho) * (знак от eps), общий знаменатель p^E q^E.
import sys, json, time
PMAX = int(sys.argv[1]); E = int(sys.argv[2])
def isprime(n):
    if n < 2: return False
    d = 2
    while d * d <= n:
        if n % d == 0: return False
        d += 1
    return True
def gp(p):
    for x in range(1, p):
        y2 = p - x * x
        if y2 <= 0: break
        y = int(round(y2 ** 0.5))
        if y * y == y2: return (x, y)
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gpow(a, k):
    r = (1, 0)
    for _ in range(k): r = gmul(r, a)
    return r
P1 = [p for p in range(5, PMAX + 1) if p % 4 == 1 and isprime(p)]
t0 = time.time(); hits = []; nvals = 0; npairs = 0
for i, p in enumerate(P1):
    for q in P1[i+1:]:
        npairs += 1
        pi, qi = gp(p), gp(q)
        D = p ** E * q ** E
        classes = {}
        for r in range(-E, E + 1):
            A = gpow(pi if r >= 0 else (pi[0], -pi[1]), abs(r))
            for s in range(-E, E + 1):
                if r == 0 and s == 0: continue
                B = gpow(qi if s >= 0 else (qi[0], -qi[1]), abs(s))
                rho = gmul(A, B); N = rho[0]**2 + rho[1]**2   # N = p^|r| q^|s|
                sq = gmul(rho, rho)                         # rho/rhobar = rho^2/N
                scale = D // N
                for eps in range(2):   # eps=1 (Im = Im rho^2/N) и eps=i (Im(i w) = Re w); ±eps дают ±Im
                    num = abs(sq[1] if eps == 0 else sq[0]) * scale
                    if num == 0 or num == D: continue
                    cls = (eps, abs(r) % 2, abs(s) % 2)
                    classes.setdefault(cls, set()).add(num)
        for cls, S in classes.items():
            L = sorted(S); nvals += len(L)
            for a_ in range(len(L)):
                x = L[a_]
                for b_ in range(a_ + 1, len(L)):
                    y = L[b_]
                    cands = [x + y, y - x, 2*y - x, 2*x - y, x + 2*y, y + 2*x, y - 2*x]
                    if (x + y) % 2 == 0: cands.append((x + y) // 2)
                    if (y - x) % 2 == 0: cands.append((y - x) // 2)
                    for c in cands:
                        if c > 0 and c != x and c != y and c in S:
                            hits.append((p, q, cls, x, y, c))
res = dict(PMAX=PMAX, E=E, pairs=npairs, values=nvals, hits=len(hits), hit_examples=[str(h) for h in hits[:10]], time_s=round(time.time()-t0, 1))
print(res); json.dump(res, open(f'my_rank2_P{PMAX}_E{E}.json', 'w'), indent=1)
