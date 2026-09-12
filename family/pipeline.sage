# Конвейер для сечения, заданного примитивной парой b^2 + h^2 = 2 n^2 (центр 1, крест B(s)=b/n, H(s)=h/n).
#   A = (h^2+n^2)/2, C = (b^2+n^2)/2, beta = 2n^2 = A + C.
#   C_s: Y_{1,2}^2 = A(t^2+1)^2 ± beta t(t^2-1),  Y_{3,4}^2 = C(t^2+1)^2 ± beta t(t^2-1).
#   X = (Y3-Y4)/(Y1-Y2):  Q: (A X^2 - C)(X^4-1) = □,  Q': (C X^2 - A)(X^4-1) = □.
#   Дополненное покрытие класса (d1,d2,d3): p^2-q^2 = d1 r^2, p^2+q^2 = d2 s^2, A p^2 - C q^2 = d3 w^2, C p^2 - A q^2 = d3 v^2.
#   C_J2: d3 W^2 = t(t^2-1)(A t - C)(C t - A), t = X^2;  z = (t-1)/(t+1) -> корни 0, oo, ±1, ±b', b' = (C-A)/(C+A);
#   бисэллиптична над k = Q(sqrt((C-A)(C+A)))  =>  rank J2(Q) = rank E1(k).
# Этапы: (1) ELS-классы дополненной системы (строго), (2) точки малой высоты, (3) E1 над k, ранговые границы,
#        (4) подъёмы найденных X и проверка вырожденности.
import functools, itertools, sys, time
print = functools.partial(print, flush=True)
b0, h0, n0 = [ZZ(x) for x in sys.argv[1:4]]
assert b0^2 + h0^2 == 2*n0^2 and gcd([b0, h0, n0]) == 1
A, C = ZZ((h0^2 + n0^2)/2), ZZ((b0^2 + n0^2)/2)
# делить на gcd(A,C) НЕЛЬЗЯ: это закрутило бы условия квадратности; требуем gcd = 1
assert gcd(A, C) == 1, ('gcd(A,C) != 1', A, C)
print(f"section (b,h,n) = ({b0},{h0},{n0}):  A = {A}, C = {C}  (C^2-A^2 = {(C^2-A^2).factor()})")
Fh = [lambda p, q: p^2 - q^2, lambda p, q: p^2 + q^2, lambda p, q: A*p^2 - C*q^2, lambda p, q: C*p^2 - A*q^2]

def sqf(n):
    n = ZZ(n); return sign(n) * prod(pr^(e % 2) for pr, e in n.abs().factor())
def is_sq_Ql(a, l):
    a = QQ(a); v = a.valuation(l)
    if v % 2: return False
    u = a / l^v
    if l == 2: return (u.numerator() * u.denominator()) % 8 == 1
    return kronecker(ZZ(u.numerator() * u.denominator()), l) == 1
def locally_solvable(ds, l, K=None):
    # строгий тест (как verify_s15/classes.sage), 4 формы
    e = 3 if l == 2 else 1
    if K is None: K = 16 if l == 2 else 10
    unknown = []
    def rec(chart, a, k):
        p, q = (a, 1) if chart == 0 else (1, a)
        und = False
        for fh, d in zip(Fh, ds):
            v = fh(p, q)
            if v != 0 and v.valuation(l) + e <= k:
                if not is_sq_Ql(v / d, l): return False
            else:
                und = True
        if not und: return True
        if k >= K: unknown.append((chart, a, k)); return False
        return any(rec(chart, a + j*l^k, k + 1) for j in range(l))
    if rec(0, 0, 0) or rec(1, 0, 1): return True
    return False if not unknown else None
def real_ok(ds):
    # нужна вещественная точка: знаки всех четырёх форм в одной точке X
    xs = [QQ(x)/8 for x in range(-80, 81)] + [QQ(10^6)]
    return any(all(fh(x, 1) == 0 or sign(fh(x, 1)) == sign(d) for fh, d in zip(Fh, ds)) for x in xs)

# (1) кандидаты: простые, делящие попарные результанты форм (gcd(p,q)=1)
res_primes = set((2*(C - A)*(C + A)*A*C).prime_factors())
cand_primes = sorted(res_primes)
units = [s_ * prod(pr for pr, e_ in zip(cand_primes, es) if e_) for s_ in (1, -1) for es in itertools.product((0, 1), repeat=len(cand_primes))]
D2 = [d for d in units if d > 0]
check_primes = sorted(set(cand_primes) | set(primes(3, 60)))
t0 = time.time()
els = []
for d1 in units:
    for d2 in D2:
        d3 = sqf(d1*d2); ds = (d1, d2, d3, d3)
        if not real_ok(ds): continue
        bad = None
        for l in check_primes:
            r = locally_solvable(ds, l)
            if r is False: bad = l; break
            if r is None: print(f"   UNKNOWN local test {ds[:3]} at {l}")
        if bad is None: els.append(ds[:3])
print(f"(1) ELS classes of the augmented system ({len(units)*len(D2)} candidates, {time.time()-t0:.0f}s): {els}")

# (2) точки малой высоты в ELS-классах
H = int(sys.argv[4]) if len(sys.argv) > 4 else 300
pts = {cls: [] for cls in els}
for q in range(1, H):
    for p in range(0, H):
        if gcd(p, q) != 1: continue
        vals = [fh(p, q) for fh in Fh]
        if 0 in vals: continue
        cls = (sqf(vals[0]), sqf(vals[1]), sqf(vals[2]))
        if cls in pts and sqf(vals[3]) == cls[2] and all(QQ(v / d).is_square() for v, d in zip(vals, cls + (cls[2],))):
            pts[cls].append(QQ(p)/q)
print(f"(2) points with max(p,q) < {H}:", {c_: v for c_, v in pts.items()})

# (3) E1 над k для каждого ELS-класса
Rt = PolynomialRing(QQ, 'T'); T = Rt.gen()
f = T*(T^2 - 1)*(A*T - C)*(C*T - A)
Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); Frz = Rz.fraction_field()
Fz = Rz(Frz((1 - z)^6) * Frz(f((1 + z)/(1 - z))))
cz = Fz.leading_coefficient(); q4 = Rz(Fz/(cz*z))
Ap, Bp = q4[2], q4[0]
bet = (C - A)/(C + A)
assert Bp == bet^2 and Ap == -(1 + bet^2), (Ap, Bp, bet)
dk = sqf((C - A)*(C + A))
k = QuadraticField(dk, 'rk'); rk = k.gen()
sb = k(bet).sqrt()
Ru = PolynomialRing(k, 'uu'); uu = Ru.gen()
print(f"(3) k = Q(sqrt({dk})), beta' = {bet}, sqrt(beta') = {sb}")
E1s = {}
for cls in els:
    D = cls[2]
    cub = Ru(D * cz * (uu^2 + Ap - 2*bet) * (uu + 2*sb))
    a3, a2, a1, a0 = [cub[j] for j in (3, 2, 1, 0)]
    E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2])
    Em = E1.global_minimal_model(semi_global=True)
    lo, hi = None, None
    if 'simon' in sys.argv:
        try:
            lo, hi, gens = Em.simon_two_descent()
        except Exception as ex:
            pass
    E1s[cls] = (E1, Em)
    print(f"   class {cls}: E1 = {Em.ainvs()}\n      j = {Em.j_invariant()}, torsion {Em.torsion_subgroup().invariants()}, Simon rank bounds [{lo},{hi}]")
    # контроль расщепления: #C_J2(F_p) = p+1 - a_P - a_Pbar на расщеплённых p
    Cq = Rt(D * f); okc = True; nn = 0
    for p in primes(11, 300):
        if kronecker(dk, p) != 1 or (2*D*A*C*(C^2 - A^2)) % p == 0: continue
        Ps = k.primes_above(p)
        if any(E1.has_bad_reduction(P) for P in Ps): continue
        Cp = Cq.change_ring(GF(p))
        cnt = sum(1 if Cp(x) == 0 else (2 if Cp(x).is_square() else 0) for x in GF(p)) + 1
        nn += 1
        if cnt != p + 1 - sum(E1.reduction(P).trace_of_frobenius() for P in Ps): okc = False
    print(f"      control #C_J2(F_p) = p+1-a_P-a_Pbar on {nn} split primes: {okc}")

# (4) подъёмы найденных X: все рациональные t из тождества (1) и число различных клеток
Bf = lambda u_: (1 + 2*u_ - u_^2)/(1 + u_^2); Hf = lambda u_: (1 - 2*u_ - u_^2)/(1 + u_^2)
bS, hS = QQ(b0)/n0, QQ(h0)/n0
def cells(t0):
    return [QQ(1), bS^2, hS^2, Bf(t0)^2, Hf(t0)^2, (hS^2 + Hf(t0)^2)/2, (hS^2 + Bf(t0)^2)/2, (bS^2 + Hf(t0)^2)/2, (bS^2 + Bf(t0)^2)/2]
aT = A*(T^2 + 1)^2; bT = (A + C)*T*(T^2 - 1); cT = C*(T^2 + 1)^2
g4 = lambda x: (C*x^2 - A)*(A*x^2 - C)
for cls, Xs in pts.items():
    for X0 in Xs:
        G = g4(X0)
        if not G.is_square(): print(f"(4) {cls} X={X0}: g4 not a square -> no lift"); continue
        h0_ = X0*G.sqrt()*2/((A + C)*(1 - X0^4))
        ts = set()
        for sg in (1, -1):
            ts |= set((T*(T^2 - 1) - sg*h0_*(T^2 + 1)^2).roots(QQ, multiplicities=False))
        lifts = []
        for t0 in ts:
            vals = [aT(t0) + bT(t0), aT(t0) - bT(t0), cT(t0) + bT(t0), cT(t0) - bT(t0)]
            if all(v.is_square() for v in vals):
                lifts.append((t0, len(set(cells(t0)))))
        print(f"(4) {cls} X={X0}: lifts (t, #distinct cells) = {sorted(lifts)}")
