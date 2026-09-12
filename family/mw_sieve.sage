# Решето Морделла–Вейля для «необъяснённых» p-адических нулей (сечение (b,h,n), класс D, точка t0).
# Неизвестная P с x(P) in Q: образ P в конечных группах E(F_q1) x E(F_q2) равен образу T + nG (T in E[2], n in Z),
# если [Gamma : <G, E[2]>] взаимно прост с порядками этих групп (проверяем l-насыщение G для всех l | порядков).
# Условия: (a) при Chabauty-простых — P в классе с необъяснённым нулём; (b) при вспомогательных q — x(P) mod q1 == x(P) mod q2.
import functools, sys, itertools
print = functools.partial(print, flush=True)
args = sys.argv[1:]
b0, h0, n0, D = [ZZ(x) for x in args[:4]]; t0_arg = QQ(args[4])
# формат Chabauty-условий: "p:N:T#,r;T#,r"  (классы с необъяснёнными нулями)
chab = []
for spec in args[5:]:
    if spec.startswith('aux='): continue
    p_, N_, cl = spec.split(':'); chab.append((int(p_), int(N_), [tuple(int(v) for v in c.split(',')) for c in cl.split(';')]))
aux_max = int([a for a in args if a.startswith('aux=')][0][4:]) if any(a.startswith('aux=') for a in args) else 400
_saved_argv = list(sys.argv)
sys.argv = [sys.argv[0]] + [str(x) for x in (b0, h0, n0, D, t0_arg)] + [str(chab[0][0])]
exec(preparse(open('chabauty.sage').read().split('# --- p-адические вложения')[0].replace(
    'assert kronecker(dk, p) == 1 and all(not E1.has_bad_reduction(P) for P in k.primes_above(p)), "p must split and be good"', 'pass')))
def red_data(q):
    Ps = k.primes_above(q)
    res = []
    for P in Ps:
        F = k.residue_field(P); Er = E1.reduction(P)
        def red1(Q_, F=F, Er=Er, P=P):
            if Q_.is_zero() or Q_[0].valuation(P) < 0: return Er(0)
            return Er([F(Q_[0]), F(Q_[1])])
        Gr = red1(G); Tr = [red1(T_) for T_ in T2tors]
        # точки T + rG считаем в E(F_q): red(T) + r*red(G)
        res.append((F, Er, red1, Gr, Tr))
    return res
def order_G(q):
    return lcm([rd[3].order() for rd in red_data(q)])
sys.argv = _saved_argv
# (a) Chabauty-условия как множества (T#, n mod N)
allowed = None; M = 1
for (p_, N_, cl) in chab:
    rds = red_data(p_)
    targets = [tuple(rd[4][ti] + r*rd[3] for rd in rds) for ti, r in cl]
    S = set()
    for ti in range(4):
        for r in range(N_):
            img = tuple(rd[4][ti] + r*rd[3] for rd in rds)
            if img in targets: S.add((ti, r))
    print(f"Chabauty prime {p_} (N={N_}): allowed (T#, n mod N) for an unknown point: {sorted(S)}")
    if allowed is None: allowed, M = S, N_
    else:
        M2 = lcm(M, N_); new = set()
        for (ti, r) in allowed:
            for (tj, s) in S:
                if ti != tj: continue
                try: new.add((ti, crt(r, s, M, N_) % M2))
                except ValueError: pass
        allowed, M = new, M2
    print(f"   after combining: modulus {M}, {len(allowed)} residue classes")
# (b) вспомогательные q: x(T + nG) mod q1 == mod q2 (как элементы F_q; бесконечность == бесконечность)
used_l = set()
for q in primes(7, aux_max):
    if not allowed: break
    if kronecker(dk, q) != 1 or any(E1.has_bad_reduction(P) for P in k.primes_above(q)) or (2*D*A*C) % q == 0: continue
    rds = red_data(q)
    Nq = lcm([rd[3].order() for rd in rds])
    grp_orders = [rd[1].order() for rd in rds]
    if max(Nq, 1) > 2000: continue
    ok_set = set()
    for ti in range(4):
        for r in range(Nq):
            pts = [rd[4][ti] + r*rd[3] for rd in rds]
            xs = [None if P_.is_zero() else ZZ(P_[0]) for P_ in pts]    # F_q -> {0..q-1}
            if xs[0] == xs[1]: ok_set.add((ti, r))
    M2 = lcm(M, Nq); new = set()
    for (ti, r) in allowed:
        for s in range(Nq):
            if (ti, s) in ok_set:
                try: new.add((ti, crt(r, s, M, Nq) % M2))
                except ValueError: pass
    if len(new) < len(allowed) * (M2 // M):
        for go in grp_orders: used_l |= set(ZZ(go).prime_factors())
        allowed, M = new, M2
        print(f"aux q={q}: N_q={Nq}, modulus {M}, {len(allowed)} classes left")
    if M > 10^12: break
used_l |= set(ZZ(c[0]).prime_factors() for c in []) if False else set()
for (p_, N_, cl) in chab:
    for rd in red_data(p_): used_l |= set(ZZ(rd[1].order()).prime_factors())
    used_l.add(p_)
sat = {l: E1.saturation([G], one_prime=ZZ(l))[1] for l in sorted(used_l)}
print("l-saturation of G for all l dividing used group orders:", sat)
print("RESULT:", "no unknown point (sieve empty)" if not allowed else f"{len(allowed)} classes survive mod {M}")
