# Верхняя граница rank J2 без Magma: 15 Richelot-изогений C_J2 (все корни рациональны) -> ищем двойственную кривую,
# бисэллиптичную НАД Q (подъём инволюции рационален) -> Jac ~_Q E_a x E_b -> rank = rank E_a + rank E_b (mwrank, proof).
# C_J2: D W^2 = t(t^2-1)(A t - C)(C t - A); корни 0, oo, 1, -1, C/A, A/C.
import functools, itertools, sys
print = functools.partial(print, flush=True)
b0, h0, n0, D = [ZZ(x) for x in sys.argv[1:5]]
A, C = ZZ((h0^2 + n0^2)/2), ZZ((b0^2 + n0^2)/2)
R = PolynomialRing(QQ, 't'); t = R.gen()
f0 = D*t*(t^2 - 1)*(A*t - C)*(C*t - A)
roots = [QQ(0), oo, QQ(1), QQ(-1), QQ(C)/A, QQ(A)/C]
def quad(r1, r2):
    if r1 is oo: r1, r2 = r2, r1
    return (t - r1) if r2 is oo else (t - r1)*(t - r2)
def matchings(lst):
    if not lst: yield []; return
    a = lst[0]
    for j in range(1, len(lst)):
        for m in matchings(lst[1:j] + lst[j+1:]): yield [(a, lst[j])] + m
def jac_quartic(g, F_):
    a_, b_, c_, d_, e_ = [g[j] for j in (4, 3, 2, 1, 0)]
    I = 12*a_*e_ - 3*b_*d_ + c_^2; J = 72*a_*c_*e_ + 9*b_*c_*d_ - 27*a_*d_^2 - 27*e_*b_^2 - 2*c_^3
    return EllipticCurve(F_, [-27*I, -27*J])
def split_over_Q(f1):
    # инволюции P^1 над Q, переставляющие корни f1 парами; нужен рациональный подъём (kappa — квадрат)
    L = f1.splitting_field('a') if f1.degree() > 0 else QQ
    RL = PolynomialRing(L, 'x'); x = RL.gen()
    rts = [r for r, _ in f1.change_ring(L).roots()]
    if f1.degree() == 5: rts = rts + [oo]
    out = []
    for m in matchings(rts):
        pairs = m
        fin = [(a_, b_) for a_, b_ in pairs if a_ is not oo and b_ is not oo]
        # M(z) = (al z + be)/(ga z - al); M(a)=b: al(a+b) + be - ga ab = 0 ; пара с oo: M(a)=oo => ga a - al = 0
        rows = []
        for a_, b_ in pairs:
            if a_ is oo or b_ is oo:
                r_ = b_ if a_ is oo else a_
                rows.append([-1, 0, r_])
            else:
                rows.append([a_ + b_, 1, -a_*b_])
        Mt = matrix(L, rows)
        ker = Mt.right_kernel().basis()
        if len(ker) != 1: continue
        al, be, ga = ker[0]
        if al^2 + be*ga == 0: continue
        v = vector(L, [al, be, ga]); piv = next(c_ for c_ in v if c_ != 0); v = v/piv
        if not all(c_ in QQ for c_ in v): continue
        al, be, ga = [QQ(c_) for c_ in v]
        Fr = R.fraction_field()
        Mz = (al*t + be)/(ga*t - al)
        kap = Fr(f1(Mz) * (ga*t - al)^6) / Fr(f1)
        if kap not in QQ: continue
        out.append(((al, be, ga), QQ(kap)))
    return out
print(f"A={A}, C={C}, D={D}")
results = []
for idx, m in enumerate(matchings(roots)):
    F = [quad(a_, b_) for a_, b_ in m]
    lead = f0 // prod(F)
    assert lead * prod(F) == f0 and lead.degree() == 0
    F[0] = lead * F[0]
    co = lambda P: [P[0], P[1], P[2]]
    delta = matrix(QQ, [co(F[0]), co(F[1]), co(F[2])]).det()
    if delta == 0: continue
    G1 = F[1].derivative()*F[2] - F[1]*F[2].derivative()
    G2 = F[2].derivative()*F[0] - F[2]*F[0].derivative()
    G3 = F[0].derivative()*F[1] - F[0]*F[1].derivative()
    f1 = delta*G1*G2*G3
    # контроль изогении на нескольких простых
    bad = set((2*f0.discriminant()*f1.discriminant()).numerator().prime_factors()) | set(QQ(f1.leading_coefficient()).numerator().prime_factors())
    okF = all(HyperellipticCurve(f0.change_ring(GF(p))).frobenius_polynomial() == HyperellipticCurve(f1.change_ring(GF(p))).frobenius_polynomial()
              for p in primes(5, 80) if p not in bad and QQ(f1.leading_coefficient()).valuation(p) == 0)
    sp = split_over_Q(f1)
    rat = [s_ for s_ in sp if s_[1].is_square()]
    results.append((idx, m, okF, sp, rat, f1))
    print(f"matching {idx:2d}: {[(str(a_), str(b_)) for a_, b_ in m]}  isogeny ctrl {okF};  Q-involutions: {len(sp)}, with rational lift: {len(rat)}")
    for (al, be, ga), kap in rat:
        # квотиенты: переносим инволюцию в z -> -z через неподвижные точки (над K), затем descend_to(QQ)
        fixpol = R(ga*t^2 - 2*al*t - be)
        K = fixpol.splitting_field('c') if not fixpol.is_constant() and len(fixpol.roots(QQ)) < 2 else QQ
        RK = PolynomialRing(K, 'w'); w = RK.gen()
        if ga != 0:
            fps = [r_ for r_, _ in fixpol.change_ring(K).roots()]
            p1, p2 = fps
        else:
            p1, p2 = K(-be/(2*al)), None
        FK = RK.fraction_field()
        zz = (p1 + p2*w)/(1 + w) if p2 is not None else p1 + w
        num = FK(f1.change_ring(K)(zz)) * FK((1 + w)^6 if p2 is not None else 1)
        Fw = RK(num.numerator()) / RK(num.denominator())
        Fw = RK(Fw)
        ev = all(Fw[j] == 0 for j in range(1, Fw.degree() + 1, 2))
        od = all(Fw[j] == 0 for j in range(0, Fw.degree() + 1, 2))
        if ev:
            Gt = RK([Fw[2*j] for j in range(Fw.degree()//2 + 1)])
            Ea = jac_quartic(Gt if Gt.degree() == 4 else Gt, K) if Gt.degree() >= 3 else None
            Eb = jac_quartic(RK(w*Gt), K)
            for nm, E_ in (("E_a", Ea), ("E_b", Eb)):
                if E_ is None: continue
                ds = E_.descend_to(QQ)
                for e in ds[:1]:
                    em = e.minimal_model()
                    try: rr = em.rank(proof=True)
                    except Exception: rr = f"bounds {em.rank_bounds()}"
                    print(f"      {nm}: {em.ainvs()}  cond {em.conductor().factor()}  rank {rr}")
        else:
            print("      (w-model not even; skipped)", "odd" if od else "")
