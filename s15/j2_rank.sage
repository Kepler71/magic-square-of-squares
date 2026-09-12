# rank J2 = Jac(C_J2), C_J2: D W^2 = f(t) = t(t^2-1)(109t-229)(229t-109)  (компонента рода 2 фактора C_{1234}).
# z = (t-1)/(t+1): t -> 1/t становится z -> -z; ветвление z = 0, oo, ±1, ±60/169.
# y^2 = c z (z^4 + A z^2 + B), B = beta^2, beta = 60/169 рационально =>
# инволюция tau: z -> beta/z, y -> gamma y / z^3, gamma^2 = beta^3, определена над k = Q(sqrt(beta)) = Q(sqrt15).
# Инварианты u = z + beta/z, v = y (1 + gamma/z^3):  v^2 = c (u^2 + A - 2beta)(u - sb)^2 (u + 2sb), sb = sqrt(beta)
# => E1: v1^2 = c (u^2 + A - 2beta)(u + 2 sb) над k;  Jac(C_J2) ~_Q Res_{k/Q} E1 => rank J2(Q) = rank E1(k).
import functools
print = functools.partial(print, flush=True)
Rt.<t> = QQ[]
f = t*(t^2 - 1)*(109*t - 229)*(229*t - 109)
Rz.<z> = QQ[]
Fz = Rz(((1 - z)^6 * f((1 + z)/(1 - z))).numerator()) if False else None
# F(z) = (1-z)^6 f((1+z)/(1-z)); W' = W (1-z)^3
Fr = Rz.fraction_field()
Fz = Rz(Fr((1 - z)^6) * Fr(f((1 + z)/(1 - z))))
print("F(z) =", Fz.factor(), "  deg", Fz.degree())
c = Fz.leading_coefficient()
q4 = Rz(Fz / (c*z))
A, B = q4[2], q4[0]
assert q4[4] == 1 and q4[3] == 0 and q4[1] == 0
beta = B.sqrt(); print("A =", A, " B =", B, " beta =", beta)
k.<r15> = QuadraticField(15)
sb = k(beta).sqrt(); print("sqrt(beta) =", sb)
Ru.<uu> = k[]
for d3t, D in [((5, 13, 65), 5*13*65*65), ((6, 26, 39), 6*26*39*39)]:
    print(f"\n===== class {d3t}: D = {D} = {D.factor()}")
    cub = Ru(D * c * (uu^2 + A - 2*beta) * (uu + 2*sb))        # (D v1)^2 = D c (...)
    a3, a2, a1, a0 = [cub[i] for i in (3, 2, 1, 0)]
    E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2]).global_minimal_model() if False else EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2])
    print("E1/k:", E1.ainvs(), "\n   j =", E1.j_invariant())
    # E1 — базовое расширение кривой над Q?  (j рационально и твист)
    print("   j in Q:", E1.j_invariant() in QQ)
    # контроль: #C_J2(F_p) = p + 1 - a_P(E1) - a_Pbar(E1) для p, расщеплённых в k
    Cq = Rt(D * f)
    ok = True; nchk = 0
    for p in primes(17, 400):
        if p in (109, 229) or kronecker(15, p) != 1 or D % p == 0: continue
        Ps = k.primes_above(p)
        if any(E1.has_bad_reduction(P) for P in Ps): continue
        Fp = GF(p)
        Cp = Cq.change_ring(Fp); cnt = sum(1 if Cp(x) == 0 else (2 if Cp(x).is_square() else 0) for x in Fp) + 1   # odd degree: 1 point at oo
        ap = sum(E1.reduction(P).trace_of_frobenius() for P in Ps)
        nchk += 1
        if cnt != p + 1 - ap: ok = False; print("   ctrl mismatch at", p)
    print(f"   control #C_J2(F_p) = p+1 - a_P - a_Pbar: {ok} on {nchk} split primes")
    try:
        lo, hi, gens = E1.simon_two_descent(lim1=5, lim3=50, limtriv=50)
        print(f"   rank E1(k) bounds (Simon 2-descent): [{lo}, {hi}]   gens found: {len(gens)}")
    except Exception as ex:
        print("   simon_two_descent failed:", ex)
    # E1 над k — если j рационально, E1 = твист базового расширения; ранг над k = rank E0 + rank E0^(15)
    if E1.j_invariant() in QQ:
        try:
            E0 = E1.descend_to(QQ)
            print("   descends to Q:", [e.ainvs() for e in E0])
            for e in E0:
                r0 = e.rank(proof=True); r15 = e.quadratic_twist(15).rank(proof=True)
                print(f"      E0 {e.ainvs()}: rank {r0}, twist by 15 rank {r15}  => rank E1(k) = {r0 + r15}")
        except Exception as ex:
            print("   descend failed:", ex)
