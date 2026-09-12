# Полный прогон одного сечения: модель E1, генератор, насыщенность, эллиптический Шаботи
# с перечислением по РАЗЛИЧНЫМ линиям (без двойного счёта классов).
# Использование: sage all.sage b h n [Mser] [prec] [pmax]
load('/home/kep/magicKube/family/independent_chabauty/chab.sage')
import sys, time

b0, h0, n0 = [ZZ(a) for a in sys.argv[1:4]]
Mser = int(sys.argv[4]) if len(sys.argv) > 4 else 16
prec = int(sys.argv[5]) if len(sys.argv) > 5 else 50
pmax = int(sys.argv[6]) if len(sys.argv) > 6 else 200
t00 = time.time()

A = ZZ((h0^2 + n0^2)/2); C = ZZ((b0^2 + n0^2)/2)
bet = QQ(C - A)/QQ(C + A)
# известная точка: u0 = b/n; u0^2 - 4beta' = (h/n)^2 => z = (b +- h)/(2n) рационально
u0 = QQ(b0)/n0
z0 = QQ(b0 + h0)/(2*n0)
assert z0 + bet/z0 == u0, "u0 = b/n не даёт z = (b+h)/2n"
t0 = (1 + z0)/(1 - z0)
Rt = PolynomialRing(QQ, 'T'); T_ = Rt.gen()
fT = T_*(T_^2 - 1)*(A*T_ - C)*(C*T_ - A)
D = sqfree(fT(t0).numerator()*fT(t0).denominator())
print(f"### сечение ({b0},{h0},{n0}): A={A}, C={C}, beta'={bet}")
print(f"известная точка: u0 = b/n = {u0}, z = {z0}, t = X^2 = {t0}, X = {sqrt(t0) if t0.is_square() else None}")
print(f"класс D = sqfree(f(t)) = {D}")

S = Section(b0, h0, n0, D)
E = S.E; k = S.k
G = S.point_from_u(u0)
assert G is not None, "известная точка не лежит на кривой класса D"
tors = E.torsion_points()
print(f"k = Q(sqrt({S.dk})), E1 = {E.ainvs()}")
print(f"корни 2-кручения: {S.roots};  x = {S.scale} * u")
print(f"torsion {E.torsion_subgroup().invariants()}; G: x={G[0]}, u={S.u_of(G)}, порядок бесконечен: {not G.has_finite_order()}")
hG = G.height()
print(f"h(G) = {hG}")

# ---------- насыщенность ----------
exact_ok = []
for l in [2, 3, 5, 7, 11, 13]:
    if any((G - Tt).division_points(l) for Tt in tors):
        print(f"!!! ОШИБКА ВЫПОЛНЕНИЯ: G делится на {l} по модулю кручения -- генератор не насыщен")
    else:
        exact_ok.append(l)
proved = {l: 0 for l in exact_ok}
for pr in primes(3, 800):
    for P in k.primes_above(pr):
        if E.has_bad_reduction(P): continue
        F = P.residue_field(); Er = E.reduction(P)
        rp = lambda pt: Er(0) if pt.is_zero() else Er([F(pt[0]), F(pt[1])])
        m = Er.cardinality(); Ab = Er.abelian_group()
        invs = [ZZ(d) for d in Ab.invariants()]
        for l in ZZ(m).prime_divisors():
            if l == pr or proved.get(l): continue
            if all(not all(ZZ(c) % gcd(l, d) == 0 for c, d in zip(Ab.discrete_log(rp(G) - rp(Tt)), invs)) for Tt in tors):
                proved[l] = pr
lmax = max([l for l in primes(2, 500) if all(q in proved for q in primes(2, l + 1))] or [1])
print(f"насыщенность: точно (division_points) для l in {exact_ok}; всего исключены l in {sorted(proved)}")
print(f"   => все l <= {lmax} исключены; для l > {lmax} потребовалась бы точка с h <= {RR(hG)/RR(lmax)^2:.5f}")

# ---------- выбор простых ----------
cands = []
for p in primes(5, pmax):
    if kronecker(S.dk, p) != 1: continue
    if E.discriminant().norm() % p == 0: continue
    Fp = GF(p); r = Fp(S.dk).sqrt()
    ok = True; ords = []; cards = []
    for sgn in (1, -1):
        red = lambda a: Fp(k(a).list()[0]) + sgn*r*Fp(k(a).list()[1])
        Eb = EllipticCurve(Fp, [red(a) for a in E.ainvs()])
        if Eb.cardinality() % p == 0: ok = False; break
        cards.append(Eb.cardinality())
        ords.append(Eb([red(G[0]), red(G[1])]).order())
    if not ok: continue
    N = lcm(ords)
    if gcd(N, p) != 1: continue
    cands.append((N, p, cards, ords))
cands.sort()
print(f"кандидаты (N, p): {[(N, p) for N, p, _, _ in cands[:8]]}")

# ---------- Шаботи по линиям ----------
def run(p, M, prec):
    ch = Chab(S, G, ZZ(p), prec, M).setup()
    known = {}
    LIM = 3*ch.N + 6
    for Tt in tors:
        for m in range(-LIM, LIM + 1):
            P = m*G + Tt
            uu = S.u_of(P)
            if uu is None or uu in QQ:
                kk = (str(ch.redpt(P, 1)), str(ch.redpt(P, -1)))
                nk = (str(-ch.redpt(P, 1)), str(-ch.redpt(P, -1)))
                known.setdefault(kk, []).append((m, uu, 2 if (m == 0 and kk == nk) else 1))
    groups = {}
    for Tt in tors:
        for m0 in range(ch.N):
            Q = m0*G + Tt
            groups.setdefault((str(ch.redpt(Q, 1)), str(ch.redpt(Q, -1))), []).append((Tt, m0))
    total = tot_known = 0; allok = True; rows = []
    for kk, members in sorted(groups.items()):
        bs = []
        for (Tt, m0) in members:
            co, gamma, br, phico = ch.theta(m0*G + Tt)
            n0_, mu, tailok, vals = strassmann(co, ch.v, ch.p, M)
            intok = all(ch.emb(c, sg).valuation() >= 0 for c in phico for sg in (1, -1))
            bndok = all(vals[j] >= j*ch.v - vpfact(j, ch.p) for j in range(M))
            if not (tailok and intok and bndok):
                allok = False
                print(f"!!! ОШИБКА ВЫПОЛНЕНИЯ p={p} линия {kk}: tail={tailok} int={intok} bound={bndok}")
            bs.append(n0_)
        if len(set(bs)) != 1:
            allok = False
            print(f"!!! ОШИБКА: несогласованные границы в склейке {kk}: {bs}")
        nb = max(bs); mk = sum(x[2] for x in known.get(kk, []))
        total += nb; tot_known += mk
        if nb: rows.append((kk, nb, mk, members, known.get(kk, [])))
    return ch, total, tot_known, allok, rows, len(groups), known

result = None
for (N, p, cards, ords) in cands[:6]:
    if N > 40:
        print(f"p={p}: N={N} слишком велик, пропуск"); continue
    t1 = time.time()
    ch, total, tot_known, allok, rows, nlines, known = run(p, Mser, prec)
    print(f"p={p}: |E(F_P)|={cards}, ord(Gbar)={ords}, N={N}, линий {nlines}, "
          f"нулей <= {total}, известных {tot_known}, проверки {allok}  [{time.time()-t1:.0f}s]")
    for kk, nb, mk, members, kn in rows:
        print(f"      линия {kk}: Strassmann {nb}, известных {mk} {kn}"
              f"{'' if nb == mk else '   <-- ЛИШНИХ ' + str(nb - mk)}")
    if total == tot_known and allok:
        result = (p, nlines, total, tot_known)
        print(f"*** ПОЛНЫЙ УЧЁТ при p = {p}")
        break
if result is None:
    print("*** полного учёта не получено на пробованных простых")
print(f"ИТОГ ({b0},{h0},{n0}): D={D}, k=Q(sqrt({S.dk})), tors={E.torsion_subgroup().invariants()}, "
      f"h(G)={hG}, sat<= l<= {lmax}, result={result}  [{time.time()-t00:.0f}s]")
