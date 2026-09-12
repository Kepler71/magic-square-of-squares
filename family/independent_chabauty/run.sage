load('/home/kep/magicKube/family/independent_chabauty/chab.sage')
import sys
b, h, n, D, u0, p = sys.argv[1:7]
prec = int(sys.argv[7]) if len(sys.argv) > 7 else 80
M = int(sys.argv[8]) if len(sys.argv) > 8 else 30
S = Section(ZZ(b), ZZ(h), ZZ(n), ZZ(D))
E = S.E
G = S.point_from_u(QQ(u0))
print(f"section ({b},{h},{n}) D={D}: E = {E.ainvs()},  k = Q(sqrt({S.dk}))")
print(f"G = {G[0]}, u(G) = {S.u_of(G)}, h(G) = {G.height()}")
ch = Chab(S, G, ZZ(p), prec, M).setup()
print(f"p = {p}: |E(F_p)| = {[ch.Ebar(s).cardinality() for s in (1,-1)]}, ord(Gbar) = {ch.orders}, N = {ch.N}")
print(f"v_p(t_R) = {[t.valuation() for t in ch.tR]}, v = {ch.v}")

tors = ch.tors
# известные точки с рациональным u
known = {}
for T in tors:
    for m in range(-ch.N, ch.N + 1):
        P = m*G + T
        uu = S.u_of(P)
        if uu is None or uu in QQ:
            known.setdefault((T, m % ch.N), []).append(m)
print("known rational-u points (T, m0) -> m:", {(str(T[0]) if not T.is_zero() else 'O', m0): v for (T, m0), v in known.items()})

total = 0
rows = []
for T in tors:
    Tn = 'O' if T.is_zero() else str(T[0])
    for m0 in range(ch.N):
        Q = m0*G + T
        co, gamma, branches, phico = ch.theta(Q)
        # контроль: коэффициенты phi(x(Q+[t])) должны быть p-целыми в обоих вложениях
        intok = all(ch.emb(c, sg).valuation() >= 0 for c in phico for sg in (1, -1))
        n0, mu, tailok, vals = strassmann(co, ch.v, ch.p, M)
        # проверка теоретической оценки на вычисленных коэффициентах
        bndok = all(vals[j] >= j*ch.v - vpfact(j, ch.p) for j in range(M))
        # нули: известные целые s и "призрак" s = -1/2
        kn = known.get((T, m0), [])
        sint = sorted((m - m0)//ch.N for m in kn)
        gh = None
        if ch.N % 2 == 0 and m0 == ch.N//2:
            half = ch.Qp(-1)/2
            val = sum(co[j]*half^j for j in range(M))
            gh = val.valuation()
        rows.append((Tn, m0, n0, mu, tailok, intok, bndok, sint, gh))
        total += n0
        print(f"  T={Tn:>12} m0={m0}: Strassmann {n0}  (mu={mu}, tail_ok={tailok}, coeffs_int={intok}, bound_ok={bndok})"
              f"  known s={sint}" + (f"  theta(-1/2) val={gh}" if gh is not None else ""))
print(f"TOTAL zero bound = {total}")
kn_total = sum(len(v) for v in known.values())
print(f"known rational-u points = {kn_total}")
