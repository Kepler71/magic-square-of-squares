# MW лучевой поверхности над K = Q(zeta24) = Q(i, √2, √3): точный базис и действие Галуа.
# Все 16 минимальных сечений (высота 3/4) лежат в подсистеме a = 36 (sections_a36.sage).
# Ранг MW над подполем k = ранг Gal(K/k)-инвариантов (все MW определено над K: целые сечения
# порождают MW — det 1/4, heights_modp.sage — и все они определены над K).
import functools, itertools, time
print = functools.partial(print, flush=True)
K.<z> = CyclotomicField(24)
I_, S2, S3 = K(-1).sqrt(), K(2).sqrt(), K(3).sqrt()
Rs.<s> = K[]
Fs = Rs.fraction_field()
P = 9*s^4 - 120*s^3 + 664*s^2 - 480*s - 1392
Q = 9*s^4 - 120*s^3 + 760*s^2 - 480*s - 2160
A4, A6 = -432*P, 3456*(s-6)*(3*s-2)*Q
E = EllipticCurve(Fs, [Fs(A4), Fs(A6)])

# --- минимальные сечения: минимальные простые подсистемы a=36 над Q, точки над K
Bq.<b,c,d,e,f,g> = PolynomialRing(QQ, order='degrevlex')
Bqt.<tq> = Bq[]
Rq.<sq> = QQ[]
A4q = -432*(9*sq^4 - 120*sq^3 + 664*sq^2 - 480*sq - 1392)
A6q = 3456*(sq-6)*(3*sq-2)*(9*sq^4 - 120*sq^3 + 760*sq^2 - 480*sq - 2160)
Xq = 36*tq^2 + b*tq + c; Yq = d*tq^3 + e*tq^2 + f*tq + g
J = Bq.ideal((Yq^2 - (Xq^3 + Bqt(A4q.change_ring(Bq))(tq)*Xq + Bqt(A6q.change_ring(Bq))(tq))).coefficients())
secs = []
for Pj in J.minimal_associated_primes():
    for pt in Pj.change_ring(Bq.change_ring(K)).variety():
        v = {str(k): val for k, val in pt.items()}
        Xs = 36*s^2 + v['b']*s + v['c']; Ys = v['d']*s^3 + v['e']*s^2 + v['f']*s + v['g']
        assert Ys^2 == Xs^3 + A4*Xs + A6
        secs.append(E(Fs(Xs), Fs(Ys)))
print("a=36 integral sections over K:", len(secs), "(ожидается 23 = 1+2+2+2+4+4+8)")

# --- высоты (как в heights_modp.sage)
r8 = [r for r, _ in (s^2 - 8).roots()]; assert len(r8) == 2
def at_inf(hh, wt): return Fs(hh(1/s) * s^wt)
A4i = at_inf(Fs(A4), 4)
def RO(R):
    num, den = R[0].numerator(), R[0].denominator(); assert den.is_square()
    return den.sqrt().degree() + max(0, (num.degree() - den.degree() - 2) // 2)
def node_fin(R, r):
    if R.is_zero() or R[0].denominator()(r) == 0: return False
    return R[1](r) == 0 and 3*R[0](r)^2 + A4(r) == 0
def node_inf(R):
    if R.is_zero(): return False
    Xi, Yi = at_inf(R[0], 2), at_inf(R[1], 3)
    if Xi.denominator()(0) == 0: return False
    return Yi(0) == 0 and 3*Xi(0)^2 + A4i(0) == 0
def h(R):
    if R.is_zero(): return QQ(0)
    cc = sum(QQ(1)/2 for r in r8 if node_fin(R, r))
    if node_inf(R): cc += QQ(3)/4 if node_inf(2*R) else QQ(1)
    return 2 + 2*RO(R) - cc
def pair(P1, P2): return (h(P1 + P2) - h(P1) - h(P2)) / 2

T  = E(Fs(36*s^2 - 240*s + 144), Fs(0))
S1 = E(Fs(72*s^2 - 96*s + 288), Fs(432*(s+2)*(s^2-12)))
hs = [h(R) for R in secs]
print("heights:", sorted(hs))
mins = [R for R, hv in zip(secs, hs) if hv == QQ(3)/4]
print("minimal sections:", len(mins), "(mod p: 16)")

# --- точный базис из минимальных сечений
basis = None
for trip in itertools.combinations(mins, 3):
    G = matrix(QQ, [[pair(x, y) for y in trip] for x in trip])
    if G.det() == QQ(1)/4:
        basis, GB = trip, G; break
assert basis is not None
print("\nEXACT basis (3 minimal sections), Gram:\n", GB, "\ndet =", GB.det())
def coords(R):
    v = GB.solve_right(vector(QQ, [pair(R, x) for x in basis]))
    assert all(x in ZZ for x in v), f"non-integral coords {v}"
    return v
print("S1 in basis coordinates:", coords(S1), "  (T: height 0, torsion)")
allc = [coords(R) for R in mins]
print("all 16 minimal sections have integral coordinates: True; distinct vectors:", len(set(tuple(x) for x in allc)))

# --- действие Галуа
G_K = K.automorphisms()
def act(sig, R):
    Xs, Ys = R[0], R[1]
    ap = lambda F_: Fs(F_.numerator().map_coefficients(sig)) / Fs(F_.denominator().map_coefficients(sig))
    return E(ap(Xs), ap(Ys))
mats = {}
for sig in G_K:
    M = matrix(ZZ, [coords(act(sig, x)) for x in basis]).transpose()   # столбцы = образы базиса
    mats[sig] = M
    assert M.transpose() * GB * M == GB, "Galois must preserve heights"
print("\nGalois acts by isometries of the lattice: True (all", len(G_K), "automorphisms)")
# подполя K и ранг инвариантов
def fixed_rank(H):
    Ms = [mats[sig] - 1 for sig in H]
    return 3 - block_matrix([[M] for M in Ms]).rank() if Ms else 3
sub = K.subfields()
print("\nMW rank over subfields of Q(i,√2,√3):")
for Fk, emb, _ in sorted(sub, key=lambda x: x[0].degree()):
    gen = emb(Fk.gen())
    H = [sig for sig in G_K if sig(gen) == gen]
    d = Fk.degree()
    name = "Q" if d == 1 else (f"Q(√{Fk.discriminant() if Fk.discriminant() % 4 == 1 else Fk.discriminant() // 4})" if d == 2 else f"deg {d}, disc {Fk.discriminant().factor()}")
    print(f"   {name:32s} rank {fixed_rank(H)}")
for nm, R in zip(["B1", "B2", "B3"], basis):
    fld = [Fk for Fk, emb, _ in sub if all(act(sig, R) == R for sig in G_K if sig(emb(Fk.gen())) == emb(Fk.gen()))]
    print(f"\n{nm}: X = {R[0]}\n    Y = {R[1]}")

# --- уточнение: подполя по их квадратичным подполям; изометрия с A3*
def qname(Fk):
    D = Fk.discriminant(); sq = D if D % 4 == 1 else D // 4
    return f"√{sq}"
print("\nMW rank over subfields (labelled by quadratic subfields):")
for Fk, emb, _ in sorted(sub, key=lambda x: x[0].degree()):
    gen = emb(Fk.gen())
    H = [sig for sig in G_K if sig(gen) == gen]
    quads = sorted(qname(F2) for F2, _, _ in Fk.subfields(2)) if Fk.degree() > 1 else []
    lab = "Q" if Fk.degree() == 1 else "Q(" + ",".join(quads if Fk.degree() == 4 else quads[:1] or [qname(Fk)]) + ")"
    if Fk.degree() == 8: lab = "Q(i,√2,√3)"
    print(f"   {lab:24s} rank {fixed_rank(H)}")
A3s = IntegralLattice("A3").dual_lattice() if False else None
GA3 = matrix(QQ, CartanMatrix(['A', 3])).inverse()          # Грам A3* в базисе фундаментальных весов
L1 = IntegralLattice(4*GB); L2 = IntegralLattice(4*GA3)
print("\nMWL isometric to A3* :", L1.is_isomorphic(L2) if hasattr(L1, 'is_isomorphic') else QuadraticForm(ZZ, 2*4*GB).is_globally_equivalent_to(QuadraticForm(ZZ, 2*4*GA3)))
# сечения базиса в базисе 1, i, √2, √3, i√2, i√3, √6, i√6
bas = [K(1), I_, S2, S3, I_*S2, I_*S3, S2*S3, I_*S2*S3]; nms = ["", "i", "√2", "√3", "i√2", "i√3", "√6", "i√6"]
Mb = matrix(QQ, [list(x) for x in bas]).transpose()
def nice(x):
    co = Mb.solve_right(vector(QQ, list(K(x))))
    return " + ".join(f"{q}{n}" for q, n in zip(co, nms) if q != 0).replace("+ -", "- ") or "0"
for nm, R in zip(["B1", "B2", "B3"], basis):
    Xs, Ys = R[0].numerator(), R[1].numerator()
    fixers = [Fk for Fk, emb, _ in sub if Fk.degree() < 8 and all(act(sig, R) == R for sig in G_K if sig(emb(Fk.gen())) == emb(Fk.gen()))]
    mind = min([Fk.degree() for Fk in fixers] + [8])
    print(f"\n{nm} (field of definition degree {mind}):")
    print("   X = " + " ; ".join(f"[s^{i}] {nice(Xs[i])}" for i in range(2, -1, -1)))
    print("   Y = " + " ; ".join(f"[s^{i}] {nice(Ys[i])}" for i in range(Ys.degree(), -1, -1)))
