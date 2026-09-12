# Точные минимальные сечения (высота 3/4) над K = Q(sqrt2, sqrt3) и точная решётка в char 0.
import functools, itertools
print = functools.partial(print, flush=True)
Kabs.<z> = NumberField(QQ["x"]([1, 0, -10, 0, 1]))          # z = sqrt2 + sqrt3
r2 = (z^3 - 9*z) / 2; r3 = (11*z - z^3) / 2
assert r2^2 == 2 and r3^2 == 3
K = Kabs
Rs.<s> = K[]
Fs = Rs.fraction_field()
P = 9*s^4 - 120*s^3 + 664*s^2 - 480*s - 1392
Q = 9*s^4 - 120*s^3 + 760*s^2 - 480*s - 2160
A4, A6 = -432*P, 3456*(s-6)*(3*s-2)*Q
E = EllipticCurve(Fs, [Fs(A4), Fs(A6)])

# минимальные сечения: орбиты размера 4 подсистемы a = 36 (см. sections_a36.sage), точки над K
Bq.<b,c,d,e,f,g> = PolynomialRing(QQ, order='degrevlex')
Bqt.<tq> = Bq[]
Xq = 36*tq^2 + b*tq + c; Yq = d*tq^3 + e*tq^2 + f*tq + g
Rq.<sq> = QQ[]
A4q = -432*(9*sq^4 - 120*sq^3 + 664*sq^2 - 480*sq - 1392)
A6q = 3456*(sq-6)*(3*sq-2)*(9*sq^4 - 120*sq^3 + 760*sq^2 - 480*sq - 2160)
J = Bq.ideal((Yq^2 - (Xq^3 + Bqt(A4q.change_ring(Bq))(tq)*Xq + Bqt(A6q.change_ring(Bq))(tq))).coefficients())
import time; t0 = time.time()
mp = J.minimal_associated_primes()
V = []
for Pj in mp:
    if Pj.vector_space_dimension() == 4:
        PK = Pj.change_ring(Bq.change_ring(K))
        pts = PK.variety()
        print("orbit of size 4: points over Q(sqrt2,sqrt3):", len(pts))
        V += [dict(a=K(36), **{str(k): v for k, v in pt.items()}) for pt in pts]
print("sections found:", len(V), f"({time.time()-t0:.1f}s)")
B = None
def nice(x):
    # запись элемента K в базисе 1, √2, √3, √6
    M = matrix(QQ, [list(K(1)), list(r2), list(r3), list(r2*r3)]).transpose()
    co = M.solve_right(vector(QQ, list(x)))
    parts = [(co[0], ""), (co[1], "√2"), (co[2], "√3"), (co[3], "√6")]
    return " + ".join(f"{q}{n}" if n else f"{q}" for q, n in parts if q != 0) or "0"

secs = []
for v in V:
    Xs = v['a']*s^2 + v['b']*s + v['c']; Ys = v['d']*s^3 + v['e']*s^2 + v['f']*s + v['g']
    assert Ys^2 == Xs^3 + A4*Xs + A6          # точная проверка в char 0
    secs.append(E(Fs(Xs), Fs(Ys)))

# высоты точно (те же формулы, что в heights_modp.sage; √2 ∈ K, так что узлы I2 рациональны над K)
r8 = [r for r, _ in (s^2 - 8).roots()]
assert len(r8) == 2
def at_inf(hh, wt): return Fs(hh(1/s) * s^wt)
A4i = at_inf(Fs(A4), 4)
def RO(R):
    num, den = R[0].numerator(), R[0].denominator()
    assert den.is_square()
    io = num.degree() - den.degree() - 2
    return den.sqrt().degree() + max(0, io // 2)
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
print("control: h(T) =", h(T), " h(S1) =", h(S1), " h(2S1) =", h(2*S1))
hs = [h(R) for R in secs]
print("heights over Q(√2,√3):", sorted(hs))
mins = [R for R, hv in zip(secs, hs) if hv == QQ(3)/4]
print("\nminimal sections (h = 3/4) over Q(√2,√3):", len(mins))
for R in mins:
    Xs, Ys = R[0].numerator(), R[1].numerator()
    print("  X =", " + ".join(f"({nice(Xs[i])})s^{i}" for i in range(2, -1, -1)))
# точный базис: S1 + два минимальных сечения с det = 1/4
free = [R for R, hv in zip(secs, hs) if hv != 0]
found = None
for P2, P3 in itertools.combinations(mins, 2):
    G = matrix(QQ, [[pair(x, y) for y in (S1, P2, P3)] for x in (S1, P2, P3)])
    if G.det() == QQ(1)/4:
        found = (P2, P3, G); break
if found:
    P2, P3, G = found
    print("\nEXACT basis {S1, P2, P3}, Gram matrix:\n", G, "\n det =", G.det())
    for nm, R in (("P2", P2), ("P3", P3)):
        Xs, Ys = R[0].numerator(), R[1].numerator()
        print(f"  {nm}: X = " + " + ".join(f"({nice(Xs[i])})s^{i}" for i in range(2, -1, -1)))
        print(f"      Y = " + " + ".join(f"({nice(Ys[i])})s^{i}" for i in range(3, -1, -1)))
else:
    print("no basis of form {S1, P2, P3} with P2, P3 minimal over Q(√2,√3)")
# ранг над Q(√2,√3): ранг Грама всех сечений, определённых над K (нижняя граница; целые сечения над K)
Gk = matrix(QQ, [[pair(x, y) for y in free] for x in free])
print("\nrank of lattice of integral sections over Q(√2,√3):", Gk.rank())
