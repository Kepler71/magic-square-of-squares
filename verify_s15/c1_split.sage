# C1 бисэллиптична: M(z) = (z - c)/(c z - 1), c = 229/109. Если f1(M(z))(cz-1)^6 = kappa f1(z) с kappa = квадрат в Q,
# то подъём y -> sqrt(kappa) y/(cz-1)^3 определён над Q и Jac(C1) ~_Q E+ x E-  => rank = rank E+ + rank E- (mwrank).
import functools
print = functools.partial(print, flush=True)
R.<z> = QQ[]
Fr = R.fraction_field()
f1 = 65*(z^2 - 1)*(z^2 + 1)*(32161*z^2 - 49922*z + 32161)
for c in [QQ(229)/109, QQ(109)/229]:
    M = (z - c)/(c*z - 1)
    kap = Fr(f1(M) * (c*z - 1)^6) / Fr(f1)
    print(f"c = {c}: kappa = {kap}  square in Q: {QQ(kap).is_square() if kap in QQ else 'not constant'}")
c = QQ(229)/109
kap = QQ(Fr(f1((z - c)/(c*z - 1)) * (c*z - 1)^6) / Fr(f1))
lam = kap.sqrt() if kap.is_square() else None
assert lam is not None
# Инварианты: z -> M(z). Квотиент над Q: u = (z + M(z)) — не годится как параметр степени 2? z + M(z) = (c z^2 - 2z + c)/(c z - 1) ... берём
# стандартно: для инволюции tau: (z,y) -> (M(z), s*lam*y/(cz-1)^3), s = ±1, квотиенты C/<tau> и C/<tau*iota> — эллиптические.
# Строим через замену z = (w + p1)/(w + p2) с p1, p2 — неподвижными точками M, чтобы M стало w -> -w; неподвижные точки над Q(sqrt(-15)).
K.<a> = QuadraticField(-15)
RK.<w> = K[]
fix = (c*z^2 - 2*z + c).change_ring(K).roots(multiplicities=False)
p1, p2 = fix
print("fixed points of M:", fix)
# z = (p1 - p2*w... ) выберем z = (p1 + p2 w)/(1 + w): w = 0 -> p1, w = oo -> p2;  M меняет w -> -w (проверим)
zz = (p1 + p2*w)/(1 + w)
FK = RK.fraction_field()
Mw = FK((zz - c)/(c*zz - 1))
wimg = FK((Mw - p1)/(p2 - Mw))           # обратная замена: w' = (M(z) - p1)/(p2 - M(z))
print("M in w-coordinate is w -> -w:", wimg == -FK(w))
Fw = RK((FK(f1.change_ring(K)(zz)) * (1 + w)^6).numerator())
Fw = Fw / FK((1 + w)^6 * f1.change_ring(K)(zz)).denominator() if False else Fw
print("F(w) degree:", Fw.degree(), " even:", all(Fw[j] == 0 for j in range(1, Fw.degree() + 1, 2)))
# если F(w) = G(w^2) (чётный), то E+: y^2 = G(t), E-: y^2 = t G(t) над K; проверим, опускаются ли они на Q
Gt = RK([Fw[2*j] for j in range(Fw.degree()//2 + 1)])
def jac_quartic(g, F_):
    a_, b_, c_, d_, e_ = [g[j] for j in (4, 3, 2, 1, 0)]
    I = 12*a_*e_ - 3*b_*d_ + c_^2; J = 72*a_*c_*e_ + 9*b_*c_*d_ - 27*a_*d_^2 - 27*e_*b_^2 - 2*c_^3
    return EllipticCurve(F_, [-27*I, -27*J])
def jac_cubic(g, F_):
    a3, a2, a1, a0 = [g[j] for j in (3, 2, 1, 0)]
    return EllipticCurve(F_, [0, a2, 0, a1*a3, a0*a3^2])
Ep = jac_cubic(Gt, K) if Gt.degree() == 3 else jac_quartic(Gt, K)
Em = jac_quartic(RK(w*Gt), K) if Gt.degree() == 3 else None
print("E+ j:", Ep.j_invariant(), " E- j:", Em.j_invariant() if Em else None)
for nm, E_ in [("E+", Ep), ("E-", Em)]:
    if E_ is None: continue
    ds = E_.descend_to(QQ)
    print(f"{nm} descends to Q:", [e.ainvs() for e in ds])
    for e in ds:
        em = e.minimal_model()
        print(f"   {em.ainvs()}  cond {em.conductor().factor()}  rank(proof=True) = {em.rank(proof=True)}  bounds {em.rank_bounds()}")
