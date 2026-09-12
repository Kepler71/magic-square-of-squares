# Задача A: геометрический род H0(k,u).
# Три независимых метода + контроль каждого на кривых с известным родом.
#   M1: Singular genus (normal.lib) через Curve(...).geometric_genus(), характеристика 0
#   M2: алгебраическое функциональное поле над GF(p) (алгоритм Хесса в Sage, независим от Singular)
#   M3: Риман-Гурвиц по проекции на k над GF(p): 2g-2 = -2n + sum(e-1) по разветвлённым местам
# Род редукции по почти всем p равен роду в характеристике 0; поэтому M2/M3 берём на нескольких p.
import sys, time

R.<k,u> = QQ[]
H0 = R(sage_eval(open('h0.txt').read(), locals={'k': k, 'u': u}))
assert H0.degree(k) == 8 and H0.degree(u) == 8

def M1(f):
    return Curve(f).geometric_genus()

def ff(f, p):
    K.<T> = FunctionField(GF(p))
    S.<Y> = K[]
    g = S(f.change_ring(GF(p))(T, Y))
    g = g / g.leading_coefficient()
    if not g.is_irreducible():
        return None, None
    L.<y> = K.extension(g)
    return K, L

def M2(f, p):
    K, L = ff(f, p)
    return None if L is None else L.genus()

def M3(f, p):
    # сумма (e-1) по всем местам L над всеми местами K, где возможна разветвлённость
    K, L = ff(f, p)
    if L is None:
        return None
    n = L.degree()
    assert p > n, "нужна ручная (tame) ситуация p > deg"
    O, Oinf = L.maximal_order(), L.maximal_order_infinite()
    Kmo = K.maximal_order()
    # кандидаты в точки ветвления конечной части: неприводимые делители дискриминанта по u
    # и старшего коэффициента по u (надмножество, лишние дадут e=1 и вклад 0)
    A.<kk, uu> = GF(p)[]
    fA = A(f.change_ring(GF(p))(kk, uu))
    Ru = K._ring                      # GF(p)[T]
    to_T = lambda h: Ru(h.univariate_polynomial()) if hasattr(h, 'univariate_polynomial') else Ru(h)
    Dk = to_T(fA.discriminant(uu))
    lc = to_T(fA.polynomial(uu).leading_coefficient())
    assert Dk != 0
    total = 0
    # decomposition(prime) -> [(prime ideal, residue degree f, ramification index e)]
    # степень места над q равна f*deg(q); вклад в дифференту (tame) (e-1)*f*deg(q)
    for q in {q for q, _ in (Dk * lc).factor()}:
        total += sum((e - 1) * fr * q.degree() for _, fr, e in O.decomposition(Kmo.ideal(K(q))))
    for _, fr, e in Oinf.decomposition():   # место k = oo имеет степень 1
        total += (e - 1) * fr
    g2 = -2 * n + total
    assert g2 % 2 == 0
    return g2 // 2 + 1

# ---------------- контрольные кривые с известным родом ----------------
import random
random.seed(int(1))
Rr = R
rnd88 = sum(ZZ(random.randint(-9, 9)) * k^i * u^j for i in range(9) for j in range(9))
controls = [
    ("Ферма x^5+y^5=1", k^5 + u^5 - 1, 6),
    ("Ферма x^7+y^7=1", k^7 + u^7 - 1, 15),
    ("гиперэлл. u^2=f8(k)", u^2 - (k^8 + 3*k^5 - 7*k^3 + k + 11), 3),
    ("genus-2 quotient Y^2=(109X^2-229)(X^4-1)", u^2 - (109*k^2 - 229)*(k^4 - 1), 2),
    ("случайная бистепени (8,8)", rnd88, 49),
    ("гиперэлл. u^2 = k^20-k+1", u^2 - (k^20 - k + 1), 9),
]
primes = [101, 103, 107]

args = sys.argv[1:]
do_m1 = 'nom1' not in args

print("=== контроль ===")
for name, f, gtrue in controls:
    t0 = time.time()
    r1 = M1(f) if do_m1 else '-'
    r2 = [M2(f, p) for p in primes]
    r3 = [M3(f, p) for p in primes]
    ok = all(v == gtrue for v in ([r1] if do_m1 else []) + r2 + r3)
    print(f"{name:45s} true={gtrue:3d}  M1={r1}  M2={r2}  M3={r3}  {'OK' if ok else 'MISMATCH'}  ({time.time()-t0:.1f}s)")
    sys.stdout.flush()

print("\n=== H0 ===")
for p in [101, 103, 107, 109, 113, 127, 131]:
    t0 = time.time()
    print(f"p={p}: M2={M2(H0, p)}  M3={M3(H0, p)}  ({time.time()-t0:.1f}s)")
    sys.stdout.flush()
if do_m1:
    t0 = time.time()
    print("M1 (char 0):", M1(H0), f"({time.time()-t0:.1f}s)")
