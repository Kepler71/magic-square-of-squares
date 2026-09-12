# Часть 2: РАНГ. Единственное место, от которого зависит исключение (11,4).
def hdr(t):
    print("\n" + "="*78); print(t); print("="*78)

M, N = ZZ(11), ZZ(4)
S = QQ(M^2+N^2)/2; B = S*M^2*N^2
R1,R2,R3 = QQ(-B), QQ(-S*M^4), QQ(-S*N^4)
Rx.<XX> = PolynomialRing(QQ)
cub = (XX-R1)*(XX-R2)*(XX-R3); c = cub.coefficients(sparse=False)
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
Emin = E.minimal_model()
print("E    =", E)
print("Emin =", Emin)
print("a-инварианты Emin =", Emin.a_invariants())
print("проводник =", Emin.conductor(), "=", factor(Emin.conductor()))
print("изоморфизм E ~ Emin:", E.is_isomorphic(Emin))

hdr("D1. PARI ellrank НАПРЯМУЮ, разные effort")
from sage.libs.pari import pari
pari.allocatemem(4*10^9, silent=True)
ep = pari.ellinit(list(Emin.a_invariants()))
for eff in [0,1,2,3]:
    try:
        alarm(600)
        res = pari.ellrank(ep, eff)
        cancel_alarm()
        print("  ellrank(effort=%d) -> %s" % (eff, res))
    except Exception as ex:
        try: cancel_alarm()
        except: pass
        print("  ellrank(effort=%d) FAILED/TIMEOUT: %s" % (eff, ex))
print("  формат: [нижняя граница ранга, верхняя граница ранга, s, найденные точки]")
print("  PARI версия:", pari.version())

hdr("D2. root number / parity")
print("  Emin.root_number() =", Emin.root_number())
print("  => по теореме о 2-чётности (Докчицер): rank + corank Sha[2^oo] нечётен"
      if Emin.root_number()==-1 else "  => чётен")

hdr("D3. mwrank: dim Sel_2 и верхняя граница")
try:
    print("  E.selmer_rank()  =", E.selmer_rank())
except Exception as ex: print("  selmer_rank FAILED:", ex)
try:
    print("  E.rank_bound()   =", E.rank_bound())
except Exception as ex: print("  rank_bound FAILED:", ex)
try:
    print("  E.rank_bounds()  =", E.rank_bounds())
except Exception as ex: print("  rank_bounds FAILED:", ex)
try:
    print("  Emin.selmer_rank() =", Emin.selmer_rank(), " Emin.rank_bound() =", Emin.rank_bound())
except Exception as ex: print("  Emin FAILED:", ex)

hdr("D4. аналитический ранг")
for alg in ['pari','sympow','rubinstein']:
    try:
        alarm(600)
        v = Emin.analytic_rank(algorithm=alg)
        cancel_alarm()
        print("  analytic_rank(%s) = %s" % (alg, v))
    except Exception as ex:
        try: cancel_alarm()
        except: pass
        print("  analytic_rank(%s) FAILED: %s" % (alg, ex))
try:
    alarm(900)
    print("  analytic_rank_upper_bound() =", Emin.analytic_rank_upper_bound())
    cancel_alarm()
except Exception as ex:
    try: cancel_alarm()
    except: pass
    print("  analytic_rank_upper_bound FAILED:", ex)

hdr("D5. Sha аналитический (BSD-условно) — индикатор размера Sha")
try:
    alarm(900)
    print("  sha.an_numerical() =", Emin.sha().an_numerical())
    cancel_alarm()
except Exception as ex:
    try: cancel_alarm()
    except: pass
    print("  an_numerical FAILED:", ex)

hdr("D6. поиск дополнительных независимых точек (если ранг 3 — они должны найтись)")
print("  gens(proof=False) =", E.gens(proof=False))
for h in [12, 16, 20]:
    try:
        alarm(600)
        pts = Emin.point_search(h, rank_bound=4)
        cancel_alarm()
        if pts:
            sat, idx, reg = Emin.saturation(pts)
            print("  point_search(%d): %d точек -> независимых после saturation: %d" % (h, len(pts), len(sat)))
        else:
            print("  point_search(%d): точек не найдено" % h)
    except Exception as ex:
        try: cancel_alarm()
        except: pass
        print("  point_search(%d) FAILED: %s" % (h, ex))
print("ГОТОВО")
