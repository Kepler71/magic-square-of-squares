# Проверка универсальной теоремы Codex/Astra (RESULTS_FOR_CLAUDE_2026-09-12_G1.md, §2):
#   E_B: Y² = X(X − (m²−n²)²)(X − (m²+n²)²)  — модель фактора B = Jac(z² = F0·F8);
#   утверждение: кручение E_B(Q) = Z/2 × Z/4 для всех допустимых взаимно простых m,n.
# Плюс сверка: совпадает ли E_B с моей моделью через QuarticCurve(F0·F8, t=0) там, где та строилась.
import sys
load('/home/kep/magicKube/corners/corner_search.sage')
RT = PolynomialRing(QQ, 'T'); T = RT.gen()
def EB(m, n):
    m, n = ZZ(m), ZZ(n)
    return EllipticCurve(QQ, [0, -((m^2-n^2)^2 + (m^2+n^2)^2), 0, (m^2-n^2)^2*(m^2+n^2)^2, 0])
bad = []; tors = {}
MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20
pairs = [(m, n) for m in range(2, MAX+1) for n in range(1, m) if gcd(m, n) == 1 and m != n]
for (m, n) in pairs:
    E = EB(m, n)
    inv = tuple(E.torsion_subgroup().invariants())
    tors[inv] = tors.get(inv, 0) + 1
    if inv != (2, 4): bad.append((m, n, inv))
print(f"пар: {len(pairs)}; распределение кручения: {tors}")
print("нарушения Z/2×Z/4:", bad if bad else "нет")
# сверка моделей на пилотах, где строилась моя кривая (нужна квадратность (m²+n²)/2 — не обязательно, но пробуем t=0)
print("\nсверка модели Codex с моей (через квартику F0·F8, точка t = 0):")
for (m, n) in [(7,1), (17,7), (23,7), (41,1), (13,3), (19,11), (11,4), (15,1)]:
    F0 = m^2 + n^2*T^2; F8 = n^2 + m^2*T^2
    try:
        qc = QuarticCurve(RT(F0*F8), QQ(0)); mine = qc.E.minimal_model()
    except Exception as ex:
        print(f"  ({m},{n}): моя модель не строится — {type(ex).__name__}"); continue
    theirs = EB(m, n).minimal_model()
    r = pari(theirs.ainvs()).ellinit().ellrank()
    print(f"  ({m},{n}): изоморфны: {mine.is_isomorphic(theirs)};  ранг E_B = [{r[0]},{r[1]}];  кручение {theirs.torsion_subgroup().invariants()}")
