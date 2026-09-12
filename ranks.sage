# Задача B: строгие ранги. rank(proof=True) — mwrank; при неудаче Sage бросает исключение,
# а не возвращает оценку. Контроль: rank_bounds, аналитический ранг, проверка генераторов и регулятора.
curves = [
    ("E1: W^2=T(T-218)(T-338)", [0, -556, 0, 73684, 0]),
    ("E: слой Бремнера lambda=2", [0, -180, 0, 5625, 0]),
    ("E5: твист E по 5", [0, -900, 0, 140625, 0]),
]
for name, a in curves:
    E = EllipticCurve(a)
    print("====", name, E.ainvs())
    print("  conductor:", E.conductor().factor(), "  j:", E.j_invariant())
    print("  torsion:", E.torsion_subgroup().invariants())
    try:
        r = E.rank(proof=True)
        print("  rank (proof=True, mwrank):", r)
    except Exception as ex:
        r = None
        print("  rank (proof=True) FAILED:", ex)
    print("  rank_bounds (2-descent):", E.rank_bounds())
    print("  analytic_rank (контроль, эвристика без BSD):", E.analytic_rank())
    G = E.gens(proof=True)
    print("  gens:", G)
    for P in G:
        assert P in E
    print("  regulator:", E.regulator())
    print("  saturated (index 1 in E(Q)/tors):", E.saturation(G)[1] == 1)
    print("  sha_an:", E.sha().an_numerical())

# Проверка утверждений записки
E = EllipticCurve([0, -180, 0, 5625, 0])
P1, P2 = E(36, 126), E(144, 252)
print("\n(36,126),(144,252) on E:", True, " height matrix det:", E.regulator_of_points([P1, P2]))
print("E(Q) == <(36,126),(144,252)> + tors ? index:", E.saturation([P1, P2])[1])
K.<s5> = QuadraticField(5)
print("rank E(Q(sqrt5)) = rank E + rank E^(5) =", E.rank() + E.quadratic_twist(5).rank())
print("E^(5) is isomorphic to given E5:", E.quadratic_twist(5).is_isomorphic(EllipticCurve([0, -900, 0, 140625, 0])))
