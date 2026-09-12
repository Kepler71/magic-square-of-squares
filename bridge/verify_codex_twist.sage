# Проверка двух поправок Codex (RESULTS_G1_RANK_BRIDGE_2026-09-12.md):
#  (1) на D04: v² = F0·F4 всегда есть точка (t,v) = (1, m²+n²);
#  (2) правильный фактор — E_b = quadratic_twist(E0, b), b = s·m²·n², а не E0 сам по себе.
# Плюс сверка рангов E_b с таблицей Codex.
RT = PolynomialRing(QQ, 'T'); T = RT.gen()
pairs = [(11,4), (13,8), (15,1), (15,8), (16,5), (19,5), (19,16)]
codex_b_class = {(11,4):274, (13,8):466, (15,1):113, (15,8):2, (16,5):562, (19,5):193, (19,16):1234}
codex_rEb      = {(11,4):1,   (13,8):3,   (15,1):1,   (15,8):2, (16,5):2,   (19,5):1,   (19,16):1}
codex_rE0      = {(11,4):1,   (13,8):2,   (15,1):2,   (15,8):1, (16,5):1,   (19,5):2,   (19,16):2}

def rk(E):
    r = pari(E.minimal_model().ainvs()).ellinit().ellrank()
    return ZZ(r[0]), ZZ(r[1])

print(f"{'(m,n)':9} {'точка D04':>11} {'класс b':>9} {'b=Codex':>8} {'E_b=twist':>10} {'r(E0)':>7} {'r(E_b)':>8} {'сходится':>9}")
for (m, n) in pairs:
    m_, n_ = QQ(m), QQ(n)
    s = (m_^2 + n_^2)/2
    F0 = m_^2 + n_^2*T^2
    F4 = s*(1 + T^2)
    F8 = n_^2 + m_^2*T^2
    # (1) точка на D04 при t = 1
    val = F0(1)*F4(1)
    pt_ok = val.is_square() and sqrt(val) == m^2 + n^2

    # (2) твист
    b = s * m_^2 * n_^2
    a = m_^2/n_^2 + 1 + n_^2/m_^2
    E0 = EllipticCurve(QQ, [0, a, 0, a, 1])
    Eb = EllipticCurve(QQ, [0, b*a, 0, b^2*a, b^3])
    bclass = QQ(b).squarefree_part()
    tw_ok = Eb.is_isomorphic(E0.quadratic_twist(bclass))

    # сверка разложения H: Y² = F0·F4·F8 через карту (t,Y) -> (b t², b Y)
    H = RT(F0*F4*F8)
    X, Y = b*T^2, b
    lhs = (Y)^2 * H          # (bY)² = b²·Y² = b²·F0F4F8
    rhs = X^3 + b*a*X^2 + b^2*a*X + b^3
    map_ok = (lhs == rhs)

    r0, r1 = rk(E0); rb0, rb1 = rk(Eb)
    ok = (rb0 == rb1 == codex_rEb[(m,n)]) and (r0 == r1 == codex_rE0[(m,n)])
    print(f"({m},{n})".ljust(9)
          + ("да" if pt_ok else "НЕТ").rjust(11)
          + f"{bclass}".rjust(9)
          + ("да" if bclass == codex_b_class[(m,n)] else "НЕТ").rjust(8)
          + ("да" if tw_ok else "НЕТ").rjust(10)
          + f"{r0}..{r1}".rjust(7) + f"{rb0}..{rb1}".rjust(8)
          + ("да" if ok else "НЕТ").rjust(9))
    if not map_ok:
        print(f"    карта (t,Y)->(bt²,bY): тождество НЕ выполнено")
