#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Дополнительные попытки сломать сертификат (11,4):
  G1. ELS-свидетели для ВСЕХ трёх квартик g1,g2,g3 поимённо во всех плохих местах.
  G2. Независимость от выбора m: Замечание 3.2(vi) разрешает ЛЮБОЙ корень
      m из z1z2z3 = m^2.  Их 8 (по знаку в каждой компоненте).  Все 8 обязаны
      дать один ответ.
  G3. Независимость от представителя: случайные собственные эквивалентности
      g -> g o gamma, gamma in SL_2(Z), для каждой из трёх квартик.
  G4. Использование g3(1,0) вместо g2(1,0) (Замечание 3.2(v)) -- та самая
      ловушка, на которой у (19,5) получалось +1 при неполном наборе мест.
  G5. Полная матрица спаривания на базисе Sel^2, её ранг и вытекающая
      граница на ранг E(Q) (сравнение с независимо известным ранг = 1).
"""
import sys, random
from fractions import Fraction as Fr
from math import isqrt

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import refute_11_4_claude as R

log = R.log

def main():
    m, n = 11, 4
    s = Fr(m*m + n*n, 2); b = s*m*m*n*n
    roots = [int(4*(-b)), int(4*(-s*m**4)), int(4*(-s*n**4))]
    C = R.Curve(roots)
    assert C.check()
    S = [2, 3, 5, 7, 11, 137]
    sel, _ = R.selmer_group(C, S, quiet=True)
    selset = set(tuple(x) for x in sel)
    d1 = (1, 274, 274)
    d2 = (2055, 137, 15)
    d3 = tuple(R.sqfree(d1[i]*d2[i]) for i in range(3))
    log("delta = %s ; g2-класс = %s ; g3-класс = %s" % (d1, d2, d3))
    assert d1 in selset and d2 in selset and d3 in selset
    g1 = R.build_quartic(C, d1, selset, verbose=False)
    g2 = R.build_quartic(C, d2, selset, verbose=False)
    g3 = R.build_quartic(C, d3, selset, verbose=False)
    for nm, g in (("g1", g1), ("g2", g2), ("g3", g3)):
        log("%s = %s" % (nm, [str(x) for x in g]))
        log("    I,J верны: %s ; класс z(g) = %s"
            % (R.quartic_I(g) == C.I and R.quartic_J(g) == C.J,
               tuple(R.sqfree(t) for t in R.z_invariant(C, g))))

    # ---------------- G1 : ELS witnesses ----------------
    log("\nG1. ЯВНЫЕ локальные точки (g(x,z) -- ненулевой квадрат в Q_v) для всех трёх:")
    places = [0] + S
    for nm, g in (("g1", g1), ("g2", g2), ("g3", g3)):
        row = []
        for p in places:
            got = None
            for z in range(0, 60):
                for x in range(-60, 61):
                    if x == 0 and z == 0: continue
                    v = R.qeval(g, Fr(x), Fr(z))
                    if v != 0 and R.is_square_Qp(v, p): got = (x, z); break
                if got: break
            if got is None and p:
                for k in range(1, 10):
                    for cand in ((1, p**k), (p**k, 1), (1, -(p**k))):
                        v = R.qeval(g, Fr(cand[0]), Fr(cand[1]))
                        if v != 0 and R.is_square_Qp(v, p): got = cand; break
                    if got: break
            row.append("%s:%s" % ("R" if p == 0 else p, got))
            assert got is not None, "нет свидетеля для %s при p=%s" % (nm, p)
        log("   %s : %s" % (nm, "  ".join(row)))
    log("   остальные p: p не делит 2*Delta(g) => y^2=g -- гладкая кривая рода 1")
    log("   над F_p, #C(F_p) >= p+1-2sqrt(p) > 0 при p>=5, точка поднимается по Гензелю.")
    log("   Delta(g1) = %s" % R.quartic_disc(g1))
    dsc = R.quartic_disc(g1)
    log("   простые делители Delta(g1): %s"
        % sorted(R.factorint(abs(R.as_int_sameclass(dsc))).keys()))

    base, _, _ = R.ct_pairing(C, g1, g2, g3, verbose=False)
    log("\nбазовое значение <delta, g2>_CT = %+d" % base)

    # ---------------- G2 : all 8 square roots m ----------------
    log("\nG2. Все 8 корней m из z1z2z3 (Замечание 3.2(vi)) обязаны дать один ответ:")
    z1 = R.z_invariant(C, g1); z2 = R.z_invariant(C, g2); z3 = R.z_invariant(C, g3)
    H1 = R.H_form(C, g1)
    a = g2[0]
    vals = []
    for sg in [(1,1,1),(1,1,-1),(1,-1,1),(1,-1,-1),(-1,1,1),(-1,1,-1),(-1,-1,1),(-1,-1,-1)]:
        mm = []
        for i in range(3):
            pr = z1[i]*z2[i]*z3[i]
            r1 = isqrt(pr.numerator); r2 = isqrt(pr.denominator)
            assert r1*r1 == pr.numerator and r2*r2 == pr.denominator
            mm.append(sg[i]*Fr(r1, r2))
        w = [z2[i]*z3[i]/mm[i] for i in range(3)]
        gam = [R.L_interp_gamma(C, [w[i]*H1[k][i] for i in range(3)]) for k in range(3)]
        # product over the full place set
        Sset = set([2,3,5,7])
        for t in [a, R.quartic_disc(g1)] + [c for c in g1 if c] + [c for c in gam if c]:
            for q in R.factorint(abs(R.as_int_sameclass(t))): Sset.add(q)
        allP = sorted(Sset | {q for q in range(2, R.PLIMIT) if R.is_prime(q)})
        pr_ = 1; minus = []
        for p in [0] + allP:
            pts = R.local_points(g1, gam, p, want=2)
            assert pts, "нет локальной точки p=%s" % p
            h = R.hilbert(a, R.gamma_eval(gam, *pts[0]), p)
            if h == -1: minus.append("R" if p == 0 else p)
            pr_ *= h
        vals.append(pr_)
        log("   знаки m = %-14s -> %+d   минусы: %s" % (str(sg), pr_, minus))
    log("   все 8 совпали: %s  %s" % (len(set(vals)) == 1,
        "OK" if len(set(vals)) == 1 and vals[0] == base else "!!! РАСХОЖДЕНИЕ"))

    # ---------------- G3 : random proper equivalences ----------------
    log("\nG3. Случайные собственные эквивалентности (g -> g o gamma, gamma in SL_2(Z)):")
    random.seed(7)
    ok3 = True
    for trial in range(6):
        def rnd():
            while True:
                p_, q_, r_, s_ = [random.randint(-3, 3) for _ in range(4)]
                if p_*s_ - q_*r_ == 1: return (p_, q_, r_, s_)
        G1 = R.subst(g1, *rnd()); G2 = R.subst(g2, *rnd()); G3 = R.subst(g3, *rnd())
        if any(R.quartic_I(g) != C.I or R.quartic_J(g) != C.J for g in (G1, G2, G3)):
            log("   попытка %d: инварианты не сохранились -- пропуск" % trial); continue
        if G2[0] == 0:
            log("   попытка %d: g2(1,0)=0 -- пропуск (Зам. 3.2(iii))" % trial); continue
        try:
            G1 = R.make_z_unit(C, G1); G2 = R.make_z_unit(C, G2); G3 = R.make_z_unit(C, G3)
            if None in (G1, G2, G3):
                log("   попытка %d: z(g) делитель нуля -- пропуск" % trial); continue
            v, res, info = R.ct_pairing(C, G1, G2, G3, verbose=False)
        except Exception as ex:
            log("   попытка %d: %s" % (trial, ex)); continue
        minus = sorted([("R" if p == 0 else p) for p in res if res[p] == -1],
                       key=lambda t: (t == "R", 0 if t == "R" else t))
        good = (v == base)
        ok3 &= good
        log("   попытка %d: a=%s -> %+d  минусы: %s  %s"
            % (trial, G2[0], v, minus, "OK" if good else "!!! РАСХОЖДЕНИЕ"))
    log("   вывод: %s" % ("значение инвариантно относительно представителя"
                          if ok3 else "!!! ЗНАЧЕНИЕ ЗАВИСИТ ОТ ПРЕДСТАВИТЕЛЯ"))

    # ---------------- G4 : g3(1,0) instead of g2(1,0) ----------------
    log("\nG4. Замечание 3.2(v): можно взять g3(1,0) вместо g2(1,0) -- это та ловушка,")
    log("    на которой у (19,5) при НЕПОЛНОМ наборе мест получалось +1.")
    gam, _, _ = R.gamma1_form(C, g1, g2, g3)
    for nm, aa in (("g2(1,0)", g2[0]), ("g3(1,0)", g3[0])):
        Sset = set([2, 3, 5, 7])
        for t in [aa, R.quartic_disc(g1)] + [c for c in g1 if c] + [c for c in gam if c]:
            for q in R.factorint(abs(R.as_int_sameclass(t))): Sset.add(q)
        naive = set([2, 3, 5, 7])
        for t in [g2[0], R.quartic_disc(g1)] + [c for c in g1 if c] + [c for c in gam if c]:
            for q in R.factorint(abs(R.as_int_sameclass(t))): naive.add(q)
        full = sorted(Sset | {q for q in range(2, R.PLIMIT) if R.is_prime(q)})
        pf = 1; pn = 1
        for p in [0] + full:
            pts = R.local_points(g1, gam, p, want=2)
            h = R.hilbert(aa, R.gamma_eval(gam, *pts[0]), p)
            pf *= h
            if p == 0 or p in naive: pn *= h
        log("    a = %-9s : полный набор -> %+d ; 'наивный' набор (от g2) -> %+d  %s"
            % (nm, pf, pn, "" if pf == pn else "<-- ЛОВУШКА ВОСПРОИЗВЕДЕНА"))
        log("      простые, добавляемые самим a: %s"
            % sorted(set(R.factorint(abs(R.as_int_sameclass(aa))).keys()) - naive))

    # ---------------- G5 : full CTP matrix and rank bound ----------------
    log("\nG5. Полная матрица спаривания на базисе Sel^2 и вытекающая граница ранга:")
    gens = []
    span = {(1, 1, 1)}
    for c in sorted(selset):
        if c in span: continue
        gens.append(c)
        span = {tuple(R.sqfree(x[i]*y[i]) for i in range(3))
                for x in span | {c} for y in span | {c}} | span | {c}
        span = set()
        # regenerate span from gens
        cur = {(1, 1, 1)}
        for g_ in gens:
            cur = cur | {tuple(R.sqfree(t[i]*g_[i]) for i in range(3)) for t in cur}
        span = cur
        if len(span) == len(selset): break
    log("   базис Sel^2 (%d элементов): %s" % (len(gens), gens))
    quart = {}
    def getq(cl):
        cl = tuple(cl)
        if cl not in quart:
            quart[cl] = R.build_quartic(C, cl, selset, verbose=False)
        return quart[cl]
    N = len(gens)
    M = [[0]*N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            A, B = gens[i], gens[j]
            Cc = tuple(R.sqfree(A[k]*B[k]) for k in range(3))
            gA, gB, gC = getq(A), getq(B), getq(Cc)
            if None in (gA, gB, gC) or gB[0] == 0:
                M[i][j] = None; continue
            v, _, _ = R.ct_pairing(C, gA, gB, gC, verbose=False)
            M[i][j] = 0 if v == 1 else 1
    for row in M: log("   %s" % row)
    sym = all(M[i][j] == M[j][i] for i in range(N) for j in range(N)
              if M[i][j] is not None and M[j][i] is not None)
    diag = all(M[i][i] in (0, None) for i in range(N))
    log("   симметрична: %s ; нулевая диагональ: %s" % (sym, diag))
    # F2-rank
    A = [row[:] for row in M]
    rank = 0; rows = N
    cols = list(range(N)); r = 0
    for c in range(N):
        piv = None
        for rr in range(r, N):
            if A[rr][c] == 1: piv = rr; break
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        for rr in range(N):
            if rr != r and A[rr][c] == 1:
                A[rr] = [(A[rr][k] ^ A[r][k]) for k in range(N)]
        r += 1
    rank = r
    log("   ранг матрицы CTP над F_2 = %d (обязан быть чётным)" % rank)
    log("   dim Sel^2 = %d, dim ker CTP = %d" % (N, N - rank))
    log("   im(E(Q)/2E(Q)) lies in ker CTP  =>  rank E(Q) <= %d - dim E(Q)[2] = %d"
        % (N - rank, N - rank - 2))
    log("   (независимо известно: rank E(Q) = 1 -- PARI ellrank и аналитический ранг;")
    log("    согласие является контролем, в сам сертификат ранг не входит)")

if __name__ == "__main__":
    main()
