#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
КАЛИБРОВКА реализации Теоремы 3.1 Фишера на ЕГО СОБСТВЕННОМ опубликованном
примере 3.4 (кривая 571a1).  Там L = Q[phi] -- кубическое ПОЛЕ, а не Q^3,
поэтому этот прогон проверяет формулу независимо от расщеплённого случая,
который используется для (11,4).

Опубликованные Фишером величины (arXiv:2208.14977, стр. 4-5):
    g1 = -11x^4 + 68x^3 z - 52x^2 z^2 - 164x z^3 - 64 z^4
    g2 =  -4x^4 - 60x^3 z - 232x^2 z^2 -  52x z^3 -  3 z^4
    g3 = -31x^4 - 78x^3 z +  32x^2 z^2 + 102x z^3 - 53 z^4
    I = 44608, J = 18842960, Delta = -2^12 * 571
    m = (1/9)(20 phi^2 - 8656 phi + 936032)
    gamma_1 = (4/9)(5x^2 - 16xz - 12z^2)
    вклад даёт только вещественное место, свидетель (15,4)
    спаривание нетривиально  =>  rank E(Q) = 0

Никакие из этих величин в расчёт не подставляются -- они ТОЛЬКО сверяются.
"""
from fractions import Fraction as Fr
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from refute_11_4_claude import (quartic_I, quartic_J, quartic_disc,
                                quartic_hessian, qeval, hilbert, is_square_Qp,
                                is_prime, factorint, as_int_sameclass)

def log(*a): print(*a); sys.stdout.flush()

# ---------------------------------------------------------------- L = Q[X]/f
class Lalg:
    """K[phi], phi^3 = 3I phi - J.  Elements are [c0,c1,c2]."""
    def __init__(self, I, J):
        self.I = Fr(I); self.J = Fr(J)
    def mul(self, a, b):
        r = [Fr(0)] * 5
        for i in range(3):
            for j in range(3):
                r[i + j] += a[i] * b[j]
        # reduce phi^3 = 3I phi - J ; phi^4 = 3I phi^2 - J phi
        r[1] += 3 * self.I * r[3]; r[0] -= self.J * r[3]; r[3] = Fr(0)
        r[2] += 3 * self.I * r[4]; r[1] -= self.J * r[4]; r[4] = Fr(0)
        return r[:3]
    def matrix(self, a):
        cols = []
        for k in range(3):
            ek = [Fr(0)] * 3; ek[k] = Fr(1)
            cols.append(self.mul(a, ek))
        return [[cols[k][r] for k in range(3)] for r in range(3)]
    def inv(self, a):
        M = self.matrix(a)
        A = [row[:] + [Fr(1) if i == j else Fr(0) for j in range(3)]
             for i, row in enumerate(M)]
        for c in range(3):
            piv = next(r for r in range(c, 3) if A[r][c] != 0)
            A[c], A[piv] = A[piv], A[c]
            pv = A[c][c]
            A[c] = [x / pv for x in A[c]]
            for r in range(3):
                if r != c and A[r][c] != 0:
                    f = A[r][c]
                    A[r] = [A[r][k] - f * A[c][k] for k in range(6)]
        # solve a * y = 1
        rhs = [Fr(1), Fr(0), Fr(0)]
        y = [sum(A[i][3 + k] * rhs[k] for k in range(3)) for i in range(3)]
        assert self.mul(a, y) == [Fr(1), Fr(0), Fr(0)]
        return y

def qmul(F, G):
    a0, a1, a2 = F; b0, b1, b2 = G
    return [a0*b0, a0*b1+a1*b0, a0*b2+a1*b1+a2*b0, a1*b2+a2*b1, a2*b2]

def run():
    g1 = [Fr(-11), Fr(68), Fr(-52), Fr(-164), Fr(-64)]
    g2 = [Fr(-4), Fr(-60), Fr(-232), Fr(-52), Fr(-3)]
    g3 = [Fr(-31), Fr(-78), Fr(32), Fr(102), Fr(-53)]
    log("=== КАЛИБРОВКА: пример 3.4 самого Фишера (кривая 571a1) ===")
    for nm, g in (("g1", g1), ("g2", g2), ("g3", g3)):
        log("  %s: I=%s J=%s Delta=%s" % (nm, quartic_I(g), quartic_J(g), quartic_disc(g)))
    I, J = quartic_I(g1), quartic_J(g1)
    assert I == 44608 and J == 18842960, "инварианты не сошлись со статьёй"
    assert quartic_disc(g1) == -(2**12) * 571
    assert quartic_I(g2) == I and quartic_J(g2) == J
    assert quartic_I(g3) == I and quartic_J(g3) == J
    log("  I, J, Delta совпали с опубликованными: 44608, 18842960, -2^12*571  OK")
    # L cubic field?
    f = lambda t: t**3 - 3*I*t + J
    roots = [t for t in range(-10000, 10001) if f(t) == 0]
    log("  рациональных корней X^3-3IX+J: %d  => L %s"
        % (len(roots), "кубическое ПОЛЕ (не Q^3)" if not roots else "расщепляется"))
    L = Lalg(I, J)
    phi = [Fr(0), Fr(1), Fr(0)]

    def z_inv(g):
        a, b, c, d, e = g
        return [(3*b*b - 8*a*c)/3, 4*a/3, Fr(0)]

    def H_form(g):
        h = quartic_hessian(g)
        G = [[h[t]/3, 4*g[t]/3, Fr(0)] for t in range(5)]
        Imphi2 = [Fr(I), Fr(0), Fr(-1)]                 # I - phi^2
        H = [G[0],
             [x/2 for x in G[1]],
             [G[2][k]/6 + Fr(2, 9)*Imphi2[k] for k in range(3)]]
        # verify G(1,0)*G = H^2 coefficientwise in L
        HH = [[Fr(0)]*3 for _ in range(5)]
        for i in range(3):
            for j in range(3):
                pr = L.mul(H[i], H[j])
                for k in range(3): HH[i+j][k] += pr[k]
        for t in range(5):
            lhs = L.mul(G[0], G[t])
            assert lhs == HH[t], "G(1,0)G = H^2 нарушено, t=%d" % t
        return H

    z1, z2, z3 = z_inv(g1), z_inv(g2), z_inv(g3)
    m_paper = [Fr(936032, 9), Fr(-8656, 9), Fr(20, 9)]
    prod = L.mul(L.mul(z1, z2), z3)
    log("  z(g1)z(g2)z(g3) = %s" % ([str(x) for x in prod],))
    log("  m (из статьи)^2 = %s" % ([str(x) for x in L.mul(m_paper, m_paper)],))
    assert L.mul(m_paper, m_paper) == prod, "m^2 != z1z2z3 -- расхождение со статьёй"
    log("  m^2 = z1 z2 z3 : OK (значение m из статьи подтверждено)")

    w = L.mul(L.mul(z2, z3), L.inv(m_paper))
    H1 = H_form(g1)
    gam = []
    for k in range(3):
        gam.append(L.mul(w, H1[k])[2])              # coefficient of phi^2
    log("  gamma_1 = %s x^2 + %s xz + %s z^2" % tuple(str(x) for x in gam))
    paper_gam = [Fr(4,9)*5, Fr(4,9)*(-16), Fr(4,9)*(-12)]
    log("  статья:    %s x^2 + %s xz + %s z^2" % tuple(str(x) for x in paper_gam))
    ok_gam = (gam == paper_gam)
    log("  СОВПАДЕНИЕ gamma_1 со статьёй: %s" % ok_gam)
    assert ok_gam

    a = g2[0]
    def gv(x, z): return gam[0]*x*x + gam[1]*x*z + gam[2]*z*z
    # places
    S = set([2, 3, 5, 7])
    for t in [a, quartic_disc(g1)] + [c for c in g1 if c] + [c for c in gam if c]:
        for q in factorint(abs(as_int_sameclass(t))): S.add(q)
    allP = sorted(S | {q for q in range(2, 2000) if is_prime(q)})
    res = {}
    for p in [0] + allP:
        pts = []
        cands = [(Fr(x), Fr(z)) for z in range(0, 25) for x in range(-25, 26)
                 if (x, z) != (0, 0)]
        if p:
            for k in range(1, 9):
                cands += [(Fr(1), Fr(p)**k), (Fr(p)**k, Fr(1)),
                          (Fr(1), Fr(p)**(-k)), (Fr(p)**(-k), Fr(1))]
        else:
            cands += [(Fr(1), Fr(1, 10**k)) for k in range(1, 30)]
        for (x, z) in cands:
            v = qeval(g1, x, z)
            if v == 0 or not is_square_Qp(v, p): continue
            if gv(x, z) == 0: continue
            pts.append((x, z))
            if len(pts) >= 4: break
        assert pts, "нет локальной точки при p=%s" % p
        vals = {hilbert(a, gv(x, z), p) for (x, z) in pts}
        assert len(vals) == 1, "локальный символ зависит от точки при p=%s" % p
        res[p] = vals.pop()
    pr = 1
    for p in res: pr *= res[p]
    minus = [("R" if p == 0 else p) for p in res if res[p] == -1]
    log("  мест проверено: %d ; места с -1: %s ; ПРОИЗВЕДЕНИЕ = %+d" % (len(res), minus, pr))
    log("  статья: вклад только от вещественного места, спаривание НЕТРИВИАЛЬНО")
    log("  свидетель статьи (15,4): g1 = %s (>0 ?), gamma_1 = %s (<0 ?)"
        % (qeval(g1, Fr(15), Fr(4)), gv(Fr(15), Fr(4))))
    ok = (pr == -1 and minus == ["R"])
    log("  ВЕРДИКТ КАЛИБРОВКИ: %s" % ("СОВПАЛО СО СТАТЬЁЙ" if ok else "!!! РАСХОЖДЕНИЕ"))
    return ok

if __name__ == "__main__":
    run()
