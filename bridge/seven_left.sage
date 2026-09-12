# Семь семейств G1, не исключённых у Codex/Astra: (11,4), (13,8), (15,1), (15,8), (16,5), (19,5), (19,16).
# У всех rank E_B = 1, поэтому перечисление через B невозможно. Но у кривой рода 5 есть и другие факторы:
#   A   = Jac(z² = F0·F4)      (рациональная точка при t = 0 есть, если (m²+n²)/2 — квадрат),
#   C48 = Jac(z² = F4·F8)      (то же условие),
#   E   = фактор рода 2 J_H ~ E²,  H: w² = F0·F4·F8;  E: v² = X³ + aX² + aX + 1, a = λ²+1+λ⁻², λ = m/n.
#       У E всегда есть точка (−1, 0) (2-кручение), так что она строится всегда.
# Ранг 0 у ЛЮБОГО фактора даёт конечность C(Q): C → H → E и C → C_ij — накрытия конечной степени.
import sys
RT = PolynomialRing(QQ, 'T'); T = RT.gen()
load('/home/kep/magicKube/corners/corner_search.sage')
def rk(E):
    r = pari(E.minimal_model().ainvs()).ellinit().ellrank(); return ZZ(r[0]), ZZ(r[1])
def EB(m, n):
    return EllipticCurve(QQ, [0, -((m^2-n^2)^2 + (m^2+n^2)^2), 0, (m^2-n^2)^2*(m^2+n^2)^2, 0])
pairs = [(11,4), (13,8), (15,1), (15,8), (16,5), (19,5), (19,16)]
print(f"{'(m,n)':10} {'rank B':>8} {'rank E (J_H ~ E²)':>20} {'rank A':>10} {'вывод':>28}")
for (m, n) in pairs:
    m_, n_ = QQ(m), QQ(n)
    lam = m_/n_; a = lam^2 + 1 + 1/lam^2
    E = EllipticCurve(QQ, [0, a, 0, a, 1]).minimal_model()
    rB = rk(EB(m, n)); rE = rk(E)
    F0 = m^2 + n^2*T^2; F4 = QQ(m^2+n^2)*(T^2+1)/2
    rA = None
    try:
        qc = QuarticCurve(RT(F0*F4), QQ(0)); rA = rk(qc.E.minimal_model())
    except Exception:
        rA = None
    zero = []
    if rB[1] == 0: zero.append('B')
    if rE[1] == 0: zero.append('E')
    if rA and rA[1] == 0: zero.append('A')
    verdict = ("ЗАКРЫТО через " + ",".join(zero)) if zero else "все факторы положительного ранга"
    print(f"({m},{n})".ljust(10) + f" {rB[0]}..{rB[1]}".rjust(8) + f" {rE[0]}..{rE[1]}".rjust(20) +
          (f" {rA[0]}..{rA[1]}".rjust(10) if rA else " нет точки".rjust(10)) + f"  {verdict}")
