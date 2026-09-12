load('/home/kep/magicKube/family/independent_chabauty/chab.sage')
import sys
b, h, n, D, u0, p = sys.argv[1:7]
prec = int(sys.argv[7]) if len(sys.argv) > 7 else 60
M = int(sys.argv[8]) if len(sys.argv) > 8 else 24
S = Section(ZZ(b), ZZ(h), ZZ(n), ZZ(D)); E = S.E; k = S.k
G = S.point_from_u(QQ(u0))
ch = Chab(S, G, ZZ(p), prec, M).setup()
Qp_ = ch.Qp; pp = ZZ(p)
print(f"N = {ch.N}, orders = {ch.orders}, v = {ch.v}")

# ---------- 1. редукции классов ----------
for T in ch.tors:
    Tn = 'O' if T.is_zero() else str(T[0])
    for m0 in range(ch.N):
        Q = m0*G + T
        r1, r2 = ch.redpt(Q, 1), ch.redpt(Q, -1)
        print(f"  class T={Tn} m0={m0}: red1={r1}, red2={r2}")

# ---------- 2. сквозная проверка ряда theta на целых s ----------
print("end-to-end test: theta(s) vs прямое вычисление sigma1 x - sigma2 x")
bad = 0
for T in ch.tors[:4]:
    for m0 in range(ch.N):
        Q = m0*G + T
        co, gamma, branches, phico = ch.theta(Q)
        for s in [0, 1, -1, 2, -2, 3]:
            P = (m0 + ch.N*s)*G + T
            ser = sum(co[j]*Qp_(s)^j for j in range(M))
            if P.is_zero():
                direct = Qp_(0)
            else:
                x1, x2 = ch.emb(P[0], 1), ch.emb(P[0], -1)
                direct = 1/(x1 - gamma) - 1/(x2 - gamma)
            d = ser - direct
            okk = (d == 0) or d.valuation() >= min(ser.valuation() if ser != 0 else 99,
                                                   direct.valuation() if direct != 0 else 99) + 10
            if not okk:
                bad += 1
                print(f"  MISMATCH T={'O' if T.is_zero() else T[0]} m0={m0} s={s}: ser={ser}, direct={direct}")
print("end-to-end mismatches:", bad)
