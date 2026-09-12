# Независимая проверка Richelot-изогении Codex (RESULT_FROM_CODEX_2026-09-11_RANK_ONE_PROVED.md §2).
# C0: Y^2 = 65 t(t^2-1)(109t-229)(229t-109)  (= C_J2 класса 5, закрутка 65);  C1: W^2 = 65(t^2-1)(t^2+1)(32161t^2-49922t+32161).
import functools
print = functools.partial(print, flush=True)
R.<t> = QQ[]
f0 = 65*t*(t^2 - 1)*(109*t - 229)*(229*t - 109)
f1 = 65*(t^2 - 1)*(t^2 + 1)*(32161*t^2 - 49922*t + 32161)
# (a) совпадение C0 с моей C_J2: D W^2 = f, D = 5*13*65*65 = 65^3 ~ 65 (mod квадраты)
f = t*(t^2 - 1)*(109*t - 229)*(229*t - 109)
print("(a) C0 = twist of my C_J2 by 65 (D = 65^3 ≡ 65 mod squares):", f0 == 65*f and (QQ(5*13*65*65)/65).is_square())
# (b) Richelot: F1 F2 F3 = f0 (F1 — «квадратичный» множитель с корнем в бесконечности)
F1 = 1622465*t; F2 = t^2 - 1; F3 = t^2 - QQ(64322)/24961*t + 1
print("(b) F1 F2 F3 == f0:", F1*F2*F3 == f0)
def coeffs(Fp): return [Fp[0], Fp[1], Fp[2]]
Dm = matrix(QQ, [coeffs(F1), coeffs(F2), coeffs(F3)])
delta = Dm.det(); print("    delta =", delta)
G1 = (F2.derivative()*F3 - F2*F3.derivative()); G2 = (F3.derivative()*F1 - F3*F1.derivative()); G3 = (F1.derivative()*F2 - F1*F2.derivative())
# стандарт (Cassels–Flynn / Bruin–Doerksen): дуальная кривая  delta * y^2 = G1 G2 G3
dual = delta * G1*G2*G3
q_, r_ = dual.quo_rem(f1)
print("    dual curve: y^2 = delta*G1*G2*G3 ; f1 divides it:", r_ == 0, " ratio =", q_, " ratio is a rational square:", q_.degree() == 0 and QQ(q_).is_square())
# (c) контроль изогении: у изогенных над Q якобианов совпадают многочлены Фробениуса над F_p при хорошей редукции
bad = set((2*f0.discriminant()*f1.discriminant()).numerator().prime_factors())
ok = True; n = 0
for p in primes(3, 400):
    if p in bad: continue
    C0p = HyperellipticCurve(f0.change_ring(GF(p))); C1p = HyperellipticCurve(f1.change_ring(GF(p)))
    if C0p.frobenius_polynomial() != C1p.frobenius_polynomial():
        ok = False; print("    mismatch at", p)
    n += 1
print(f"(c) Frobenius polynomials of Jac(C0), Jac(C1) equal at {n} good primes < 400:", ok)
# (d) корни f1: ±1, ±i, ρ, 1/ρ с ρ in Q(i) — алгебра Q x Q x Q(i) x Q(i) (малые группы классов)
print("(d) factorization of f1 over Q(i):", f1.change_ring(QuadraticField(-1, 'i')).factor())
# (e) рациональные точки C1 малой высоты (для контроля: Jac(C1)(Q) содержит точку бесконечного порядка)
pts = pari.hyperellratpoints(pari(f1), 10^4)
print("(e) rational points on C1 (height <= 1e4):", [(QQ(P[0]), QQ(P[1])) for P in pts][:12])
