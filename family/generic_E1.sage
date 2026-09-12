# Общий слой E1 над F = Q(s)(sqrt(-8s(s^2-1))) — символически, по той же схеме, что family/pipeline.sage.
K.<s> = FunctionField(QQ)
Ry.<Y> = K[]
L.<r> = K.extension(Y^2 - (-8*s^3 + 8*s))       # r^2 = -8s(s^2-1);  k_s = Q(sqrt(rest)) при специализации
A = s^4 + 2*s^3 + 2*s^2 - 2*s + 1
C = s^4 - 2*s^3 + 2*s^2 + 2*s + 1
bet = (C - A)/(C + A)
print("beta' = (C-A)/(C+A) =", bet)
# sqrt(beta'): beta' = (C-A)/(C+A); (C-A)(C+A) = (s^2+1)^2 * rest  =>  sqrt(beta') = r / ((s^2+1) * ... )
num = C - A
sb2 = bet
# sqrt(bet) = (C-A)/sqrt((C-A)(C+A)) = (C-A)/((s^2+1) r) * ... проверим
cand = L(num) / (L(s^2 + 1) * r)
print("кандидат sqrt(beta'):", "квадрат совпал:", cand^2 == L(bet))
sb = cand
Ru.<uu> = PolynomialRing(L)
Ap = -(1 + bet^2)
cz_sym = None
# cz и кубика — как в pipeline: cub = D * cz * (uu^2 + Ap - 2 bet) * (uu + 2 sb)
# D и cz — рациональные множители, на структуру кривой влияют только как квадратичная закрутка; берём D*cz = 1
cub = (uu^2 + L(Ap) - 2*L(bet)) * (uu + 2*sb)
a3, a2, a1, a0 = [cub[j] for j in (3, 2, 1, 0)]
E1 = EllipticCurve(L, [0, a2, 0, a1*a3, a0*a3^2])
print("E1 над F: ainvs =", E1.ainvs())
j = E1.j_invariant()
print("j(E1) в K (т.е. рационально по s, без r):", j in K)
if j in K: print("j =", K(j).factor() if K(j) != 0 else 0)
print("2-кручение: корни 2-деления")
for e_, m_ in E1.two_division_polynomial().roots():
    print("   ", e_)
