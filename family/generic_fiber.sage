# Общий слой семейства сечений: s — переменная. b = 1+2s-s^2, h = 1-2s-s^2, n = 1+s^2 (примитивная тройка при рациональном s).
R.<s> = QQ[]
b = 1 + 2*s - s^2; h = 1 - 2*s - s^2; n = 1 + s^2
assert b^2 + h^2 == 2*n^2
A = (h^2 + n^2)/2; C = (b^2 + n^2)/2
print("A =", A.factor()); print("C =", C.factor())
D2 = (C - A)*(C + A)
print("C^2 - A^2 =", D2.factor(), " степень", D2.degree())
# поле k_s = Q(sqrt(C^2-A^2)); квадратная часть
sq = prod([f^(e//2) for f, e in D2.factor()]); rest = D2 // sq^2
print("квадратная часть:", sq.factor() if sq != 1 else 1, "; свободная часть:", rest.factor(), " степень", rest.degree())
# кривая v^2 = rest(s): род и рациональные точки
if rest.degree() >= 1:
    try:
        Cv = HyperellipticCurve(rest) if rest.degree() >= 3 else None
        print("род кривой v^2 = rest:", Cv.genus() if Cv else 0)
    except Exception as ex:
        print("кривая:", ex)
    if rest.degree() == 2:
        a2, a1, a0 = rest[2], rest[1], rest[0]
        disc = a1^2 - 4*a2*a0
        print("коника v^2 = quadratic: дискриминант", disc, " старший коэффициент", a2, " квадрат?", a2.is_square() if a2 in QQ else None)
