# Исключения (как verify_s15/exceptions.sage) для всех сечений: t = 0, ±1, oo; w = 0; X = ±1; нули a±b, c±b.
secs = [(17,7,13), (7,1,5), (23,7,17), (31,17,25), (41,1,29), (47,23,37), (49,31,41), (73,17,53), (71,49,61), (89,23,65), (79,47,65)]
T = PolynomialRing(QQ, 'T').gen()
for (b, h, n) in secs:
    A, C = ZZ((h^2 + n^2)/2), ZZ((b^2 + n^2)/2); B = A + C
    a = A*(T^2 + 1)^2; bb = B*T*(T^2 - 1); c = C*(T^2 + 1)^2
    Ys = [a + bb, a - bb, c + bb, c - bb]
    ok = True; notes = []
    for t0 in (0, 1, -1):
        vals = [Y(t0) for Y in Ys]
        if all(QQ(v).is_square() for v in vals): ok = False; notes.append(f"t={t0} gives a point!")
    lead = [Y.leading_coefficient() for Y in Ys]
    if all(QQ(v).is_square() for v in lead): ok = False; notes.append("t=oo gives a point!")
    if (QQ(C)/A).is_square(): ok = False; notes.append("w=0 possible (C/A square)!")
    g4X1 = (C - A)*(A - C)
    if g4X1 >= 0: ok = False; notes.append("X=±1 not excluded!")
    zeros = [Y.roots(QQ, multiplicities=False) for Y in Ys]
    if any(zeros): notes.append(f"rational zeros of a±b/c±b: {zeros} (допустимы, Y_i = 0)")
    sqA, sqC = QQ(A).is_square(), QQ(C).is_square()
    print(f"({b},{h},{n}) A={A}{' (□)' if sqA else ''} C={C}{' (□)' if sqC else ''}: exceptions excluded: {ok} {notes}")
