good=[3, 5, 13, 37, 43, 197, 277, 283, 397, 557, 643, 683, 757, 797, 827, 907, 947, 997, 1237, 1453, 1693, 1867, 1933, 1987, 2267, 2357, 2477, 2557, 2917]
bad=[]
for l in good:
    E=EllipticCurve(QQ,[0,1+l^2,0,l^2,0])   # x(x+1)(x+l^2)
    r=pari(E.minimal_model().a_invariants()).ellinit().ellrank()
    if int(r[1])!=0: bad.append((l,r[:2]))
print("E_{1,ℓ} для 29 простых: ранг не 0 у", bad if bad else "никого — у всех 29 ранг 0 (PARI)")
