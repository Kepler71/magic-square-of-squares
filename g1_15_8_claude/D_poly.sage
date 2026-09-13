S.<x>=QQ[]
F0=289*(1+x^2)-322*x; F8=289*(1+x^2)+322*x; L=49+529*x^2; U=529+49*x^2; F4=289*(1+x^2)
D=F0*F8*L
print("D =", D)
print("совпадает с многочленом из записки Codex:", D == 44182609*x^6+37608911*x^4+47287151*x^2+4092529)
print("Magma:", f"RationalPointsGenus2(HyperellipticCurve(Polynomial({list(D)})));")
H=HyperellipticCurve(D); print("род:", H.genus())
pts=[(x0, sqrt(D(x0))) for x0 in [QQ(a)/b for b in range(1,60) for a in range(-200,201)] if D(x0).is_square()]
print("точки с малой высотой (x=a/b, |a|≤200, b<60):", sorted(set(pts))[:10], "… плюс две на бесконечности, если старший коэф. — квадрат:", QQ(D.leading_coefficient()).is_square())
