good=[3, 5, 13, 37, 43, 197, 277, 283, 397, 557, 643, 683, 757, 797, 827, 907, 947, 997, 1237, 1453, 1693, 1867, 1933, 1987, 2267, 2357, 2477, 2557, 2917]
R.<p>=QQ[]
bad=[]
for l in good:
    Eml=EllipticCurve(QQ,[0,1+l^2,0,l^2,0])
    for k,cs in ((QQ(1)/l,(1,QQ(1)/l,1+QQ(1)/l)), (QQ(l-1)/l,(1,-(1-QQ(l-1)/l),QQ(l-1)/l))):
        f=prod(1+c*p for c in cs); A,B,C=f[3],f[2],f[1]
        Es=EllipticCurve(QQ,[0,B,0,A*C,A^2])
        if not Es.is_isomorphic(Eml): bad.append((l,k))
print("кривые наклонов 1/ℓ и (ℓ−1)/ℓ, не изоморфные E_{1,ℓ} над Q:", bad if bad else "нет — все 58 изоморфны")
