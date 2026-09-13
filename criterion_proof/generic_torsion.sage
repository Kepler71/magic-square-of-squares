# Общий класс клеток {a, b, a+b}: точки кручения кривой Y^2=(1+ap)(1+bp)(1+(a+b)p) над Q(a,b)
F.<a,b> = QQ[]
K = F.fraction_field(); a,b = K(a),K(b)
R.<p> = K[]
f = (1+a*p)*(1+b*p)*(1+(a+b)*p); A,Bc,Cc = f[3],f[2],f[1]
E = EllipticCurve(K,[0,Bc,0,A*Cc,A^2])
P0 = E(0,A); T = [E(-A/c,0) for c in (a,b,a+b)]
pts = [P0] + T + [P0+t for t in T]
for P in pts:
    X,V = P.xy(); pv = X/A
    opp = {c: 1 - c*pv for c in (a,b,a+b)}      # противоположные клетки 1 - c p
    own = {c: 1 + c*pv for c in (a,b,a+b)}
    print(f"p = {pv}:  клетки 1+cp = {list(own.values())},  противоположные 1-cp = {list(opp.values())}")
