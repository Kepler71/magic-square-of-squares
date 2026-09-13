# Семь квартичных классов клеток: Y^2 = prod(1 + c_i p). Точка p=0 (Y=1) -> кубическая модель по формулам:
# r = -1/c1, x = 1/(p-r), w = Y/(p-r)^2, w^2 = g(x) = c1*prod_{i>=2}((1-c_i/c1)x + c_i), X = A x, V = A w.
# Перечисляем 8 точек кручения (О, три 2-точки, ±P0, ±(P0+T)) как сечения над Q(k), переводим в p, смотрим клетки.
K.<k> = FunctionField(QQ)
ALL = {'1+p':1,'1-p':-1,'1+kp':k,'1-kp':-k,'1+(1+k)p':1+k,'1-(1+k)p':-(1+k),'1+(1-k)p':1-k,'1-(1-k)p':-(1-k)}
classes = {
 'Якоби {±(1+k),±(1-k)}': [1+k,-(1+k),1-k,-(1-k)],
 '{±1,±k}':               [K(1),K(-1),k,-k],
 '{±k,±(1-k)}':           [k,-k,1-k,-(1-k)],
 '{±k,±(1+k)}':           [k,-k,1+k,-(1+k)],
 '{±1,±(1-k)}':           [K(1),K(-1),1-k,-(1-k)],
 '{±1,±(1+k)}':           [K(1),K(-1),1+k,-(1+k)],
 '{k,1+k,1-k,-k}':        [k,1+k,1-k,-k],
 '{1,1+k,-1,-(1-k)}':     [K(1),1+k,K(-1),-(1-k)],
}
Rx.<x> = K[]
for name, cs in classes.items():
    c1 = cs[0]; r = -1/c1
    g = c1*prod((1 - c/c1)*x + c for c in cs[1:])
    A,B,C,D = g[3],g[2],g[1],g[0]
    E = EllipticCurve(K,[0,B,0,A*C,A^2*D])
    x0 = 1/(0 - r); w0 = 1/(0 - r)^2
    P0 = E(A*x0, A*w0)
    roots = [ -c/(1 - c/c1) for c in cs[1:] ]           # корни g: (1-c/c1)x + c = 0
    T = [E(A*rt, 0) for rt in roots]
    print(f"\n=== {name}: порядок P0 = {P0.order(algorithm='generic_small')}")
    pts = [('P0',P0)] + [(f'T{i+1}',t) for i,t in enumerate(T)] + [(f'P0+T{i+1}',P0+t) for i,t in enumerate(T)]
    for lab,P in pts:
        if P == E(0): print(f"   {lab}: O"); continue
        X,V = P.xy(); xx = X/A
        if xx == 0: print(f"   {lab}: x=0 → p = ∞"); continue
        pv = r + 1/xx
        bad = [(nm, 1 + c*pv) for nm,c in ALL.items() if (1 + c*pv) in (K(2),K(-1),K(0),K(3),K(-2)) or ((1+c*pv).numerator().degree()==0 and (1+c*pv).denominator().degree()==0 and not QQ(1+c*pv).is_square())]
        consts = [(nm, v) for nm,v in bad]
        print(f"   {lab}: p = {pv};  клетки-константы не квадраты: {consts if consts else 'НЕТ — разобрать вручную'}")
