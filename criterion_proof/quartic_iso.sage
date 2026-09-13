# Символьная проверка: кривая квартичного класса изоморфна над Q(k) кривой E_{m,n}(k) из таблицы (k = r/s, масштаб по s).
R.<k> = QQ[]; K = R.fraction_field(); k = K(k)
def Emn(m,n): return EllipticCurve(K,[0,m^2+n^2,0,m^2*n^2,0])
def quart(cs):
    Rx.<x> = K[]; c1=cs[0]
    g = c1*prod((1 - c/c1)*x + c for c in cs[1:]); A,B,C,D = g[3],g[2],g[1],g[0]
    return EllipticCurve(K,[0,B,0,A*C,A^2*D])
def is_sq(f):
    f = K(f); num, den = f.numerator(), f.denominator()
    t = num*den                      # квадрат ⇔ num*den квадрат в Q[k] (с учётом константы)
    fa = t.factor()
    return fa.unit().is_square() and all(e % 2 == 0 for _,e in fa)
table = {   # класс -> (m,n) как функции k после деления на s: r->k, s->1
 'Якоби {±(1+k),±(1-k)}': ([1+k,-(1+k),1-k,-(1-k)], (k,1)),
 '{±1,±k}':               ([K(1),K(-1),k,-k],        (1-k,1+k)),
 '{±k,±(1-k)}':           ([k,-k,1-k,-(1-k)],        (1,2*k-1)),
 '{±k,±(1+k)}':           ([k,-k,1+k,-(1+k)],        (1,2*k+1)),
 '{±1,±(1-k)}':           ([K(1),K(-1),1-k,-(1-k)],  (k,2-k)),
 '{±1,±(1+k)}':           ([K(1),K(-1),1+k,-(1+k)],  (k,2+k)),
 '{k,1+k,1-k,-k}':        ([k,1+k,1-k,-k],           (1,2*k)),
 '{1,1+k,-1,-(1-k)}':     ([K(1),1+k,K(-1),-(1-k)],  (k,2)),
}
for name,(cs,(m,n)) in table.items():
    E1 = quart(cs); E2 = Emn(m,n)
    j_eq = (E1.j_invariant() == E2.j_invariant())
    u2 = (E2.c6()/E1.c6()) / (E2.c4()/E1.c4())
    print(f"{name:26} j совпадают: {j_eq};  u² = (c6'/c6)/(c4'/c4) — квадрат в Q(k): {is_sq(u2)}")
