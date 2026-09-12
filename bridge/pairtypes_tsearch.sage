secs=[(17,7,13),(7,1,5),(23,7,17),(71,49,61)]
R.<T>=QQ[]
def forms(kind,b0,h0,n0):
    if kind=='edge':
        A=(h0^2+n0^2)/2; C=(b0^2+n0^2)/2; B=2*n0^2
    else:
        A=h0^2; C=b0^2; B=4*n0^2
    p=(T^2+1)^2; q=T*(T^2-1)
    return [A*p+B*q, A*p-B*q, C*p+B*q, C*p-B*q]
def cells9(kind,b0,h0,n0,t):
    bS=QQ(b0)/n0; hS=QQ(h0)/n0
    Bf=(1+2*t-t^2)/(1+t^2); Hf=(1-2*t-t^2)/(1+t^2)
    if kind=='edge':   # centre 1; pairs (b,h)=bS^2,hS^2 ; (d,f)=Bf^2,Hf^2 ; corners = means
        return [1,bS^2,hS^2,Bf^2,Hf^2,(hS^2+Hf^2)/2,(hS^2+Bf^2)/2,(bS^2+Hf^2)/2,(bS^2+Bf^2)/2]
    else:              # centre 1; corners (a,i)=bS^2,hS^2 ; (c,g)=Bf^2,Hf^2 ; edges = x+y-1
        return [1,bS^2,hS^2,Bf^2,Hf^2,bS^2+Bf^2-1,hS^2+Hf^2-1,hS^2+Bf^2-1,bS^2+Hf^2-1]
for kind in ('edge','corner'):
    for (b0,h0,n0) in secs:
        F=forms(kind,b0,h0,n0); hits=[]
        for qq in range(1,61):
            for pp in range(-60,61):
                if gcd(pp,qq)!=1: continue
                t=QQ(pp)/qq
                if all(QQ(f(t)).is_square() for f in F):
                    hits.append((t,len(set(cells9(kind,b0,h0,n0,t)))))
        print(f"{kind:7s} ({b0},{h0},{n0}): rational t with |num|,den<=60 making all four free cells squares -> (t,#distinct cells): {sorted(hits)}")
