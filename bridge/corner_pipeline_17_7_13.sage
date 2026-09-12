# CORNER analogue of family/pipeline.sage steps (1)-(3).
# Only difference from the edge case: A = h0^2, C = b0^2  (edge case: A=(h0^2+n0^2)/2, C=(b0^2+n0^2)/2).
# Everything below is copied structurally from family/pipeline.sage, which uses ONLY A and C in steps (1)-(3).
import functools, itertools, sys, time
print = functools.partial(print, flush=True)
b0, h0, n0 = ZZ(17), ZZ(7), ZZ(13)
assert b0^2 + h0^2 == 2*n0^2 and gcd([b0,h0,n0])==1
A, C = h0^2, b0^2                    # <<< CORNER-fixed constants
print(f"CORNER section (b,h,n)=({b0},{h0},{n0}):  A={A}, C={C}, beta=4n^2={4*n0^2} (NOT A+C={A+C})")
assert gcd(A,C)==1
Fh=[lambda p,q: p^2-q^2, lambda p,q: p^2+q^2, lambda p,q: A*p^2-C*q^2, lambda p,q: C*p^2-A*q^2]
def sqf(n):
    n=ZZ(n); return sign(n)*prod(pr^(e%2) for pr,e in n.abs().factor())
def is_sq_Ql(a,l):
    a=QQ(a); v=a.valuation(l)
    if v%2: return False
    u=a/l^v
    if l==2: return (u.numerator()*u.denominator())%8==1
    return kronecker(ZZ(u.numerator()*u.denominator()),l)==1
def locally_solvable(ds,l,K=None):
    e=3 if l==2 else 1
    if K is None: K=16 if l==2 else 10
    unknown=[]
    def rec(chart,a,k):
        p,q=(a,1) if chart==0 else (1,a)
        und=False
        for fh,d in zip(Fh,ds):
            v=fh(p,q)
            if v!=0 and v.valuation(l)+e<=k:
                if not is_sq_Ql(v/d,l): return False
            else: und=True
        if not und: return True
        if k>=K: unknown.append((chart,a,k)); return False
        return any(rec(chart,a+j*l^k,k+1) for j in range(l))
    if rec(0,0,0) or rec(1,0,1): return True
    return False if not unknown else None
def real_ok(ds):
    xs=[QQ(x)/8 for x in range(-80,81)]+[QQ(10^6)]
    return any(all(fh(x,1)==0 or sign(fh(x,1))==sign(d) for fh,d in zip(Fh,ds)) for x in xs)
res_primes=set((2*(C-A)*(C+A)*A*C).prime_factors()); cand=sorted(res_primes)
units=[s_*prod(pr for pr,e_ in zip(cand,es) if e_) for s_ in (1,-1) for es in itertools.product((0,1),repeat=len(cand))]
D2=[d for d in units if d>0]
check_primes=sorted(set(cand)|set(primes(3,60)))
t0=time.time(); els=[]
for d1 in units:
    for d2 in D2:
        d3=sqf(d1*d2); ds=(d1,d2,d3,d3)
        if not real_ok(ds): continue
        bad=None
        for l in check_primes:
            r=locally_solvable(ds,l)
            if r is False: bad=l; break
            if r is None: print(f"   UNKNOWN local test {ds[:3]} at {l}")
        if bad is None: els.append(ds[:3])
print(f"(1) ELS classes ({len(units)*len(D2)} candidates, {time.time()-t0:.0f}s): {els}")
H=300
pts={cls:[] for cls in els}
for q in range(1,H):
    for p in range(0,H):
        if gcd(p,q)!=1: continue
        vals=[fh(p,q) for fh in Fh]
        if 0 in vals: continue
        cls=(sqf(vals[0]),sqf(vals[1]),sqf(vals[2]))
        if cls in pts and sqf(vals[3])==cls[2] and all(QQ(v/d).is_square() for v,d in zip(vals,cls+(cls[2],))):
            pts[cls].append(QQ(p)/q)
print(f"(2) points with max(p,q)<{H}:", {c_:v for c_,v in pts.items()})
Rt=PolynomialRing(QQ,'T'); T=Rt.gen()
f=T*(T^2-1)*(A*T-C)*(C*T-A)
Rz=PolynomialRing(QQ,'z'); z=Rz.gen(); Frz=Rz.fraction_field()
Fz=Rz(Frz((1-z)^6)*Frz(f((1+z)/(1-z))))
cz=Fz.leading_coefficient(); q4=Rz(Fz/(cz*z))
Ap,Bp=q4[2],q4[0]
bet=(C-A)/(C+A); dk=sqf((C-A)*(C+A))
print(f"    identity check Bp==bet^2 and Ap==-(1+bet^2): {Bp==bet^2 and Ap==-(1+bet^2)}   bet={bet}")
k=QuadraticField(dk,'rk'); sb=k(bet).sqrt()
Ru=PolynomialRing(k,'uu'); uu=Ru.gen()
print(f"(3) k = Q(sqrt({dk})),  beta'={bet}")
for cls in els:
    D=cls[2]
    cub=Ru(D*cz*(uu^2+Ap-2*bet)*(uu+2*sb))
    a3,a2,a1,a0=[cub[j] for j in (3,2,1,0)]
    E1=EllipticCurve(k,[0,a2,0,a1*a3,a0*a3^2])
    Em=E1.global_minimal_model(semi_global=True)
    print(f"   class {cls}: E1={Em.ainvs()}  j={Em.j_invariant()}  torsion={Em.torsion_subgroup().invariants()}")
