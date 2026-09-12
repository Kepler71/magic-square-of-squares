from math import gcd
def leg(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1

def cellsF(m,n,u,v,p):
    A=(m*m+n*n)%p
    inv2=pow(2,p-2,p)
    c0=(m*m*v*v+n*n*u*u)%p
    c8=(n*n*v*v+m*m*u*u)%p
    c4=(A*(u*u+v*v)*inv2)%p
    c2=(c4-2*m*n*u*v)%p
    c6=(c4+2*m*n*u*v)%p
    return c0,c2,c4,c6,c8

def has_Qp_point(m,n,p,which):
    """which: indices into (c0,c2,c4,c6,c8). Smooth-point test over F_p (Hensel).
       Returns (has_smooth_point, has_singular_branch)"""
    pts=[(u,1) for u in range(p)]+[(1,0)]
    sing=False
    for u,v in pts:
        cs=cellsF(m,n,u,v,p)
        vals=[cs[i] for i in which]
        if any(x==0 for x in vals): sing=True; continue
        if all(leg(x,p)==1 for x in vals): return True,sing
    return False,sing

print("=== 9. odd-prime local check (independent re-derivation of the side claim) ===")
print("S3 = seven-squares system {c0,c4,c8};  S5 = nine-cell system {c0,c2,c4,c6,c8}")
for (m,n) in [(13,8),(15,8),(19,16),(16,5)]:
    bad3=[];bad5=[];amb=[]
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113]:
        if p%2==0 or (m*m-n*n)%p==0 or (m*m+n*n)%p==0 or (m*n)%p==0: 
            continue
        h3,s3=has_Qp_point(m,n,p,[0,2,4])
        h5,s5=has_Qp_point(m,n,p,[0,1,2,3,4])
        if not h3: bad3.append((p,'sing' if s3 else 'clean'))
        if not h5: bad5.append((p,'sing' if s5 else 'clean'))
    print(f"({m},{n}): S3 no smooth F_p point at {bad3} | S5 no smooth F_p point at {bad5}")
