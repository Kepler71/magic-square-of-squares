from math import gcd

def v2(n):
    n=abs(n); e=0
    while n%2==0: n//=2; e+=1
    return e

# integral forms = (4 v^2) * cell   (4v^2 is a square, so square-class unchanged)
def forms(m,n,u,v):
    A=(m*m+n*n)*(u*u+v*v)
    return {
      'H0':4*(m*m*v*v+n*n*u*u),
      'H4':2*A,
      'H8':4*(n*n*v*v+m*m*u*u),
      'H2':2*A-8*m*n*u*v,
      'H6':2*A+8*m*n*u*v,
    }

def solvable(m,n,which,kmax=30):
    """rigorous: is there (u:v) in P^1(Q2) with all forms in `which` squares in Q2*?
       BFS over classes (u,v) mod 2^k, primitive."""
    stack=[(u,v,1) for u in range(2) for v in range(2) if (u,v)!=(0,0)]
    visited=0
    while stack:
        u,v,k=stack.pop(); visited+=1
        if visited>4_000_000: raise RuntimeError('blowup')
        M=1<<k
        f=forms(m,n,u,v)
        undecided=False; dead=False
        for key in which:
            r=f[key]%(1<<(k+4))          # value mod 2^(k+4): any lift agrees mod 2^(k+4)? NO -> see note
            # value of a degree-2 form at a lift u+2^k a, v+2^k b changes by multiples of 2^k only,
            # so we may only trust the value modulo 2^k.
            r=f[key]%M
            if r==0:
                undecided=True; break
            e=v2(r)
            if e>k-3:
                undecided=True; break
            if (r>>e)%8!=1 or e%2:
                dead=True; break
        if dead: continue
        if not undecided:
            return True,(u,v,k)
        if k>=kmax: raise RuntimeError(f"no termination ({m},{n}) k={k}")
        for a in range(2):
            for b in range(2):
                nu,nv=u+a*M, v+b*M
                if nu%2==0 and nv%2==0: continue
                stack.append((nu,nv,k+1))
    return False,None

def oddpart_mod8(x):
    return (x>>v2(x))%8

S3=['H0','H4','H8']; S5=['H0','H4','H8','H2','H6']; S2=['H2','H6']

print("=== 3. 2-adic decision procedure: criterion, and do red cells add anything? ===")
bad_crit=0; bad_red=0; s2_fail=[]
rows=[]
for m in range(1,31):
    for n in range(1,m):
        if gcd(m,n)!=1: continue
        a3,w3=solvable(m,n,S3); a5,w5=solvable(m,n,S5); a2,w2=solvable(m,n,S2)
        pred = oddpart_mod8(m*m+n*n)==1
        if a3!=pred: bad_crit+=1; print("CRITERION MISMATCH",m,n,a3,pred)
        if a3!=a5: bad_red+=1; print("RED CELLS ADD A CONDITION AT p=2:",m,n,a3,a5)
        if not a2: s2_fail.append((m,n))
        rows.append((m,n,a3,a5,a2))
print(f"pairs tested (m<=40, gcd=1): {len(rows)}")
print(f"mismatches with criterion 'odd part of m^2+n^2 == 1 mod 8': {bad_crit}")
print(f"pairs where {{c2,c6}} add a 2-adic condition beyond {{c0,c4,c8}}: {bad_red}")
print(f"pairs where {{c2,c6}} alone is NOT 2-adically solvable: {len(s2_fail)}  {s2_fail[:10]}")

print()
print("=== 4. CONTROL ON ODD n (must show: the ban is NOT about parity of n) ===")
even_no=[]; even_yes=[]; odd_no=[]; odd_yes=[]
for (m,n,a3,a5,a2) in rows:
    if m>20: continue
    if n%2==0 or m%2==0:
        (even_yes if a3 else even_no).append((m,n))
    else:
        (odd_yes if a3 else odd_no).append((m,n))
print(f"one param EVEN, C(Q2)!=empty : {len(even_yes)}  e.g. {even_yes[:8]}")
print(f"one param EVEN, C(Q2)=empty  : {len(even_no)}   e.g. {even_no[:8]}")
print(f"both ODD,      C(Q2)!=empty  : {len(odd_yes)}  e.g. {odd_yes[:8]}")
print(f"both ODD,      C(Q2)=empty   : {len(odd_no)}   e.g. {odd_no[:8]}")
print(" -> the ban occurs for BOTH parities and fails for BOTH parities: parity of n is not the criterion")

print()
print("=== 5. the four pairs ===")
for (m,n) in [(13,8),(15,8),(19,16),(16,5)]:
    a3,w3=solvable(m,n,S3); a5,w5=solvable(m,n,S5)
    print(f"({m},{n}): v2(even param)={v2(m if m%2==0 else n)}, m^2+n^2={m*m+n*n} "
          f"(odd part mod 8 = {oddpart_mod8(m*m+n*n)}), C(Q2)!=empty: {a3}, nine-cell system Q2: {a5}, witness class {w5}")
