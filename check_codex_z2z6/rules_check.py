"""Claude, 13.09: независимая проверка THEOREMS_Z2Z6 §2 (отбрасывания 1–3) и границы через Γ.
Свой точный решатель над Q_p (criterion_proof/selmer_exact.py); код Codex не используется.
E(a,b): y²=x(x+M)(x+N), M=b³(2a+b), N=a³(a+2b); корни e=(0,−M,−N)."""
import sys, json, random
from math import gcd, isqrt
from itertools import combinations
sys.path.insert(0,'/home/kep/magicKube/criterion_proof')
from selmer_exact import locally_solvable, primes_of, DepthExhausted
def legendre(x,q):
    x%=q
    return 0 if x==0 else (1 if pow(x,(q-1)//2,q)==1 else -1)
def v(x,p):
    e=0
    while x%p==0: x//=p; e+=1
    return e
def comp_data(M,N,i):
    e=[0,-M,-N]; j,k=[t for t in range(3) if t!=i]
    h=e[i]-e[j]; kk=e[i]-e[k]; return h,kk,h+kk,h*kk
def cands(B):
    P=primes_of(B); out=[]
    for r in range(len(P)+1):
        for sub in combinations(P,r):
            d=1
            for q in sub: d*=q
            out+= [d,-d]
    return out
def real_ok(d,A,B):
    # W² = dU⁴+AU²V²+(B/d)V⁴ имеет вещественную точку?
    c4,c2,c0=d,A,B//d
    if c4>0 or c0>0: return True
    # оба ≤0: нужно, чтобы квадратичная по t=U²/V² форма была >0 где-то при t≥0: max при t*=-c2/(2c4)
    return c2>0 and c2*c2>4*c4*c0
def codex_reject(d,h,k,A,B):
    """Возвращает список (правило, место) — какие отбрасывания Codex срабатывают."""
    rej=[]
    if A<0 and B>0 and d<0: rej.append((1,'inf'))
    for q in primes_of(abs(h-k)):
        if q==2 or B%q==0: continue
        if legendre(-h,q)==1 and legendre(d,q)==-1: rej.append((2,q))
    for p in primes_of(abs(B)):
        if p==2 or A%p==0: continue
        if v(abs(B),p)%2==0 and legendre(A,p)==-1 and d%p==0: rej.append((3,p))
    return rej
def quartic(d,A,B): return [B//d,0,A,0,d]
random.seed(20260913)
pairs=[(a,b) for b in range(2,41) for a in range(1,b) if gcd(a,b)==1]
extra=[]
while len(extra)<80:
    a,b=random.randint(1,400),random.randint(1,400)
    if a<b and gcd(a,b)==1: extra.append((a,b))
viol=[]; fired=0; checked=0; exact_rank0=0; res={}
for (a,b) in pairs+extra:
    M=b**3*(2*a+b); N=a**3*(a+2*b)
    bad=sorted(set([2,3]+primes_of(a*b*(b-a)*(a+b)*(a+2*b)*(2*a+b))))
    S_exact=[]; S_codex=[]
    for i in range(3):
        h,k,A,B=comp_data(M,N,i)
        Se=[]; Sc=[]
        for d in cands(B):
            rej=codex_reject(d,h,k,A,B)
            Q=quartic(d,A,B)
            for rule,place in rej:
                fired+=1
                ok = real_ok(d,A,B) if place=='inf' else locally_solvable(Q,place)
                if ok: viol.append((a,b,i,d,rule,place))
            if not rej: Sc.append(d)
            if (a,b) in pairs:
                try:
                    if real_ok(d,A,B) and all(locally_solvable(Q,p) for p in bad): Se.append(d)
                except DepthExhausted:
                    Se.append(d)   # не решено — оставляем (осторожно)
        S_exact.append(Se); S_codex.append(Sc)
    checked+=1
    if (a,b) in pairs:
        def issq(x): return x>0 and isqrt(x)**2==x
        G=[(x,y,z) for x in S_exact[0] for y in S_exact[1] for z in S_exact[2] if issq(x*y*z)]
        Gc=[(x,y,z) for x in S_codex[0] for y in S_codex[1] for z in S_codex[2] if issq(x*y*z)]
        # содержится ли точный набор в наборе Codex (отбрасывания необходимы)?
        if not set(G)<=set(Gc): viol.append((a,b,'Gamma_exact_not_subset'))
        res[f"{a},{b}"]=dict(G_exact=len(G),G_codex_odd=len(Gc))
        if len(G)==4: exact_rank0+=1
print("пар:",checked,"срабатываний отбрасываний:",fired,"нарушений:",len(viol))
print("по сетке 489: покомпонентный точный локальный Γ даёт ранг 0 у",exact_rank0,
      "; только нечётные правила Codex:",sum(1 for r in res.values() if r['G_codex_odd']==4))
json.dump(dict(violations=viol,res=res),open('rules_check.json','w'))
