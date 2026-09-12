from fractions import Fraction as Q
from pathlib import Path
import json
A0,A2,A4,A6=4092529,47287151,37608911,44182609
models=[(A4,A2*A6,A0*A6*A6),(A2,A4*A0,A6*A0*A0)]
Ns={2:(20,20),7:(2,6),17:(8,8),23:(6,2)}
def val(z,p):
 z=Q(z)
 if z==0:return 10**6
 n,d=z.numerator,z.denominator;v=0
 while n%p==0:n//=p;v+=1
 while d%p==0:d//=p;v-=1
 return v
def square(z,p):
 z=Q(z)
 if z==0:return True
 v=val(z,p)
 if v%2:return False
 u=z/Q(p)**v
 if p==2:return u.numerator*pow(u.denominator,-1,8)%8==1
 return pow(u.numerator*pow(u.denominator,-1,p)%p,(p-1)//2,p)==1

def lam(i,x,vy,p):
 a,b,c=models[i];aa=val(3*x*x+2*a*x+b,p);bb=val(2,p)+vy
 if aa<=0 or bb<=0:return Q(max(0,-val(x,p)))
 if p!=2:
  N=Ns[p][i];n=min(Q(bb),Q(N,2));return -n*(N-n)/N
 cc=val(3*x**4+4*a*x**3+6*b*x*x+12*c*x+4*a*c-b*b,p)
 return -Q(2*bb,3) if cc>=3*bb else -Q(cc,4)
rows=[]
for p in [2,7,17,23]:
 seen={};count=0;byk={}
 for k in range(-4,5):
  for u in range(1,p*p+1):
   if u%p==0:continue
   x=Q(u)*Q(p)**k;center=289*(1+x*x)
   forms=[center-322*x,center,center+322*x,49+529*x*x,529+49*x*x]
   if not all(square(z,p) for z in forms):continue
   f=(49+529*x*x)*(83521*x**4+63358*x*x+83521);vf=val(f,p);assert vf%2==0
   vy=vf//2
   r1=lam(0,A6*x*x,val(A6,p)+vy,p)
   r2=lam(1,Q(A0)/(x*x),val(A0,p)+vy-3*k,p)
   w=-r1+r2-2*k
   expect=Q(0) if p==2 else (Q(-2) if k>0 else Q(-4,3)) if p==7 else (Q(2) if k<0 else Q(4,3)) if p==23 else (Q(0) if k==0 else Q(-2) if k>0 else Q(2))
   assert w==expect,(p,x,k,w,expect)
   assert p!=17 or abs(k)!=1
   count+=1;seen.setdefault(str(w),str(x));byk.setdefault(str(k),set()).add(str(w))
 rows.append({'p':p,'samples_passing_all_five':count,'w_coefficients_witnesses':seen,'by_vx':{k:sorted(v) for k,v in byk.items()}})
 print(rows[-1],flush=True)
(Path(__file__).parent/'height_coefficients_check.json').write_text(json.dumps(rows,indent=2)+'\n')
