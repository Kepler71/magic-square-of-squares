"""Recheck Fisher certificate using Fraction/SymPy, without Sage/PARI Hilbert symbols."""
from fractions import Fraction as Q
from math import gcd,isqrt
from pathlib import Path
import json
import sympy as S
out=Path(__file__).resolve().parent
d=json.loads((out/'pairing_40_43.json').read_text())
gs=[list(map(Q,g)) for g in d['quartics']];I=Q(d['I']);J=Q(d['J'])
def invs(g):
 a,b,c,e,f=g
 return 12*a*f-3*b*e+c*c,72*a*c*f-27*a*e*e-27*b*b*f+9*b*c*e-2*c**3
def val(g,x,z):return sum(c*x**(4-i)*z**i for i,c in enumerate(g))
def valuation(x,p):
 n=x.numerator;den=x.denominator;v=0
 while n%p==0:n//=p;v+=1
 while den%p==0:den//=p;v-=1
 return v,Q(n,den)
def residue(x,m):return x.numerator*pow(x.denominator,-1,m)%m
def is_square_local(x,p):
 if not x:return True
 v,u=valuation(x,p)
 return v%2==0 and (residue(u,8)==1 if p==2 else pow(residue(u,p),(p-1)//2,p)==1)
def hilbert(a,b,place):
 if place=='real':return -1 if a<0 and b<0 else 1
 p=int(place);aa,u=valuation(a,p);bb,v=valuation(b,p)
 if p==2:
  u=residue(u,8);v=residue(v,8)
  exponent=((u-1)//2)*((v-1)//2)+aa*((v*v-1)//8)+bb*((u*u-1)//8)
 else:
  exponent=aa*bb*((p-1)//2)
  if pow(residue(u,p),(p-1)//2,p)!=1:exponent+=bb
  if pow(residue(v,p),(p-1)//2,p)!=1:exponent+=aa
 return -1 if exponent%2 else 1
def squareQ(q):return q>=0 and isqrt(q.numerator)**2==q.numerator and isqrt(q.denominator)**2==q.denominator
def primes(q):return set(S.factorint(abs(q.numerator)))|set(S.factorint(q.denominator))
assert all(invs(g)==(I,J) for g in gs)
phi=S.symbols('phi');f=phi**3-3*S.Rational(I)*phi+S.Rational(J)
rem=lambda z:S.rem(z,f,phi)
zeta=lambda g:rem((4*g[0]*phi+3*g[1]**2-8*g[0]*g[2])/3)
zs=list(map(zeta,gs));disc=(4*I**3-J**2)/27
checked=0
for run in d['runs']:
 i,j=(0,1) if run['direction']=='forward' else (1,0)
 a,b,c,e,h=gs[i]
 hess=[3*b*b-8*a*c,4*(b*c-6*a*e),2*(2*c*c-24*a*h-3*b*e)]
 H=[(4*phi*a+hess[0])/3,(4*phi*b+hess[1])/6,(4*phi*c+hess[2])/18+2*(I-phi**2)/9]
 m=sum(S.Rational(v)*phi**q for q,v in enumerate(run['m_coefficients']))
 assert rem(m*m-zs[0]*zs[1]*zs[2])==0
 fac=rem(m*S.invert(zs[i],f,phi))
 gamma=[Q(rem(fac*v).coeff(phi,2)) for v in H]
 assert gamma==list(map(Q,run['raw_gamma']))
 normalized=list(map(Q,run['gamma']))
 scale=next(n/g for n,g in zip(normalized,gamma) if g)
 assert scale!=0 and normalized==[scale*g for g in gamma]
 aa=Q(run['a'])
 assert squareQ(gs[j][0]/aa) or squareQ(gs[2][0]/aa)
 required={2,3,5,7}|primes(disc)|primes(aa)
 for coeff in gs[i]+normalized:required|=set(S.factorint(coeff.denominator))
 content=0
 for coeff in normalized:content=gcd(content,abs(coeff.numerator))
 required|=set(S.factorint(content))
 assert required <= {int(t) for t in run['places'] if t!='real'}
 total=1
 for loc in run['local']:
  x,z=Q(loc['x']),Q(loc['z']);gv=sum(normalized[q]*x**(2-q)*z**q for q in range(3))
  vv=val(gs[i],x,z)
  assert vv==Q(loc['quartic_value']) and gv==Q(loc['gamma_value']) and gv!=0 and vv!=0
  assert vv>0 if loc['place']=='real' else is_square_local(vv,int(loc['place']))
  hh=hilbert(aa,gv,loc['place']);assert hh==loc['hilbert'];total*=hh;checked+=1
 assert total==-1 and run['pairing_bit']==1
for section in d['everywhere_local_solubility']:
 g=gs[section['quartic_index']]
 required={2}|primes(disc)
 for c in g:required|=set(S.factorint(c.denominator))
 assert required <= {int(z['place']) for z in section['local_points'] if z['place']!='real'}
 for loc in section['local_points']:
  vv=val(g,Q(loc['x']),Q(loc['z']));assert vv==Q(loc['quartic_value']) and vv!=0
  assert vv>0 if loc['place']=='real' else is_square_local(vv,int(loc['place']))
# Exact Jacobian/original Weierstrass equation identity.
x,y=S.symbols('x y');u,r,s,t=map(S.Rational,d['jacobian_to_original_isomorphism'])
def equation(ai,x,y):
 a1,a2,a3,a4,a6=map(S.Rational,ai)
 return y*y+a1*x*y+a3*y-x**3-a2*x*x-a4*x-a6
assert S.expand(equation(d['quartic_jacobian_model'],u*u*x+r,u**3*y+s*u*u*x+t)-u**6*equation(d['curve_ainvs'],x,y))==0
report={'all_passed':True,'hilbert_symbols_recomputed':checked,'quartics_local_solubility_checked':3,
        'algebra_checked':'I,J; m^2=z1*z2*z3; gamma formula; proportional normalization; leading coefficient square class; complete bad-prime coverage; exact Jacobian isomorphism',
        'theoretical_inputs':'Fisher theorem and good-reduction local solubility; eclib Selmer dimension 4 (not independently recomputed here).',
        'pairing_rank_lower_bound':2,'rational_two_torsion_dimension':2,'rank_upper_bound':0}
(out/'pairing_40_43_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
