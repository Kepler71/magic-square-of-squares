"""Exact rational G1 lifting validator. Standard Python only."""
from fractions import Fraction as Q
from math import isqrt,lcm,gcd
from functools import reduce

def square_root(q):
 q=Q(q)
 if q<0:return None
 a,b=isqrt(q.numerator),isqrt(q.denominator)
 return Q(a,b) if a*a==q.numerator and b*b==q.denominator else None

def cells(m,n,t):
 t=Q(t);s=Q(m*m+n*n,2);c=s*(1+t*t)
 return [m*m+n*n*t*t,(m*t+n)**2,c-2*m*n*t,
         (m*t-n)**2,c,(m+n*t)**2,
         c+2*m*n*t,(m-n*t)**2,n*n+m*m*t*t]

def validate(m,n,t):
 t=Q(t);c=cells(m,n,t)
 lines=[c[0:3],c[3:6],c[6:9]]+[[c[i],c[i+3],c[i+6]] for i in range(3)]+[[c[0],c[4],c[8]],[c[2],c[4],c[6]]]
 sums=[sum(v) for v in lines];assert len(set(sums))==1 and sums[0]==3*c[4]
 roots=[square_root(v) for v in c];sq=all(v is not None for v in roots)
 rec={'m':m,'n':n,'t':str(t),'cells':list(map(str,c)),'nonsquare_indices':[i for i,v in enumerate(roots) if v is None],'all_squares':sq,'positive':all(v>0 for v in c),'distinct':len(set(c))==9,'all_eight_sums_equal':True,'valid_magic_square_of_distinct_positive_squares':sq and all(v>0 for v in c) and len(set(c))==9}
 if sq:
  mult=lcm(*(v.denominator for v in roots));ir=[int(v*mult) for v in roots]
  common=reduce(gcd,ir);ir=[v//common for v in ir]
  ic=[v*v for v in ir];assert all(Q(v)==c[i]*(Q(mult,common)**2) for i,v in enumerate(ic))
  rec['primitive_integer_roots']=ir;rec['primitive_integer_cells']=ic
 return rec

def quartic_data(m,n,gate):
 s=Q(m*m+n*n,2);d=Q((m*m-n*n)**2);h=Q(4*m*m*n*n)
 if gate=='minus':return s,2*s*s-h,[Q(0),d,-h]
 if gate=='plus':return Q(m*n),Q(m**4+n**4),[Q(0),d,d+h]
 raise ValueError(gate)

def to_elliptic(m,n,gate,t,w):
 t,w=Q(t),Q(w);a,c,er=quartic_data(m,n,gate)
 assert w*w==a*a*t**4+c*t*t+a*a
 if t==0:
  assert w in [a,-a]
  return None if w==a else (Q(0),Q(0))
 X=2*a*(w+a)/(t*t)+c;Y=2*a*X/t
 assert Y*Y==X*(X-er[1])*(X-er[2]);return X,Y

def from_elliptic(m,n,gate,P):
 a,c,er=quartic_data(m,n,gate)
 if P is None:return {'chart':'finite','t':Q(0),'w':a}
 X,Y=map(Q,P);assert Y*Y==X*(X-er[1])*(X-er[2])
 if not Y:
  if X==0:return {'chart':'finite','t':Q(0),'w':-a}
  assert X in er[1:];return {'chart':'infinity','w_over_t2':(X-c)/(2*a)}
 t=2*a*X/Y;w=(X-c)*t*t/(2*a)-a
 assert w*w==a*a*t**4+c*t*t+a*a
 assert to_elliptic(m,n,gate,t,w)==(X,Y)
 return {'chart':'finite','t':t,'w':w}
