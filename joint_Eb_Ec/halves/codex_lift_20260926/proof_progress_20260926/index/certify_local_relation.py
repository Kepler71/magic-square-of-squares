import json
from sage.all import *
S=PolynomialRing(QQ,'c');c=S.gen();K=S.fraction_field();b=QQ(24)
E=EllipticCurve(QQ,[-b*b,0])
W=PolynomialRing(QQ,'x');x=W.gen()
x3=E.multiplication_by_m(3,x_only=True)
X=K(1-c)
Y=K((1+b*b)**2/(4*(1-b*b)))
Z=K(x3.numerator()(1+c)/x3.denominator()(1+c))
f=lambda u:u**3-b*b*u
d=(Y-Z)**2
Phi=(d*(X+Y+Z)-f(Y)-f(Z))**2-4*f(Y)*f(Z)
num=Phi.numerator()
fac=[(q,int(e)) for q,e in num.factor()]
q=[q for q,e in fac if q.degree()==19][0]
q=S(q*q.denominator());q=S(q/gcd([ZZ(cf) for cf in q.coefficients()]))
polys={'Phi_denominator':Phi.denominator(),'x3_denominator':Z.denominator(), 'x2_x3_difference':(Y-Z).numerator(), 'cell_collision':c*(c-b)*(c+b)*(c-2*b)*(c+2*b)*(2*c-b)*(2*c+b)}
assert all(q.gcd(S(v)).degree()==0 for v in polys.values())
out={'b':str(b),'factor_degrees':[(int(t.degree()),e) for t,e in fac], 'q_coefficients':[str(t) for t in q.list()], 'exact_gcds':{name:str(q.gcd(S(v))) for name,v in polys.items()},'places':[]}
for p,c0 in [(2,8),(3,6)]:
 vf=q(c0).valuation(p);vd=q.derivative()(c0).valuation(p)
 assert vf>2*vd
 F=Qp(p,100);cc=F(c0)
 for it in range(8):cc-=q(cc)/q.derivative()(cc)
 EE=EllipticCurve(F,[-b*b,0]);Ps=[]
 for j in [-1,0,1]:
  xx=1+j*cc;yy=f(xx).sqrt();Ps.append(EE(xx,yy))
 vals=[1+i*b+j*cc for i in [-1,0,1] for j in [-1,0,1]]
 assert all(v.is_square() for v in vals)
 near=[]
 for s in [-1,1]:
  for t in [-1,1]:
   T=Ps[0]+s*2*Ps[1]+t*3*Ps[2]
   near.append({'coeffs':[1,s*2,t*3], 'is_zero_to_precision':bool(T.is_zero()),'x_valuation':None if T.is_zero() else str(T[0].valuation())})
 out['places'].append({'p':p,'c0':c0,'v_q':int(vf),'v_derivative':int(vd),'unique_root_modulus':p**int(vf-vd),'root_approx':str(cc),'all_nine_locally_squares':True,'relation_approximations':near})
print(json.dumps(out,indent=2));open('local_relation_certificate.json','w').write(json.dumps(out,indent=2))
