from sage.all import *
for lab in ['37a1','389a1']:
    E=EllipticCurve(lab); P=E.gens()[0]
    for n in [5,10,20]:
        Q=n*P; x=Q[0]; hx=log(max(abs(x.numerator()),abs(x.denominator())))
        print(lab,n,'ĥ_Sage',float(Q.height()),'h(x)',float(hx),'ratio',float(Q.height()/hx))
