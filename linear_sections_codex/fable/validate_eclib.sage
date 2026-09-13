# Печатаем маркеры и даём eclib писать в stdout (C++-вывод не перехватывается redirect_stdout); разбор — parse_eclib.py
import random,sys
from sage.libs.eclib.interface import mwrank_EllipticCurve
random.seed(int(1))
pairs=[(m,n) for n in range(2,150) for m in range(1,n) if gcd(m,n)==1]
for m,n in random.sample(pairs,150):
    M=EllipticCurve([0,m*m+n*n,0,m*m*n*n,0]).minimal_model()
    Ep=EllipticCurve([0,-2*(m*m+n*n),0,(m*m-n*n)**2,0]).minimal_model()
    print('###PAIR',m,n,'EPRIME',','.join(map(str,Ep.ainvs())),flush=True)
    sys.stdout.flush()
    ec=mwrank_EllipticCurve(list(map(int,M.ainvs())),verbose=True); ec.two_descent(verbose=True,second_descent=False)
    sys.stdout.flush()
print('###END')
