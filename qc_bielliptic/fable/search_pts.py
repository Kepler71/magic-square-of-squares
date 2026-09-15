# поиск рациональных точек на C: w^2=(x^2-a^2)(x^2-b^2)(x^2-c^2), x=m/n, |m|,|n|<=B (точная арифметика)
import sys; from math import gcd, isqrt
a,b,c,B=map(int,sys.argv[1:5]); found=set()
for n in range(1,B+1):
    for m in range(0,B+1):
        if gcd(m,n)!=1: continue
        v=(m*m-a*a*n*n)*(m*m-b*b*n*n)*(m*m-c*c*n*n)*n**6
        if v<0: continue
        r=isqrt(v)
        if r*r==v: found.add((m,n))
print('точки x=m/n (с точностью до знака x):',sorted(found))
