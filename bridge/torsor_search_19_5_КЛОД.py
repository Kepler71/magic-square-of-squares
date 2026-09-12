from math import isqrt, gcd
# torsor D for class (1,193,193):
#   z1^2 - 193 z2^2 = (e2-e1) w^2 = -23410128 w^2
#   z1^2 - 193 z3^2 = (e3-e1) w^2 =  1621200 w^2
# 193 | z1  =>  z1 = 193*a :
#   A = 193 a^2 + 121296 w^2 = z2^2
#   B = 193 a^2 -   8400 w^2 = z3^2
def issq(n):
    if n < 0: return False
    r = isqrt(n)
    return r*r == n
found=[]
H=200000
for w in range(1, 4001):
    w2 = w*w
    c1 = 121296*w2
    c2 = 8400*w2
    for a in range(0, 4001):
        if gcd(a,w)!=1: continue
        t = 193*a*a
        B = t - c2
        if B < 0: continue
        if not issq(B): continue
        if issq(t + c1):
            found.append((a,w))
            print("HIT", a, w)
print("searched a,w <= 4000 ; hits:", found)
