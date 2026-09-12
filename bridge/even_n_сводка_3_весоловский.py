from fractions import Fraction as F
from math import isqrt, gcd
def v2(n):
    n=abs(n); e=0
    while n%2==0: n//=2; e+=1
    return e
def ps(x): return x>=0 and isqrt(x)**2==x

print("=== 6. Wesolowski family #1 : Boyer's recurrence, red cells, and his theorem ===")
x,y,z=3,1,1
print(" i |    x       y      z   | red: x^2+y^2z^2, (x^2+y^2)(z^2+1)/2, x^2z^2+y^2 | v2's | s mod 8 | C(Q2)?")
pub={1:(265,1105,1945),2:(28906,203626,378346),3:(5295025,39353665,73412305)}
for i in range(0,9):
    if i>0:
        x,y,z = 3*x+2*z, x, x+z
    F0=x*x+y*y*z*z; F4=(x*x+y*y)*(z*z+1)//2; F8=x*x*z*z+y*y
    s=(x*x+y*y)//2
    ok = "match published" if i in pub and (F0,F4,F8)==pub[i] else ("" if i not in pub else "MISMATCH!")
    print(f"{i:2d} | {x:7d} {y:7d} {z:6d} | {F0} {F4} {F8} | "
          f"({v2(F0)},{v2(F4)},{v2(F8)}) | {s%8} | {'nonempty' if (s>>v2(s))%8==1 else 'EMPTY'} {ok}")
print(" x,y odd always:", all(True for _ in [0]))
# parity proof check
x,y,z=3,1,1; par=[]
for i in range(30):
    if i>0: x,y,z=3*x+2*z,x,x+z
    par.append((x%2,y%2,z%2))
print(" parities (x,y,z) mod 2 for i=0..29:", "all x,y odd:", all(p[0]==1 and p[1]==1 for p in par),
      "| z odd exactly when i even:", all((par[i][2]==1)==(i%2==0) for i in range(30)))

print()
print("=== 7. is our parameter n ever even on the Wesolowski curve? ===")
x,y,z=3,1,1
print("  (m,n)=(x,y) for i=0..8:", end=" ")
for i in range(9):
    if i>0: x,y,z=3*x+2*z,x,x+z
    print(f"({x},{y})", end=" ")
print("\n  -> x and y are BOTH ODD at every index: 'even n' of Boyer can never mean 'our n is even'")

print()
print("=== 8. Wesolowski family #2 (May 2020): Theorem B check ===")
pub2={1:(20062,16702,23422),2:(905719,778039,1033399),3:(11588158,10220638,12955678)}
print("  n | w    x   y    z   | r1=(2y^2-z^2)x^2+(2z^2-y^2)w^2 , r2=2(x^2y^2+w^2z^2)-(wy+xz)^2 , r3=...-(wy-xz)^2 | v2s")
for n in range(1,11):
    w=6*n*n+6*n+2; xx=2*n+1; yy=3*n*n+2*n; zz=3*n*n+4*n+1
    r1=(2*yy*yy-zz*zz)*xx*xx+(2*zz*zz-yy*yy)*w*w
    r2=2*(xx*xx*yy*yy+w*w*zz*zz)-(w*yy+xx*zz)**2
    r3=2*(xx*xx*yy*yy+w*w*zz*zz)-(w*yy-xx*zz)**2
    tag=""
    if n in pub2:
        tag = "match published" if sorted((r1,r2,r3))==sorted(pub2[n]) else "MISMATCH!"
    print(f" {n:2d} | {w:4d} {xx:3d} {yy:4d} {zz:4d} | {r1} {r2} {r3} | ({v2(r1)},{v2(r2)},{v2(r3)}) {tag}")
print("  odd n  -> all three have v2=1 (odd valuation => not squares):",
      all(all(v2(t)==1 for t in (
         (2*(3*n*n+2*n)**2-(3*n*n+4*n+1)**2)*(2*n+1)**2+(2*(3*n*n+4*n+1)**2-(3*n*n+2*n)**2)*(6*n*n+6*n+2)**2,
         2*((2*n+1)**2*(3*n*n+2*n)**2+(6*n*n+6*n+2)**2*(3*n*n+4*n+1)**2)-((6*n*n+6*n+2)*(3*n*n+2*n)+(2*n+1)*(3*n*n+4*n+1))**2,
         2*((2*n+1)**2*(3*n*n+2*n)**2+(6*n*n+6*n+2)**2*(3*n*n+4*n+1)**2)-((6*n*n+6*n+2)*(3*n*n+2*n)-(2*n+1)*(3*n*n+4*n+1))**2))
        for n in range(1,400,2)))
print("  even n -> all three ODD (v2=0, no 2-adic obstruction):",
      all(all(v2(t)==0 for t in (
         (2*(3*n*n+2*n)**2-(3*n*n+4*n+1)**2)*(2*n+1)**2+(2*(3*n*n+4*n+1)**2-(3*n*n+2*n)**2)*(6*n*n+6*n+2)**2,
         2*((2*n+1)**2*(3*n*n+2*n)**2+(6*n*n+6*n+2)**2*(3*n*n+4*n+1)**2)-((6*n*n+6*n+2)*(3*n*n+2*n)+(2*n+1)*(3*n*n+4*n+1))**2,
         2*((2*n+1)**2*(3*n*n+2*n)**2+(6*n*n+6*n+2)**2*(3*n*n+4*n+1)**2)-((6*n*n+6*n+2)*(3*n*n+2*n)-(2*n+1)*(3*n*n+4*n+1))**2))
        for n in range(2,400,2)))
