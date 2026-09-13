# Отождествление 14 хороших классов с E_{m,n} и их квадратичными кручениями, на нескольких k.
import json
T=json.load(open('full_table_2_48.json')); rows=T['rows']
GOODC=['1,k,1-k','1,k,1+k','1,k,-(1-k)','1,1-k,-k','1,-k,-(1+k)','1,1+k,-k']
GOODQ=['k,1-k,-k,-(1-k)','k,1+k,1-k,-k','1+k,1-k,-(1+k),-(1-k)','1,k,-1,-k','1,1-k,-1,-(1-k)','1,1+k,-1,-(1-k)','k,1+k,-k,-(1+k)','1,1+k,-1,-(1+k)']
def Emn(m,n): return EllipticCurve([0,m*m+n*n,0,m*m*n*n,0]).minimal_model()
for ks in ['5/7','3/11','7/12','11/13']:
    k=QQ(ks); r,s=k.numerator(),k.denominator()
    cand={'r':r,'s':s,'s-r':s-r,'r+s':r+s,'2r':2*r,'2s':2*s,'2(s-r)':2*(s-r),'2(r+s)':2*(r+s),'1':1}
    names=list(cand)
    print('== k =',ks)
    for kind,lst in (('C',GOODC),('Q',GOODQ)):
        for cls in lst:
            key=rows[ks][kind+'|'+cls]['key']; E=EllipticCurve([QQ(z) for z in key.split(',')])
            found=[]
            for i in range(len(names)):
                for j in range(i+1,len(names)):
                    a,b=cand[names[i]],cand[names[j]]
                    if gcd(a,b)!=1 or a==b: continue
                    F=Emn(a,b); D=E.is_quadratic_twist(F)
                    if D!=0: found.append((names[i],names[j],int(D)))
            print(f'{kind}|{cls:24s}', found)
