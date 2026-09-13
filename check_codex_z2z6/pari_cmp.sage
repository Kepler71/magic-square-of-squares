import json
res=json.load(open('rules_check.json'))['res']
bad=[]; n0=0; mine0=0
for key,r in res.items():
    a,b=map(int,key.split(',')); M=b^3*(2*a+b); N=a^3*(a+2*b)
    E=EllipticCurve([0,M+N,0,M*N,0]); lo,hi=E.pari_curve().ellrank()[:2]
    if hi==0: n0+=1
    bound=int(log(r['G_exact'],2))-2
    if r['G_exact']==4: mine0+=1
    if bound<lo: bad.append((a,b,bound,lo))
    if r['G_exact']==4 and hi!=0: bad.append((a,b,'mine0_pari_hi',hi))
print("PARI верхняя граница 0:",n0,"; мой покомпонентный точный Γ даёт 0:",mine0,"; противоречий:",bad)
