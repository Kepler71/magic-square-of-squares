from itertools import combinations
import json

def prime_factors(n):
    n=abs(n); out=[]; p=2
    while p*p<=n:
        if n%p==0:
            out.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:out.append(n)
    return out

def classes(n):
    ps=[p for p in prime_factors(n) if p%4==1]
    products=[1]
    for p in ps:products += [p*x for x in products]
    return sorted(d for d in products if d%24==1)

out={}
for r,s in [(126,451),(73,362),(265,298)]:
    out[f'{r}/{s}']={
        'minus_prime_support':prime_factors(r-s),
        'plus_prime_support':prime_factors(r+s),
        'T_classes':classes(r-s),
        'L_classes':classes(r+s)}
assert out['126/451']['T_classes']==[1]
assert out['126/451']['L_classes']==[1,577]
assert out['73/362']['T_classes']==[1]
assert out['73/362']['L_classes']==[1,145]
assert out['265/298']['T_classes']==[1]
assert out['265/298']['L_classes']==[1]
print(json.dumps(out,indent=2))
with open('class_lists.json','w') as f:json.dump(out,f,indent=2)
