# 48/163: второй независимый путь. Кривые из 6 клеток: две пары и две одиночные клетки разных значений
import itertools, json
from chab_all_lib import magma
r,s=48,163; vals=[r,s-r,s,s+r]
cases=[]
for pa in itertools.combinations(vals,2):
    rest=[v for v in vals if v not in pa]
    for sg in (1,-1):                       # знак второй одиночной; первая фиксирована +, т.к. z->-z
        cases.append(("6",[rest[0],sg*rest[1]],list(pa)))
# повтор 5-клеточных с таймаутом
for c,a,b in [(48,115,163),(48,115,211),(48,163,211),(115,48,211),(163,48,115),(163,48,211),(163,115,211),(211,48,115),(211,48,163),(211,115,163),(115,163,211)]:
    cases.append(("5",[c],[a,b]))
out=open('extra_48.jsonl','a')
for kind,singles,pair in cases:
    f="*".join(f"(1+({c})*z)" for c in singles)+"*"+"*".join(f"(1-{a}^2*z^2)" for a in pair)
    code=f"""P<z>:=PolynomialRing(Rationals());
C:=HyperellipticCurve({f});
pts,flag:=RationalPointsGenus2(C);
print "PROVEN:", flag; print "ZVALS:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in pts];"""
    res=magma(code); rec=dict(kind=kind,singles=singles,pair=pair,out=res)
    out.write(json.dumps(rec,ensure_ascii=False)+"\n"); out.flush()
    print(kind,singles,pair,res[:120].replace("\n"," | "),flush=True)
