# Три открытых наклона: кривые рода 2 из 5 клеток = одна клетка 1±cz и две пары (1−a²z²)(1−b²z²)
import itertools, json
from chab_all_lib import magma
out=open('rest3.jsonl','a')
for r,s in [(79,110),(11,142),(48,163)]:
    vals=[r,s-r,s,s+r]
    for single in vals:
        for sign in (1,):  # z -> -z переводит +c в -c при тех же парах
            rest=[v for v in vals if v!=single]
            for a,b in itertools.combinations(rest,2):
                code=f"""P<z>:=PolynomialRing(Rationals());
C:=HyperellipticCurve((1+({sign*single})*z)*(1-{a}^2*z^2)*(1-{b}^2*z^2));
pts,flag:=RationalPointsGenus2(C);
print "PROVEN:", flag; print "ZVALS:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in pts];"""
                res=magma(code)
                rec=dict(slope=f"{r}/{s}",single=sign*single,pair=[a,b],out=res)
                out.write(json.dumps(rec,ensure_ascii=False)+"\n"); out.flush()
                print(rec['slope'],sign*single,a,b,('PROVEN: true' in res) and 'ДОКАЗАНО' or res[:60].replace('\n',' '),flush=True)
