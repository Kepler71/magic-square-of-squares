import json, time, re
from chab_all_lib import magma
out=open('confirm15.jsonl','a')
for line in open('confirm15.txt'):
    sl,c,a,b=line.split(); 
    code=f"""P<z>:=PolynomialRing(Rationals());
C:=HyperellipticCurve((1+{c}*z)*(1-{a}^2*z^2)*(1-{b}^2*z^2));
pts,flag:=RationalPointsGenus2(C);
print "PROVEN:", flag; print "Z:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in pts];
J:=Jacobian(C); pp:=Points(C : Bound:=1000); Q:=[q : q in [J!(R-S) : R in pp, S in pp | R ne S] | Order(q) eq 0];
S2,idx:=Chabauty(Q[1]); print "CHAB2:", idx;"""
    while True:
        t0=time.time(); res=magma(code); dt=time.time()-t0
        if 'too many connections' in res or res.strip()=='':
            print('лимит — пауза 1800 с',flush=True); time.sleep(1800); continue
        break
    out.write(json.dumps(dict(slope=sl,c=c,pair=[a,b],out=res[:500]),ensure_ascii=False)+"\n"); out.flush()
    print(sl,res[:160].replace('\n',' | '),flush=True); time.sleep(max(0,30-dt))
print('готово',flush=True)
