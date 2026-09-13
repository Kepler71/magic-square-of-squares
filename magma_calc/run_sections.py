import json, subprocess, time, re, sys
secs=json.load(open('sections.json'))
out=open('sections_results.txt','a')
def ask(code):
    open('_q.m','w').write(code)
    r=subprocess.run(['curl','-s','-m','150','-A','Mozilla/5.0','--data-urlencode','input@_q.m',
                      'https://magma.maths.usyd.edu.au/xml/calculator.xml'],capture_output=True,text=True)
    body=re.sub(r'<[^>]*>','\n',r.stdout)
    lines=[l for l in body.split('\n') if l.strip() and not re.fullmatch(r'(60|50000|\d{9,}|2\.29-10|[\d.]+|[\d.]+MB)',l.strip())]
    return '\n'.join(lines)
for s in secs:
    b,h,n=s['bhn']; A,C=s['A'],s['C']
    for d3 in s['d3']:
        code=f"""P<t>:=PolynomialRing(Rationals());
H:=HyperellipticCurve({d3}*t*(t^2-1)*({A}*t-{C})*({C}*t-{A}));
H2,phi:=ReducedMinimalWeierstrassModel(H);
pts2,ok:=RationalPointsGenus2(H2);
pts:=[Inverse(phi)(P) : P in pts2];
print "NPTS:", #pts;
print "PTS:", pts;
print "COMPLETE:", ok;
"""
        t0=time.time(); res=ask(code); dt=time.time()-t0
        rec={'section':[b,h,n],'d3':d3,'seconds':round(dt,1),'output':res}
        out.write(json.dumps(rec,ensure_ascii=False)+'\n'); out.flush()
        print(f"({b},{h},{n}) d3={d3}: {dt:.0f} с | "+res.replace('\n',' | ')[:220], flush=True)
