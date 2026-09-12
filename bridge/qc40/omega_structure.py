from sage.all import *
from pathlib import Path
import json
out=Path(__file__).resolve().parent
src=(out/'vendor/qc_g2_bielliptic_sage109.sage').read_text()
old='    Wlist = list(itertools.product(*Wqprimelist))\n    Omega = []\n    for i in Wlist:\n        Omega.append(sum(list(i)))\n    Omega = f7(Omega)\n    return Omega'
assert src.count(old)==1
src=src.replace(old,'    return [(q, v) for q, v in zip(bad_primes, Wqprimelist)]')
exec(preparse(src),globals())
R=PolynomialRing(QQ,'x');x=R.gen();f=(49+529*x*x)*(83521*x**4+63358*x*x+83521)
rows=Omega_set(f,11,20)
data={'local_factors':[{'q':int(q),'count':len(v),'values':list(map(str,v))} for q,v in rows],'cartesian_count':int(prod(len(v) for q,v in rows)),'status':'standard overestimate before restrictions from the five square forms; no Cartesian expansion'}
(out/'omega_structure.json').write_text(json.dumps(data,indent=2)+'\n');print({r['q']:r['count'] for r in data['local_factors']});print(data['cartesian_count'])
